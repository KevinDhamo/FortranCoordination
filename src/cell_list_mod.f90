module cell_list_mod
    use types_mod
    use config_mod
    use error_mod
    use omp_lib
    implicit none
    private

    ! Module variables
    type(cell_type), allocatable :: cells(:,:,:)
    integer :: n_cells_x, n_cells_y, n_cells_z
    real(dp) :: cell_size
    real(dp) :: max_cutoff
    integer :: last_cell_update = -1
    real(dp), allocatable :: prev_coords(:,:)
    integer :: num_cell_resets = 0
    integer :: num_cell_updates = 0
    
    ! Cell statistics variables
    integer :: total_empty_cells = 0
    integer :: max_atoms_in_cell = 0
    integer :: min_atoms_in_cell = 0
    real(dp) :: avg_atoms_per_cell = 0.0_dp
    integer :: stats_collection_count = 0

    ! Public procedures and variables
    public :: initialize_cell_list
    public :: update_cell_list
    public :: get_neighboring_cells
    public :: cleanup_cell_list
    public :: get_cell_grid_dims
    public :: get_cell
    public :: get_cell_statistics
    public :: num_cell_resets, num_cell_updates
    public :: total_empty_cells, max_atoms_in_cell, min_atoms_in_cell, avg_atoms_per_cell

contains

    subroutine initialize_cell_list(box_length, cutoffs, n_atoms)
        real(dp), intent(in) :: box_length(3)
        real(dp), intent(in) :: cutoffs(:)
        integer, intent(in) :: n_atoms
        integer :: ix, iy, iz, alloc_stat
        
        ! Find maximum cutoff distance
        max_cutoff = maxval(cutoffs)
        
        ! Calculate cell size and grid dimensions
        cell_size = max_cutoff * CELL_SIZE_FACTOR
        n_cells_x = max(1, floor(box_length(1)/cell_size))
        n_cells_y = max(1, floor(box_length(2)/cell_size))
        n_cells_z = max(1, floor(box_length(3)/cell_size))
        
        ! Adjust cell size to fit box exactly
        cell_size = min(box_length(1)/n_cells_x, &
                       box_length(2)/n_cells_y, &
                       box_length(3)/n_cells_z)

        ! Allocate cell grid
        if (allocated(cells)) deallocate(cells)
        allocate(cells(n_cells_x, n_cells_y, n_cells_z), stat=alloc_stat)
        call check_allocation(alloc_stat, "cell grid")

        ! Initialize cells
        !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(ix,iy,iz)
        do iz = 1, n_cells_z
            do iy = 1, n_cells_y
                do ix = 1, n_cells_x
                    call cells(ix,iy,iz)%init(MAX_ATOMS_PER_CELL)
                end do
            end do
        end do
        !$OMP END PARALLEL DO

        ! Allocate previous coordinates array
        if (allocated(prev_coords)) deallocate(prev_coords)
        allocate(prev_coords(n_atoms,3), stat=alloc_stat)
        call check_allocation(alloc_stat, "previous coordinates")
        
        if (VERBOSE) then
            write(*,'(A)') ' Cell list initialization:'
            write(*,'(A,F10.4)') '   Cell size: ', cell_size
            write(*,'(A,3I6)') '   Grid dimensions:', n_cells_x, n_cells_y, n_cells_z
            write(*,'(A,I0)') '   Total cells: ', n_cells_x * n_cells_y * n_cells_z
        end if
    end subroutine initialize_cell_list

    subroutine update_cell_list(coords, atom_types, n_atoms, box_length, frame, include_mask)
        real(dp), intent(in) :: coords(:,:)
        integer, intent(in) :: atom_types(:)
        integer, intent(in) :: n_atoms
        real(dp), intent(in) :: box_length(3)
        integer, intent(in) :: frame
        logical, intent(in) :: include_mask(:)
        
        integer :: ix, iy, iz, i
        real(dp) :: max_displacement, dx, dy, dz, disp_sq
        logical :: force_update
        
        ! Check if update is needed
        force_update = (last_cell_update == -1)
        
        if (.not. force_update .and. mod(frame, CELL_UPDATE_FREQ) /= 0) then
            ! Check maximum atomic displacement
            max_displacement = 0.0_dp
            !$OMP PARALLEL DO PRIVATE(dx,dy,dz,disp_sq) REDUCTION(max:max_displacement)
            do i = 1, n_atoms
                if (.not. include_mask(atom_types(i))) cycle
                
                dx = coords(i,1) - prev_coords(i,1)
                dy = coords(i,2) - prev_coords(i,2)
                dz = coords(i,3) - prev_coords(i,3)
                
                ! Apply minimum image convention
                dx = dx - nint(dx/box_length(1)) * box_length(1)
                dy = dy - nint(dy/box_length(2)) * box_length(2)
                dz = dz - nint(dz/box_length(3)) * box_length(3)
                
                disp_sq = dx*dx + dy*dy + dz*dz
                max_displacement = max(max_displacement, sqrt(disp_sq))
            end do
            !$OMP END PARALLEL DO
            
            force_update = (max_displacement > REBUILD_THRESHOLD * cell_size)
        end if
        
        if (.not. force_update) then
            num_cell_updates = num_cell_updates + 1
            return
        end if
        
        ! Reset all cells
        !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(ix,iy,iz)
        do iz = 1, n_cells_z
            do iy = 1, n_cells_y
                do ix = 1, n_cells_x
                    cells(ix,iy,iz)%n_atoms = 0
                end do
            end do
        end do
        !$OMP END PARALLEL DO
        num_cell_resets = num_cell_resets + 1
        
        ! Assign atoms to cells
        !$OMP PARALLEL DO PRIVATE(i) SCHEDULE(dynamic)
        do i = 1, n_atoms
            if (.not. include_mask(atom_types(i))) cycle
            call add_atom_to_cell(i, coords(i,:), box_length)
        end do
        !$OMP END PARALLEL DO
        
        ! Update previous coordinates
        prev_coords = coords
        last_cell_update = frame
        
        ! Collect statistics on cell usage
        call collect_cell_statistics()
    end subroutine update_cell_list

    subroutine add_atom_to_cell(atom_idx, pos, box_length)
        integer, intent(in) :: atom_idx
        real(dp), intent(in) :: pos(3)
        real(dp), intent(in) :: box_length(3)
        integer :: ix, iy, iz
        
        ! Calculate cell indices
        ix = floor(pos(1)/cell_size) + 1
        iy = floor(pos(2)/cell_size) + 1
        iz = floor(pos(3)/cell_size) + 1
        
        ! Apply periodic boundary conditions
        ix = modulo(ix-1, n_cells_x) + 1
        iy = modulo(iy-1, n_cells_y) + 1
        iz = modulo(iz-1, n_cells_z) + 1
        
        ! Resize cell if needed
        !$OMP CRITICAL
        if (cells(ix,iy,iz)%n_atoms >= size(cells(ix,iy,iz)%atoms)) then
            call cells(ix,iy,iz)%resize(2 * size(cells(ix,iy,iz)%atoms))
        end if
        
        ! Add atom to cell
        cells(ix,iy,iz)%n_atoms = cells(ix,iy,iz)%n_atoms + 1
        cells(ix,iy,iz)%atoms(cells(ix,iy,iz)%n_atoms) = atom_idx
        !$OMP END CRITICAL
    end subroutine add_atom_to_cell

    subroutine get_neighboring_cells(ix, iy, iz, neighbor_cells, n_neighbors)
        integer, intent(in) :: ix, iy, iz
        integer, intent(out) :: neighbor_cells(27,3)
        integer, intent(out) :: n_neighbors
        integer :: dx, dy, dz, nx, ny, nz
        
        n_neighbors = 0
        do dz = 0, 1
            do dy = -1, 1
                do dx = -1, 1
                    ! Skip cells that are not touching the current cell
                    if (dx == -1 .and. dy == -1 .and. dz == 0) cycle
                    if (dx == -1 .and. dy == 1 .and. dz == 0) cycle
                    if (dx == -1 .and. dz == 1) cycle
                    
                    ! Apply periodic boundary conditions
                    nx = modulo(ix+dx-1, n_cells_x) + 1
                    ny = modulo(iy+dy-1, n_cells_y) + 1
                    nz = modulo(iz+dz-1, n_cells_z) + 1
                    
                    n_neighbors = n_neighbors + 1
                    neighbor_cells(n_neighbors,1) = nx
                    neighbor_cells(n_neighbors,2) = ny
                    neighbor_cells(n_neighbors,3) = nz
                end do
            end do
        end do
    end subroutine get_neighboring_cells

    subroutine get_cell_grid_dims(nx, ny, nz)
        integer, intent(out) :: nx, ny, nz
        nx = n_cells_x
        ny = n_cells_y
        nz = n_cells_z
    end subroutine get_cell_grid_dims

    subroutine get_cell(ix, iy, iz, cell)
        integer, intent(in) :: ix, iy, iz
        type(cell_type), intent(out) :: cell
        cell = cells(ix, iy, iz)
    end subroutine get_cell
    
    subroutine collect_cell_statistics()
        integer :: ix, iy, iz, empty_count, total_cells, total_atoms
        integer :: local_min_atoms, local_max_atoms
        
        ! Initialize counters
        empty_count = 0
        total_cells = n_cells_x * n_cells_y * n_cells_z
        total_atoms = 0
        local_min_atoms = huge(local_min_atoms)
        local_max_atoms = 0
        
        ! Count statistics
        do iz = 1, n_cells_z
            do iy = 1, n_cells_y
                do ix = 1, n_cells_x
                    ! Track empty cells
                    if (cells(ix,iy,iz)%n_atoms == 0) then
                        empty_count = empty_count + 1
                    end if
                    
                    ! Track min/max atoms in a cell
                    local_max_atoms = max(local_max_atoms, cells(ix,iy,iz)%n_atoms)
                    if (cells(ix,iy,iz)%n_atoms > 0) then
                        local_min_atoms = min(local_min_atoms, cells(ix,iy,iz)%n_atoms)
                    end if
                    
                    ! Track total atoms for average calculation
                    total_atoms = total_atoms + cells(ix,iy,iz)%n_atoms
                end do
            end do
        end do
        
        ! Update global statistics (running average)
        stats_collection_count = stats_collection_count + 1
        
        ! Update max and min (global maxima and minima)
        if (stats_collection_count == 1) then
            ! First collection, initialize values
            min_atoms_in_cell = local_min_atoms
            max_atoms_in_cell = local_max_atoms
            total_empty_cells = empty_count
            avg_atoms_per_cell = real(total_atoms, dp) / max(1, total_cells - empty_count)
        else
            ! Update running statistics
            min_atoms_in_cell = min(min_atoms_in_cell, local_min_atoms)
            max_atoms_in_cell = max(max_atoms_in_cell, local_max_atoms)
            
            ! Update running average of empty cells and atoms per cell
            total_empty_cells = nint((total_empty_cells * (stats_collection_count - 1) + &
                                    empty_count) / real(stats_collection_count, dp))
            
            ! Calculate average atoms per non-empty cell
            avg_atoms_per_cell = ((avg_atoms_per_cell * (stats_collection_count - 1)) + &
                                real(total_atoms, dp) / max(1, total_cells - empty_count)) / &
                                real(stats_collection_count, dp)
        end if
    end subroutine collect_cell_statistics
    
    function get_cell_statistics() result(stats)
        real(dp) :: stats(4)
        stats(1) = real(total_empty_cells, dp)
        stats(2) = avg_atoms_per_cell
        stats(3) = real(min_atoms_in_cell, dp)
        stats(4) = real(max_atoms_in_cell, dp)
    end function get_cell_statistics

    subroutine cleanup_cell_list()
        integer :: ix, iy, iz
        
        if (allocated(cells)) then
            !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(ix,iy,iz)
            do iz = 1, n_cells_z
                do iy = 1, n_cells_y
                    do ix = 1, n_cells_x
                        call cells(ix,iy,iz)%cleanup()
                    end do
                end do
            end do
            !$OMP END PARALLEL DO
            deallocate(cells)
        end if
        
        if (allocated(prev_coords)) deallocate(prev_coords)
    end subroutine cleanup_cell_list

end module cell_list_mod