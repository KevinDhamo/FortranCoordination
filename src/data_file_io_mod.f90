module data_file_io_mod
    use types_mod
    use config_mod
    use error_mod
    implicit none
    private

    ! Public procedures
    public :: read_data_file

contains
    subroutine read_data_file(filename, atom_info, n_atoms, n_types, box_length)
        character(len=*), intent(in) :: filename
        type(atom_type_info), allocatable, intent(out) :: atom_info(:)
        integer, intent(out) :: n_atoms, n_types
        real(dp), intent(out) :: box_length(3)
        
        integer :: data_unit, io_stat, i
        character(256) :: line
        integer :: type_id
        real(dp) :: mass
        character(20) :: element_name
        real(dp) :: xlo, xhi, ylo, yhi, zlo, zhi
        logical :: found_atoms = .false.
        logical :: found_types = .false.
        logical :: found_masses = .false.
        logical :: found_box = .false.
        
        ! Check if file exists
        inquire(file=filename, exist=found_atoms)
        if (.not. found_atoms) then
            call handle_error("Data file not found: " // trim(filename), ERR_FILE_IO)
        end if
        
        ! Open data file
        open(newunit=data_unit, file=filename, status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) call handle_error("Failed to open data file", ERR_FILE_IO)
        
        ! First pass: read header information
        do
            read(data_unit, '(A)', iostat=io_stat) line
            if (io_stat /= 0) exit
            
            ! Skip empty lines and comments
            if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
            
            ! Parse essential header information
            if (index(line, 'atoms') > 0 .and. index(line, 'atom types') == 0) then
                read(line, *, iostat=io_stat) n_atoms
                found_atoms = (io_stat == 0)
            else if (index(line, 'atom types') > 0) then
                read(line, *, iostat=io_stat) n_types
                found_types = (io_stat == 0)
            else if (index(line, 'xlo xhi') > 0) then
                read(line, *, iostat=io_stat) xlo, xhi
                found_box = (found_box .or. io_stat == 0)
            else if (index(line, 'ylo yhi') > 0) then
                read(line, *, iostat=io_stat) ylo, yhi
                found_box = (found_box .or. io_stat == 0)
            else if (index(line, 'zlo zhi') > 0) then
                read(line, *, iostat=io_stat) zlo, zhi
                found_box = (found_box .or. io_stat == 0)
            else if (index(line, 'Masses') > 0) then
                found_masses = .true.
                exit
            end if
        end do
        
        ! Verify we found all required information
        if (.not. (found_atoms .and. found_types .and. found_box .and. found_masses)) then
            call handle_error("Missing required information in data file", ERR_FILE_IO)
        end if
        
        ! Calculate box dimensions
        box_length(1) = xhi - xlo
        box_length(2) = yhi - ylo
        box_length(3) = zhi - zlo
        
        ! Allocate atom_info array
        if (allocated(atom_info)) deallocate(atom_info)
        allocate(atom_info(n_types), stat=io_stat)
        call check_allocation(io_stat, "atom_info")
        
        ! Initialize atom_info with defaults
        do i = 1, n_types
            atom_info(i)%type_id = i
            atom_info(i)%include = .true.
            write(atom_info(i)%name, '(A,I0)') 'Type', i
            atom_info(i)%mass = 0.0_dp
        end do
        
        ! Skip blank line after Masses header
        read(data_unit, *)
        
        ! Read Masses section
        do i = 1, n_types
            read(data_unit, '(A)', iostat=io_stat) line
            if (io_stat /= 0) call handle_error("Error reading Masses section", ERR_FILE_IO)
            
            ! Parse mass line: type_id mass # element_name
            read(line, *, iostat=io_stat) type_id, mass
            if (io_stat /= 0) cycle
            
            ! Store mass
            if (type_id > 0 .and. type_id <= n_types) then
                atom_info(type_id)%mass = mass
                
                ! Try to extract element name from comment if present
                if (index(line, '#') > 0) then
                    line = adjustl(line(index(line, '#')+1:))
                    read(line, *, iostat=io_stat) element_name
                    if (io_stat == 0) atom_info(type_id)%name = trim(adjustl(element_name))
                end if
            end if
        end do
        
        close(data_unit)
        
        if (VERBOSE) then
            write(*,'(A)') ' System information from data file:'
            write(*,'(A,I0)') '   Number of atoms: ', n_atoms
            write(*,'(A,I0)') '   Number of atom types: ', n_types
            write(*,'(A)') '   Box dimensions:'
            write(*,'(A,F10.4)') '     X length: ', box_length(1)
            write(*,'(A,F10.4)') '     Y length: ', box_length(2)
            write(*,'(A,F10.4)') '     Z length: ', box_length(3)
            write(*,'(A)') '   Atom types:'
            do i = 1, n_types
                write(*,'(A,A,A,I0,A,F10.4)') '     Type ', &
                    trim(atom_info(i)%name), ' (', i, ') mass: ', atom_info(i)%mass
            end do
        end if
    end subroutine read_data_file

end module data_file_io_mod