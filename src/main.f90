program coordination_analysis
    use types_mod
    use config_mod
    use cell_list_mod
    use coordination_mod
    use io_mod
    use error_mod
    use benchmark_mod
    use omp_lib
    implicit none

    ! System variables
    integer :: n_atoms = 0
    integer :: n_types = 0
    type(atom_type_info), allocatable :: atom_info(:)
    real(dp) :: box_length(3)
    integer, parameter :: MAX_FRAMES = 10000
    
    ! Atomic data
    real(dp), allocatable :: coords(:,:)
    integer, allocatable :: atom_types(:)
    character(len=8), allocatable :: elements(:)
    
    ! Analysis variables
    integer :: frame = 0
    integer :: total_frames = 0
    logical :: eof = .false.
    real(dp) :: start_time, end_time
    logical, allocatable :: include_mask(:)
    
    ! Local variables
    integer :: io_stat, num_threads
    logical :: benchmark_mode = .false.
    character(len=32) :: arg
    
    ! Check command line arguments for benchmark mode
    if (command_argument_count() > 0) then
        call get_command_argument(1, value=arg)
        benchmark_mode = (trim(arg) == '--benchmark')
    end if
    
    ! Initialize OpenMP (only if not in benchmark mode)
    if (.not. benchmark_mode) then
        num_threads = omp_get_max_threads()
        call omp_set_num_threads(num_threads)
        if (VERBOSE) write(*,'(A,I4,A)') ' Using', num_threads, ' OpenMP threads'
    end if
    
    ! Initialize configuration
    call initialize_config()
    call read_config()
    
    ! Read system information
    call read_data_file(DATA_FILE, atom_info, n_atoms, n_types, box_length)
    
    ! Allocate arrays with error checking
    allocate(coords(n_atoms,3), atom_types(n_atoms), elements(n_atoms), &
             include_mask(n_types), stat=io_stat)
    if (io_stat /= 0) call handle_error("Failed to allocate arrays", ERR_ALLOCATION)
    
    ! Initialize arrays
    coords = 0.0_dp
    atom_types = 0
    elements = ''
    include_mask = atom_info%include
    
    ! Count frames
    if (VERBOSE) write(*,'(A)') ' Counting frames in trajectory file...'
    total_frames = count_trajectory_frames(INPUT_FILE)
    total_frames = min(MAX_FRAMES, total_frames)
    
    if (benchmark_mode) then
        ! Run benchmark mode
        call initialize_benchmark()
        call run_benchmark(coords, atom_types, elements, box_length, &
                         n_atoms, n_types, atom_info, total_frames)
        call cleanup_benchmark()
    else
        ! Normal analysis mode
        ! Set up analysis
        call initialize_coordination(n_atoms, n_types, atom_info)
        call initialize_cell_list(box_length, get_cutoffs(), n_atoms)
        
        ! Initialize I/O
        call initialize_io(INPUT_FILE, OUTPUT_FILE, total_frames)
        
        ! Process trajectory
        start_time = omp_get_wtime()
        
        do while (.not. eof .and. frame < total_frames)
            call read_trajectory_frame(coords, atom_types, elements, box_length, &
                                    frame, eof, atom_info)
            if (eof) exit
            
            call calculate_coordination(coords, atom_types, n_atoms, box_length, &
                                     frame, include_mask)
            
            call write_output_frame(frame, elements, atom_types, coord_numbers, n_atoms)
            call update_progress_(frame)
        end do
        
        ! Finalize
        end_time = omp_get_wtime()
        
        ! Write final report before cleanup
        call write_final_report(frame, n_atoms, end_time - start_time, atom_info, &
                              n_types, coord_numbers, atom_types)
        
        ! Cleanup analysis
        call cleanup_io()
        call cleanup_coordination()
        call cleanup_cell_list()
    end if
    
    ! Final cleanup
    call cleanup_config()
    
    if (allocated(coords)) deallocate(coords)
    if (allocated(atom_types)) deallocate(atom_types)
    if (allocated(elements)) deallocate(elements)
    if (allocated(include_mask)) deallocate(include_mask)
    if (allocated(atom_info)) deallocate(atom_info)

contains
    function get_cutoffs() result(cutoffs)
        real(dp), allocatable :: cutoffs(:)
        integer :: i
        
        allocate(cutoffs(n_pairs))
        do i = 1, n_pairs
            cutoffs(i) = pairs(i)%cutoff
        end do
    end function get_cutoffs

end program coordination_analysis