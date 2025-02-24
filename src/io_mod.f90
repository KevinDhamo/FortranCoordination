module io_mod
    use types_mod
    use trajectory_io_mod, only: init_trajectory_reader, read_trajectory_frame, &
                                count_trajectory_frames, cleanup_trajectory_reader
    use data_file_io_mod, only: read_data_file
    use output_io_mod, only: init_output, write_output_frame, &
                            write_final_report, cleanup_output
    use progress_mod, only: init_progress, update_progress_, &
                           finalize_progress
    implicit none
    private

    ! Re-export only necessary procedures
    public :: read_data_file          ! From data_file_io_mod
    public :: read_trajectory_frame   ! From trajectory_io_mod
    public :: write_output_frame      ! From output_io_mod
    public :: write_final_report      ! From output_io_mod
    public :: count_trajectory_frames ! From trajectory_io_mod
    public :: update_progress_       ! From progress_mod
    
    ! Initialize and cleanup procedures
    public :: initialize_io
    public :: cleanup_io

contains
    subroutine initialize_io(trajectory_file, output_file, n_frames)
        character(len=*), intent(in) :: trajectory_file, output_file
        integer, intent(in) :: n_frames
        
        call init_trajectory_reader(trajectory_file)
        call init_output(output_file)
        call init_progress(n_frames)
    end subroutine initialize_io

    subroutine cleanup_io()
        call cleanup_trajectory_reader()
        call cleanup_output()
        call finalize_progress()
    end subroutine cleanup_io

end module io_mod