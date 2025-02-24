module progress_mod
    use types_mod, only: dp
    use config_mod, only: PROGRESS_BAR_WIDTH
    use omp_lib
    implicit none
    private

    ! Public procedures
    public :: init_progress
    public :: display_progress
    public :: update_progress_
    public :: finalize_progress

    ! Module variables
    integer, save :: total_frames = 0
    real(dp), save :: start_time = 0.0_dp
    real(dp), allocatable, save :: frame_times(:)

contains
    subroutine init_progress(n_frames)
        integer, intent(in) :: n_frames
        integer :: alloc_stat

        total_frames = n_frames
        start_time = omp_get_wtime()

        if (allocated(frame_times)) deallocate(frame_times)
        allocate(frame_times(n_frames), stat=alloc_stat)
        if (alloc_stat /= 0) then
            write(*,*) "Warning: Could not allocate frame timing array"
            return
        end if
        frame_times = 0.0_dp
    end subroutine init_progress

    subroutine display_progress(frame, elapsed_time)
        integer, intent(in) :: frame
        real(dp), intent(in) :: elapsed_time
        integer :: progress, width
        character(len=30) :: time_string
        character(len=20) :: frame_string
        real(dp) :: frames_per_second, estimated_remaining
        
        ! Calculate progress width
        width = PROGRESS_BAR_WIDTH - 1
        
        ! Handle edge case where total_frames might be zero
        if (total_frames <= 0) then
            progress = 0
        else
            progress = min(width, nint(real(frame) / real(total_frames) * width))
        end if
        
        ! Calculate performance metrics
        if (elapsed_time > 0.0_dp) then
            frames_per_second = real(frame) / elapsed_time
            if (total_frames > 0) then
                estimated_remaining = (real(total_frames) - real(frame)) / frames_per_second
            else
                estimated_remaining = 0.0_dp
            end if
        else
            frames_per_second = 0.0_dp
            estimated_remaining = 0.0_dp
        end if
        
        ! Format frame counter
        if (total_frames > 0) then
            write(frame_string,'(I0,A,I0)') frame, '/', total_frames
        else
            write(frame_string,'(I0)') frame
        end if
        frame_string = adjustr(frame_string(1:min(11,len_trim(frame_string))))
        
        ! Format time
        write(time_string,'(F8.3,A)') elapsed_time, 's'
        
        ! Display progress bar
        write(*,'(a1,1x,A11,1x)',advance='no') achar(13), frame_string
        write(*,'(A)',advance='no') '['
        write(*,'(A)',advance='no') repeat('#', progress)
        write(*,'(A)',advance='no') repeat(' ', width - progress)
        write(*,'(A)',advance='no') ']'
        
        ! Display timing information
        write(*,'(1x,A,A)',advance='no') trim(adjustl(time_string))
        if (frames_per_second > 0.0_dp) then
            write(*,'(A,F6.2,A)',advance='no') ' [', frames_per_second, ' frames/s]'
        end if
        if (estimated_remaining > 0.0_dp) then
            write(*,'(A,F6.1,A)',advance='no') ' [ETA: ', estimated_remaining, 's]'
        end if
        
        call flush(6)
    end subroutine display_progress

    subroutine update_progress_(frame)
        integer, intent(in) :: frame
        real(dp) :: current_time
        
        current_time = omp_get_wtime()
        if (allocated(frame_times)) then
            frame_times(frame + 1) = current_time
        end if
        
        call display_progress(frame, current_time - start_time)
    end subroutine update_progress_

    subroutine finalize_progress()
        if (allocated(frame_times)) deallocate(frame_times)
    end subroutine finalize_progress

end module progress_mod