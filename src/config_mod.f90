module config_mod
    use types_mod
    use error_mod
    implicit none
    private

    ! Public parameters and procedures
    public :: initialize_config
    public :: read_config
    public :: get_config_value
    public :: cleanup_config

    ! File names and paths
    character(len=*), parameter, public :: INPUT_FILE = 'trajectory.dump'
    character(len=*), parameter, public :: OUTPUT_FILE = 'coordination_numbers.dat'
    character(len=*), parameter, public :: DATA_FILE = 'lammps.data'
    character(len=*), parameter, public :: CONFIG_FILE = 'analysis_config.txt'

    ! Default parameters
    real(dp), parameter, public :: DEFAULT_CUTOFF = 3.0_dp
    integer, parameter, public :: PROGRESS_BAR_WIDTH = 50
    logical, parameter, public :: VERBOSE = .true.

    ! Cell list configuration
    integer, parameter, public :: CELL_UPDATE_FREQ = 5
    real(dp), parameter, public :: CELL_SIZE_FACTOR = 1.0_dp
    integer, parameter, public :: MAX_ATOMS_PER_CELL = 50
    real(dp), parameter, public :: REBUILD_THRESHOLD = 0.1_dp

    ! Type definition for configuration values - must be defined before use
    type :: config_value_type
        character(len=32) :: key = ''
        character(len=128) :: value = ''
    end type config_value_type

    ! Module variables
    type(config_value_type), allocatable :: config_data(:)

contains

    subroutine initialize_config()
        integer :: n_config_items = 20  ! Initial size
        integer :: alloc_stat
        
        if (allocated(config_data)) deallocate(config_data)
        allocate(config_data(n_config_items), stat=alloc_stat)
        call check_allocation(alloc_stat, "config_data")
        
        ! Initialize all entries to empty strings
        config_data%key = ''
        config_data%value = ''
    end subroutine initialize_config

    subroutine read_config()
        integer :: io_stat, unit
        character(len=256) :: line
        character(len=32) :: key
        character(len=128) :: value
        logical :: file_exists
        
        inquire(file=CONFIG_FILE, exist=file_exists)
        if (.not. file_exists) return
        
        open(newunit=unit, file=CONFIG_FILE, status='old', action='read', iostat=io_stat)
        if (io_stat /= 0) then
            call handle_error("Failed to open configuration file", ERR_FILE_IO)
            return
        end if
        
        do
            read(unit, '(A)', iostat=io_stat) line
            if (io_stat /= 0) exit
            
            ! Skip comments and empty lines
            if (line(1:1) == '#' .or. len_trim(line) == 0) cycle
            
            ! Parse key-value pair
            call parse_config_line(line, key, value)
            if (len_trim(key) > 0) then
                call store_config_value(key, value)
            end if
        end do
        
        close(unit)
    end subroutine read_config

    subroutine parse_config_line(line, key, value)
        character(len=*), intent(in) :: line
        character(len=*), intent(out) :: key
        character(len=*), intent(out) :: value
        integer :: eq_pos
        
        eq_pos = index(line, '=')
        if (eq_pos > 0) then
            key = adjustl(line(:eq_pos-1))
            value = adjustl(line(eq_pos+1:))
        else
            key = ''
            value = ''
        end if
    end subroutine parse_config_line

    subroutine store_config_value(key, value)
        character(len=*), intent(in) :: key
        character(len=*), intent(in) :: value
        integer :: i, empty_slot
        
        ! Find existing key or first empty slot
        empty_slot = -1
        do i = 1, size(config_data)
            if (config_data(i)%key == key) then
                config_data(i)%value = value
                return
            else if (empty_slot == -1 .and. len_trim(config_data(i)%key) == 0) then
                empty_slot = i
            end if
        end do
        
        ! Store in empty slot if found
        if (empty_slot > 0) then
            config_data(empty_slot)%key = key
            config_data(empty_slot)%value = value
        else
            call handle_error("Configuration storage full", ERR_INVALID_PARAM, fatal=.false.)
        end if
    end subroutine store_config_value

    function get_config_value(key, default) result(value)
        character(len=*), intent(in) :: key
        character(len=*), intent(in) :: default
        character(len=128) :: value
        integer :: i
        
        value = default
        if (.not. allocated(config_data)) return
        
        do i = 1, size(config_data)
            if (config_data(i)%key == key) then
                value = config_data(i)%value
                return
            end if
        end do
    end function get_config_value

    subroutine cleanup_config()
        if (allocated(config_data)) deallocate(config_data)
    end subroutine cleanup_config

end module config_mod