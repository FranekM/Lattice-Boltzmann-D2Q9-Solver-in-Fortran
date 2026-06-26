module lbm_config
    use lbm_d2q9_consts
    use lbm_d2q9_types

    implicit none
    private


    public :: experiment_config
    public :: field_config
    public :: configs
    public :: config_count

    public :: read_configs
    public :: read_arguments


    type logger_config
        character(32) :: type = ""

        character(256) :: file = ""
        integer :: Nf = 0
        integer :: x = 0
        integer :: y = 0
    end type

    type gaussian_config
        real(dp) :: eps = 0.0_dp
        real(dp) :: val_0 = 0.0_dp
        real(dp) :: sigma = 0.0_dp
        integer :: xc = 0
        integer :: yc = 0
    end type

    type sine_config
        real(dp) :: val_0 = 0.0_dp
        real(dp) :: eps_x = 0.0_dp
        real(dp) :: eps_y = 0.0_dp
        integer :: mode_x = 0
        integer :: mode_y = 0
        integer :: nx = 0
        integer :: ny = 0
    end type

    type constant_config
        real(dp) :: value = 0.0_dp
    end type

    type field_config
        character(32) :: type = ""

        type(gaussian_config) :: gauss
        type(sine_config) :: sine
        type(constant_config) :: const
    end type

    type experiment_config
        character(64) :: family = ""
        character(64) :: name = ""

        ! simulation paramters
        integer :: Nx = 0
        integer :: Ny = 0
        integer :: Nt = 0
        real(dp) :: dpy = 0.0_dp
        real(dp) :: dpx = 0.0_dp

        ! logging paramters
        type(logger_config) :: log

        ! fields paramters
        type(field_config) :: rho
        type(field_config) :: ux
        type(field_config) :: uy

        ! boundry conditions
        character(32) :: bc_list(16) = ""
        integer :: bc_places(16) = 0
        real(dp) :: bc_speed(16) = 0.0_dp
        integer :: bc_count = 0

        ! omegas
        integer :: No = 0
        real(dp) :: omegas(25) = 0.0_dp

        ! python display paramters
        character(16) :: disp_x = ""
        character(16) :: disp_y = ""
    end type


    type(experiment_config) :: configs(100)
    integer :: config_count = 0

contains

    subroutine read_configs(selected_configs)
        character(*), intent(in) :: selected_configs(:)

        character(256) :: path
        integer :: i, unit, stat
        type(experiment_config) :: cfg
        namelist /experiment/ cfg

        do i = 1, size(selected_configs)
            write(path, '(A,A,A)') 'config/', trim(selected_configs(i)), '.in'

            open(               &
                newunit = unit, &
                file    = path, &
                status  = 'old' &
            )

            do
                cfg = experiment_config()
                read(                       &
                    unit,                   &
                    nml     = experiment,   &
                    iostat  = stat          &
                )
                if (stat /= 0) exit
                config_count = config_count + 1
                configs(config_count) = cfg
            end do

            close(unit)
        end do
    end subroutine


    subroutine read_arguments(args)
		character(:), allocatable, intent(out) :: args(:)
		integer :: i, argc, length

		argc = command_argument_count()
		allocate(character(128) :: args(argc))

		do i = 1, argc
			call get_command_argument(i, args(i))
		end do
	end subroutine

end module lbm_config