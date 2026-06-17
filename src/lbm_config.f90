module lbm_config
    use lbm_d2q9_consts

    implicit none
    private

    public :: experiment_config
    public :: configs
    public :: config_count
    public :: read_configs

    type experiment_config
        character(64) :: name   = ""
        character(5)  :: type   = ""
        character(5)  :: output = ""

        integer :: Nx = 0
        integer :: Ny = 0

        integer  :: n    = 0
        real(dp) :: eps = 0.0_dp
        real(dp) :: rho = 0.0_dp

        integer  :: No         = 0
        real(dp) :: omegas(25) = 0.0_dp

        integer :: Nt     = 0
        integer :: frames = 0

        character(256) :: filename = ""

        character(16) :: disp_x = ""
        character(16) :: disp_y =""
    end type

    type(experiment_config) :: configs(100)
    integer :: config_count = 0

contains

    subroutine read_configs(path)
        character(*), intent(in) :: path
        integer :: unit, stat
        type(experiment_config) :: cfg
        namelist /experiment/ cfg

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
    end subroutine read_configs

end module lbm_config