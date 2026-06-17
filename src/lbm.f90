program LatticeBoltzmannMethod_D2Q9
    use lbm_d2q9_consts
    use lbm_simulation
    use lbm_config
    use lbm_swd
    use lbm_dwd
    use lbm_io

    implicit none

    integer :: i, o
    integer :: log_step
    real(dp) :: omega
    procedure(logger_func), pointer :: logger => null()

    call read_configs('config.in')

    do i = 1, config_count
        print *, 'Starting experiment: ', trim(configs(i)%name)

        if (.not. prepare_output_file(configs(i)%filename)) cycle

        do o = 1, configs(i)%No
            omega = configs(i)%omegas(o)
            print *, ' - simulating with omega = ', omega

            select case(configs(i)%type)
                case("swd")
                    if (configs(i)%output == "fld") logger => swd_field_logger
                    if (configs(i)%output == "prf") logger => swd_profile_logger
                    if (configs(i)%output == "amp") logger => swd_amplitude_logger
                    log_step = configs(i)%Nt / configs(i)%frames

                    call swd_simulate(                           &
                        grid_size_x     = configs(i)%Nx,         &
                        grid_size_y     = configs(i)%Ny,         &
                        set_n           = configs(i)%n,          &
                        set_eps         = configs(i)%eps,        &
                        set_rho         = configs(i)%rho,        &
                        set_omega       = omega,                 &
                        Nt              = configs(i)%Nt,         &
                        log_step        = log_step,              &
                        logger          = logger,                &
                        filename        = configs(i)%filename    &
                    )

                case("dwd")
                    if (configs(i)%output == "fld") logger => dwd_field_logger
                    if (configs(i)%output == "prf") logger => dwd_profile_logger
                    if (configs(i)%output == "amp") logger => dwd_amplitude_logger
                    log_step = configs(i)%Nt / configs(i)%frames

                    call dwd_simulate(                           &
                        grid_size_x     = configs(i)%Nx,         &
                        grid_size_y     = configs(i)%Ny,         &
                        set_n           = configs(i)%n,          &
                        set_eps         = configs(i)%eps,        &
                        set_rho         = configs(i)%rho,        &
                        set_omega       = omega,                 &
                        Nt              = configs(i)%Nt,         &
                        log_step        = log_step,              &
                        logger          = logger,                &
                        filename        = configs(i)%filename    &
                    )
            end select
        end do
        print *, 'Experiment: ', trim(configs(i)%name), ' finished'
    end do

end program LatticeBoltzmannMethod_D2Q9