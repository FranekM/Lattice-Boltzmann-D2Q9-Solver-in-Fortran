program LatticeBoltzmannMethod_D2Q9
    use lbm_d2q9_consts
    use lbm_d2q9_types
    use lbm_simulation
    use lbm_faktory
    use lbm_loggers
    use lbm_fields
    use lbm_config
    use lbm_utils

    implicit none

    integer :: i, j, k
    real(dp) :: omega

    class(logger), allocatable :: log
    class(scalar_field), allocatable :: rho
    class(scalar_field), allocatable :: uy
    class(scalar_field), allocatable :: ux
    type(boundry_condition) :: bc_list(16)

    character(:), allocatable :: selected_configs(:)

    call read_arguments(selected_configs)
    call read_configs(selected_configs)

    do i = 1, config_count
        if (.not. prepare_output_file(configs(i)%log%file)) cycle

        print *, 'Starting experiment: ', trim(configs(i)%name), ' from family: ', trim(configs(i)%family)

        do j = 1, configs(i)%No
            omega = configs(i)%omegas(j)
            print *, ' > omega = ', omega

            if (allocated(log)) deallocate(log)
            if (allocated(rho)) deallocate(rho)
            if (allocated(ux))  deallocate(ux)
            if (allocated(uy))  deallocate(uy)
            do k = 1, size(bc_list)
                nullify(bc_list(k)%apply)
                bc_list(k)%places = 0
            end do

            call build_logger(log, configs(i)%log%type, omega, i)
            call build_field(rho, configs(i)%rho, i)
            call build_field(ux, configs(i)%ux, i)
            call build_field(uy, configs(i)%uy, i)
            call build_bc_list(bc_list, i)

            call init_simulation(                       &
                grid_size_y     = configs(i)%Ny,        &
                grid_size_x     = configs(i)%Nx,        &
                rho             = rho,                  &
                uy              = uy,                   &
                ux              = ux,                   &
                omega           = omega,                &
                dpy             = configs(i)%dpy,       &
                dpx             = configs(i)%dpx,       &
                bc_list         = bc_list,              &
                bc_count        = configs(i)%bc_count,  &
                set_logger      = log                   &
            )

            call run_simulation(                            &
                Nt       = configs(i)%Nt,                   &
                log_step = configs(i)%Nt / configs(i)%log%Nf&
            )
        end do

        print *, 'Experiment: ', trim(configs(i)%name), ' finished'
    end do

end program LatticeBoltzmannMethod_D2Q9