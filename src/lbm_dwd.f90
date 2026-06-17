module lbm_dwd
    use lbm_d2q9_consts
    use lbm_d2q9_solver
    use lbm_d2q9_boundry_conditions
    use lbm_simulation

    implicit none
    private

    public :: dwd_simulate

    public :: dwd_profile_logger
    public :: dwd_field_logger
    public :: dwd_amplitude_logger

    integer :: dwd_n
    real(dp) :: dwd_eps, dwd_rho, dwd_omega

contains

    real(dp) function dwd_rho_func(y, x)
        integer, intent(in) :: y, x
        real(dp) :: dim_x, k

        dim_x = real((E-W)+1, dp)
        k = 2.0_dp * PI / dim_x
        dwd_rho_func = dwd_rho + dwd_eps * sin((dwd_n * k) * (x-W))
    end function dwd_rho_func

    real(dp) function dwd_ux_func(y, x)
        integer, intent(in) :: y, x

        dwd_ux_func = 0.0_dp
    end function dwd_ux_func

    real(dp) function dwd_uy_func(y, x)
        integer, intent(in) :: y, x

        dwd_uy_func = 0.0_dp
    end function dwd_uy_func


	subroutine dwd_amplitude_logger(log_unit, t)
        integer, intent(in) :: log_unit, t
        real(dp) :: amplitude
        integer :: x, y

        y = nint(S + real((N-S+1) / 2, dp))
        x = nint(W + real((E-W+1) / (4.0_dp * dwd_n), dp))

        amplitude = abs(measure_rho(y, x))

        write(log_unit, '(ES12.5,A,I0,A,ES20.10)') &
                        dwd_omega,      char(9),   &
                        t,              char(9),   &
                        amplitude
    end subroutine dwd_amplitude_logger

    subroutine dwd_profile_logger(log_unit, t)
        integer, intent(in) :: log_unit, t
        real(dp) :: rho, ux, uy
        integer :: x, y_mid

        y_mid = W + (E-W)/(4*dwd_n)

        do x = W, E
            rho = measure_rho(y_mid, x)
            ux = measure_u(CX, y_mid, x, rho)
            uy = measure_u(CY, y_mid, x, rho)

            write(log_unit, '(I0,A,I0,A,ES20.10,A,ES20.10,A,ES20.10)') &
                            t,      char(9),       &
                            x,      char(9),       &
                            rho,    char(9),       &
                            uy,     char(9),       &
	                        ux
        end do
    end subroutine dwd_profile_logger

    subroutine dwd_field_logger(log_unit, t)
		integer, intent(in) :: log_unit, t
        real(dp) :: rho, uy, ux
		integer :: x, y

		do x = W, E
			do y = S, N
				rho = measure_rho(y, x)
				uy = measure_u(CY, y, x, rho)
				ux = measure_u(CX, y, x, rho)

				write(log_unit, '(I0,A,I0,A,I0,A,ES20.10,A,ES20.10,A,ES20.10)')   &
								t, char(9),                                       &
								y, char(9),                                       &
								x, char(9),                                       &
								rho, char(9),                                     &
								uy, char(9),                                      &
								ux
			end do
		end do

	end subroutine dwd_field_logger


    subroutine dwd_simulate(&
        grid_size_x,        &
        grid_size_y,        &
        set_n,              &
        set_eps,            &
        set_rho,            &
        set_omega,          &
        Nt,                 &
        log_step,           &
        logger,             &
        filename            &
    )
	    integer, intent(in) :: grid_size_x, grid_size_y, set_n, Nt, log_step
	    real(dp), intent(in) :: set_eps, set_rho, set_omega
        procedure(logger_func), optional, pointer :: logger
	    character(*), intent(in), optional :: filename
		integer :: i

        dwd_n = set_n
        dwd_eps = set_eps
        dwd_rho = set_rho
        dwd_omega = set_omega

	    call init_simulation(grid_size_x, grid_size_y, dwd_omega, dwd_rho_func, dwd_uy_func, dwd_ux_func, logger, filename)
	    call add_boundry_condition(periodic_bc, ALL_WALLS, ALL_CORNERS)
	    call run_simulation(Nt, log_step)
    end subroutine dwd_simulate

end module lbm_dwd