module lbm_swd
    use lbm_d2q9_consts
    use lbm_d2q9_solver
    use lbm_d2q9_boundry_conditions
    use lbm_simulation

    implicit none
    private

    public :: swd_simulate

    public :: swd_amplitude_logger
    public :: swd_profile_logger
    public :: swd_field_logger

    integer :: swd_n
    real(dp):: swd_eps, swd_rho, swd_omega

contains

    real(dp) function swd_rho_func(y, x)
        integer, intent(in) :: y, x

        swd_rho_func = swd_rho
    end function swd_rho_func

    real(dp) function swd_ux_func(y, x)
        integer, intent(in) :: y, x
        real(dp) :: dim_y, k

        dim_y = real(N-S+1, dp)
        k = 2.0_dp * PI / dim_y
        swd_ux_func = swd_eps * sin((swd_n * k) * (y-S))
    end function swd_ux_func

    real(dp) function swd_uy_func(y, x)
        integer, intent(in) :: y, x

        swd_uy_func = 0.0_dp
    end function swd_uy_func


    subroutine swd_amplitude_logger(log_unit, t)
        integer, intent(in) :: log_unit, t
        real(dp) :: rho, amplitude
        integer :: x, y

        x = nint(W + real((E-W+1) / 2, dp))
        y = nint(S + real((N-S+1) / (4.0_dp * swd_n), dp))

        rho = measure_rho(y, x)
        amplitude = abs(measure_u(CX, y, x, rho))

        write(log_unit, '(ES12.5,A,I0,A,ES20.10)') &
                        swd_omega,      char(9),   &
                        t,              char(9),   &
                        amplitude
    end subroutine swd_amplitude_logger

	subroutine swd_profile_logger(log_unit, t)
        integer, intent(in) :: log_unit, t
        real(dp) :: rho, ux, uy
        integer :: y, x_mid

        x_mid = W + (E-W) / (4*n)

        do y = S, N
            rho = measure_rho(y, x_mid)
            ux = measure_u(CX, y, x_mid, rho)
            uy = measure_u(CY, y, x_mid, rho)

            write(log_unit, '(I0,A,I0,A,ES20.10,A,ES20.10,A,ES20.10)') &
                            t,      char(9),       &
                            y,      char(9),       &
                            rho,    char(9),       &
                            uy,     char(9),       &
	                        ux
        end do
    end subroutine swd_profile_logger

	subroutine swd_field_logger(log_unit, t)
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
	end subroutine swd_field_logger


    subroutine swd_simulate(&
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

	    swd_n = set_n
	    swd_eps = set_eps
	    swd_rho = set_rho
	    swd_omega = set_omega

	    call init_simulation(grid_size_y, grid_size_x, swd_omega, swd_rho_func, swd_uy_func, swd_ux_func, logger, filename)
	    call add_boundry_condition(periodic_bc, ALL_WALLS, ALL_CORNERS)
	    call run_simulation(Nt, log_step)
    end subroutine swd_simulate

end module lbm_swd