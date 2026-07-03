module lbm_loggers
	use lbm_d2q9_types
	use lbm_d2q9_consts
	use lbm_d2q9_solver

	implicit none
	private


	public :: point_logger
	public :: profile_x_logger
	public :: profile_y_logger
	public :: field_logger


	type, extends(logger) :: point_logger
	    real(dp) :: omega
	    integer :: y, x
	contains
	    procedure :: log => point_log
	end type

	type, extends(logger) :: profile_x_logger
		real(dp) :: omega
		integer :: y
	contains
		procedure :: log => profile_x_log
	end type

	type, extends(logger) :: profile_y_logger
		real(dp) :: omega
		integer :: x
	contains
		procedure :: log => profile_y_log
	end type

	type, extends(logger) :: field_logger
		real(dp) :: omega
	contains
		procedure :: log => field_log
	end type

contains


	subroutine point_log(this, log_unit, t)
		class(point_logger), intent(inout) :: this
        integer, intent(in) :: log_unit, t
        real(dp) :: rho, uy, ux
		integer :: y, x, gy, gx     !! local and global indices

		gy = this%y
		gx = this%x
		if (.not. is_global_y_in(gy)) return
		if (.not. is_global_x_in(gx)) return
		y = global_to_local_y(gy)
		x = global_to_local_x(gx)

        rho = measure_rho(y, x)
        uy = measure_u(CY, y, x, rho)
        ux = measure_u(CX, y, x, rho)

        write(log_unit, '(ES12.5,A,I0,A,I0,A,I0,A,ES20.10,A,ES20.10,A,ES20.10)')  &
                        this%omega, char(9),        &
                        t,          char(9),        &
                        gy,         char(9),        &
                        gx,         char(9),        &
                        rho,        char(9),        &
                        uy,         char(9),        &
                        ux
    end subroutine

	subroutine profile_x_log(this, log_unit, t)
		class(profile_x_logger), intent(inout) :: this
        integer, intent(in) :: log_unit, t
        real(dp) :: rho, ux, uy
        integer :: y, x, gy, gx        !! local and global indices

		gy = this%y
		if (.not. is_global_y_in(gy)) return
		y = global_to_local_y(gy)

        do x = W, E
            rho = measure_rho(y, x)
            ux = measure_u(CX, y, x, rho)
            uy = measure_u(CY, y, x, rho)

            gx = local_to_global_x(x)
            write(log_unit, '(ES12.5,A,I0,A,I0,A,I0,A,ES20.10,A,ES20.10,A,ES20.10)') &
                            this%omega, char(9),       &
                            t,          char(9),       &
                            gy,         char(9),       &
                            gx,         char(9),       &
                            rho,        char(9),       &
                            uy,         char(9),       &
	                        ux
        end do
    end subroutine

	subroutine profile_y_log(this, log_unit, t)
		class(profile_y_logger), intent(inout) :: this
        integer, intent(in) :: log_unit, t
        real(dp) :: rho, ux, uy
        integer :: y, x, gy, gx        !! local and global indices

		gx = this%x
		if (.not. is_global_x_in(gx)) return
		x = global_to_local_x(gx)

        do y = S, N
            rho = measure_rho(y, x)
            ux = measure_u(CX, y, x, rho)
            uy = measure_u(CY, y, x, rho)

	        gy = local_to_global_y(y)
            write(log_unit, '(ES12.5,A,I0,A,I0,A,I0,A,ES20.10,A,ES20.10,A,ES20.10)') &
                            this%omega, char(9),       &
                            t,          char(9),       &
                            gy,         char(9),       &
                            gx,         char(9),       &
                            rho,        char(9),       &
                            uy,         char(9),       &
	                        ux
        end do
    end subroutine

	subroutine field_log(this, log_unit, t)
		class(field_logger), intent(inout) :: this
		integer, intent(in) :: log_unit, t
        real(dp) :: rho, uy, ux
		integer :: y, x, gy, gx   !! local and global indices

		do y = S, N
			do x = W, E
 				rho = measure_rho(y, x)
				uy = measure_u(CY, y, x, rho)
				ux = measure_u(CX, y, x, rho)

				gy = local_to_global_y(y)
				gx = local_to_global_x(x)
				write(log_unit, '(ES12.5,A,I0,A,I0,A,I0,A,ES20.10,A,ES20.10,A,ES20.10)') &
                                this%omega, char(9),        &
								t,          char(9),        &
								gy,         char(9),        &
								gx,         char(9),        &
								rho,        char(9),        &
								uy,         char(9),        &
								ux
			end do
		end do
	end subroutine


end module lbm_loggers