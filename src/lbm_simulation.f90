module lbm_simulation
	use lbm_d2q9_boundry_conditions
	use lbm_d2q9_consts
	use lbm_d2q9_solver
	use lbm_d2q9_types

	implicit none
	private

    public :: init_simulation
    public :: run_simulation

	class(logger), allocatable :: log
	logical :: logging_enabled = .false.

contains

	logical function is_logging_step(t, log_step)
		integer, intent(in) :: t, log_step

		is_logging_step =           &
			logging_enabled         &
			.and. (log_step > 0)    &
			.and. (mod(t, log_step) == 0)
	end function


    subroutine init_simulation( &
        grid_size_y,    &
        grid_size_x,    &
        rho,            &
        uy,             &
        ux,             &
	    omega,          &
	    dpy,            &
	    dpx,            &
	    bc_list,        &
        bc_count,       &
        set_logger      &
    )
	    integer, intent(in) :: grid_size_y, grid_size_x, bc_count
		class(scalar_field), intent(in) :: rho, uy, ux
	    real(dp), intent(in) :: omega, dpy, dpx
	    type(boundry_condition), intent(in) :: bc_list(*)
	    class(logger), intent(in), optional :: set_logger

	    logging_enabled = present(set_logger)
	    if (logging_enabled) log = set_logger

        call init_lbm(grid_size_y, grid_size_x, omega, dpy, dpx, bc_list, bc_count)
        call set_initial_condition(rho, uy, ux)
    end subroutine


    subroutine run_simulation(Nt, log_step)
        integer, intent(in) :: Nt, log_step
        integer :: log_unit, t

        real(dp) :: rho_w
        integer :: x

        if (logging_enabled) open(   &
            newunit  = log_unit,     &
            file     = log%file,     &
            status   = 'unknown',    &
            position = 'append'      &
        )

        do t = 0, Nt
            if (is_logging_step(t, log_step)) then
                call log%log(log_unit, t)
            end if

            call streaming_step()
            call communication_step()
            call collision_step()

        end do

        if (logging_enabled) close(log_unit)
    end subroutine

end module lbm_simulation