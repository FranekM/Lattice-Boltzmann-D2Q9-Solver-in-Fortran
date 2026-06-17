module lbm_simulation
	use lbm_d2q9_consts
	use lbm_d2q9_solver
	use lbm_d2q9_boundry_conditions

	implicit none
	private

    public :: init_simulation
    public :: run_simulation
	public :: logger_func

	abstract interface
		subroutine logger_func(log_unit, t)
			integer, intent(in) :: log_unit, t
		end subroutine logger_func
	end interface

    real(dp) :: omega

	character(:), allocatable :: filename
	procedure(logger_func), pointer :: logger => null()
	logical :: logging_enabled = .false.

contains

	logical function is_logging_step(t, log_step)
		integer, intent(in) :: t, log_step

		is_logging_step =           &
			logging_enabled         &
			.and. (log_step > 0)    &
			.and. (mod(t, log_step) == 0)
	end function is_logging_step


    subroutine init_simulation( &
        grid_size_y,    &
        grid_size_x,    &
        set_omega,      &
        rho_func,       &
        uy_func,        &
        ux_func,        &
        set_logger,     &
        set_filename    &
    )
	    integer, intent(in) :: grid_size_y, grid_size_x
	    real(dp), intent(in) :: set_omega
		procedure(scalar_field) :: rho_func, uy_func, ux_func
	    procedure(logger_func), pointer, optional :: set_logger
	    character(*), intent(in), optional :: set_filename

	    nullify(logger)
		if (allocated(filename)) deallocate(filename)

	    if (present(set_filename)) filename = set_filename
	    if (present(set_logger)) logger => set_logger

	    omega = set_omega
	    logging_enabled = associated(logger) .and. allocated(filename)

        call init_lbm(1, 1, grid_size_y, grid_size_x, 1.0_dp)
        call set_initial_condition(rho_func, uy_func, ux_func)
	    call clear_boundry_conditions()
    end subroutine init_simulation


    subroutine run_simulation(Nt, log_step)
        integer, intent(in) :: Nt, log_step
        integer :: log_unit, t

        if (logging_enabled) open(   &
            newunit  = log_unit,     &
            file     = filename,     &
            status   = 'unknown',    &
            position = 'append'      &
        )

        do t = 0, Nt-1
            if (is_logging_step(t, log_step)) then
                call logger(log_unit, t)
            end if

            call collision_step(omega)
            call streaming_step()
            call communication_step()
        end do

        if (logging_enabled) close(log_unit)
    end subroutine run_simulation

end module lbm_simulation