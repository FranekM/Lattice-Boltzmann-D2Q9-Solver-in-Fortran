module lbm_simulation
	use lbm_d2q9_boundry_conditions
	use lbm_d2q9_consts
	use lbm_d2q9_solver
	use lbm_d2q9_types

	implicit none
	private

    public :: init_simulation
    public :: run_simulation
	public :: end_simulation

	class(logger), allocatable :: log
	logical :: logging_enabled = .false.
	character(256) :: log_file

	real :: start, finish

contains

	logical function is_debug_step(t, debug_step)
		integer, intent(in) :: t, debug_step

		is_debug_step = .false.

		if (this_image() /= 1) return
		if (t <= 0)            return
		if (debug_step <= 0)   return

		is_debug_step = (mod(t, debug_step) == 0)
	end function

	logical function is_logging_step(t, log_step)
		integer, intent(in) :: t, log_step

		is_logging_step = .false.

		if (.not. logging_enabled) return
		if (log_step <= 0)         return

		is_logging_step = (mod(t, log_step) == 0)
	end function

	function format_filename() result(f)
		character(256) :: f

		write(f, '(A,I0,A)') trim(log%file), this_image(), trim(log%ext)
	end function

	subroutine debug_progress(Nt, t)
		integer, intent(in) :: Nt, t
		real :: estimate

		call cpu_time(finish)
		estimate = ((real(Nt) / real(t)) - 1.0) * (finish - start)
		print '(A,I0,A,I0,A,F0.2,A,A,F0.2,A)',          &
			" - ",            t, "/", Nt,               &
			" in: ",          (finish - start),  "s,",  &
			" est. end in: ", estimate,          "s"
	end subroutine


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
	    class(logger), allocatable, intent(in) :: set_logger

	    type(boundry_condition) :: internal_bc

	    logging_enabled = allocated(set_logger)
	    if (logging_enabled) then
		    log = set_logger
			log_file = format_filename()
	    end if

	    ! TODO: move to configuration
	    internal_bc%apply    => image_sync_bc
	    internal_bc%postproc => null()
	    internal_bc%places   = 0
	    internal_bc%speed    = 0.0_dp

        call init_lbm(grid_size_y, grid_size_x, omega, dpy, dpx)
        call set_initial_condition(rho, uy, ux)
	    call set_boundry_condition(bc_list, bc_count, internal_bc)
    end subroutine


    subroutine run_simulation(Nt, log_step, debug_step)
        integer, intent(in) :: Nt, log_step, debug_step
        integer :: log_unit, t

        real(dp) :: rho_w
        integer :: x

        if (logging_enabled) open(   &
            newunit  = log_unit,     &
            file     = log_file,     &
            status   = 'unknown',    &
            position = 'append'      &
        )

        call cpu_time(start)

        do t = 0, Nt
            if (is_logging_step(t, log_step)) &
                call log%log(log_unit, t)
            if (is_debug_step(t, debug_step)) &
                call debug_progress(Nt, t)

            call streaming_step()
            sync all
            call communication_step()
            call collision_step()

        end do

        if (logging_enabled) close(log_unit)
    end subroutine

	subroutine end_simulation()

		if (allocated(log)) deallocate(log)
		logging_enabled = .false.
		log_file = ""
		finish = 0
		start = 0

		call end_lbm()
	end subroutine

end module lbm_simulation