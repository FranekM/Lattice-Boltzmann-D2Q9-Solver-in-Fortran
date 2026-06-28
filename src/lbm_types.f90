module lbm_d2q9_types
	use iso_fortran_env, only: dp => real64

	implicit none


	! SCALAR FIELDS
	type, abstract :: scalar_field
	contains
	    procedure(eval_interface), deferred :: eval
	end type

	abstract interface
	    real(dp) function eval_interface(this, y, x)
	        import scalar_field, dp

	        class(scalar_field), intent(in) :: this
	        integer, intent(in) :: y, x
	    end function
	end interface


	! LOGGERS
	type, abstract :: logger
		character(:), allocatable :: file
	contains
		procedure(log_interface), deferred :: log
	end type

	abstract interface
		subroutine log_interface(this, log_unit, t)
			import logger, dp

			class(logger), intent(inout) :: this
			integer, intent(in) :: log_unit, t
		end subroutine
	end interface


	! BOUNDRY CONDITIONS
	abstract interface
        subroutine bc_func(places)
            integer, intent(in) :: places
        end subroutine

		subroutine pp_func(places, speed)
			import dp
			integer, intent(in) :: places
			real(dp), intent(in) :: speed
        end subroutine
    end interface

    type boundry_condition
        procedure(bc_func), pointer, nopass :: apply => null()
        procedure(pp_func), pointer, nopass :: postproc => null()
        integer :: places = 0
	    real(dp) :: speed = 0.0_dp
    end type


end module lbm_d2q9_types