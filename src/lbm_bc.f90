module lbm_d2q9_boundry_conditions
	use lbm_d2q9_consts
	use lbm_d2q9_solver
	use lbm_d2q9_types

	implicit none
	private

	public :: periodic_bc
	public :: bounce_back_bc
	public :: flow_bc
	public :: make_bc

contains

	function mask(xs, x) result(m)
        integer, intent(in) :: xs, x
        logical :: m
        m = iand(xs, x) /= 0
    end function mask


	function make_bc(func, places) result(bc)
		procedure(bc_func), intent(in), pointer :: func
		integer, intent(in) :: places
		type(boundry_condition) :: bc

		bc%apply => func
		bc%places = places
	end function make_bc


    subroutine periodic_bc(places)
        integer, intent(in) :: places

!        print *, "E", mask(places, E_WALL)
!        print *, "N", mask(places, N_WALL)
!        print *, "W", mask(places, W_WALL)
!        print *, "S", mask(places, S_WALL)
!        print *, "NE", mask(places, NE_CORNER)
!        print *, "NW", mask(places, NW_CORNER)
!        print *, "SW", mask(places, SW_CORNER)
!        print *, "SE", mask(places, SE_CORNER)

        if (mask(places, W_WALL)) then
	        lattice(QE,  S:N, W) = lattice(QE,  S:N, gE)
            lattice(QNE, S:N, W) = lattice(QNE, S:N, gE)
            lattice(QSE, S:N, W) = lattice(QSE, S:N, gE)
        end if

        if (mask(places, E_WALL)) then
	        lattice(QW,  S:N, E) = lattice(QW,  S:N, gW)
	        lattice(QNW, S:N, E) = lattice(QNW, S:N, gW)
	        lattice(QSW, S:N, E) = lattice(QSW, S:N, gW)
        end if

        if (mask(places, S_WALL)) then
	        lattice(QN,  S, W:E) = lattice(QN,  gN, W:E)
            lattice(QNE, S, W:E) = lattice(QNE, gN, W:E)
            lattice(QNW, S, W:E) = lattice(QNW, gN, W:E)
        end if

        if (mask(places, N_WALL)) then
	        lattice(QS,  N, W:E) = lattice(QS,  gS, W:E)
            lattice(QSW, N, W:E) = lattice(QSW, gS, W:E)
            lattice(QSE, N, W:E) = lattice(QSE, gS, W:E)
		end if

        if (mask(places, NE_CORNER)) &
			lattice(QNE, S, W) = lattice(QNE, gN, gE)

        if (mask(places, NW_CORNER)) &
	        lattice(QNW, S, E) = lattice(QNW, gN, gW)

        if (mask(places, SW_CORNER)) &
	        lattice(QSW, N, E) = lattice(QSW, gS, gW)

        if (mask(places, SE_CORNER)) &
			lattice(QSE, N, W) = lattice(QSE, gS, gE)

    end subroutine periodic_bc


	subroutine bounce_back_bc(places)
		integer, intent(in) :: places

		if (mask(places, W_WALL)) then
			lattice(QE,  S:N,  W) = lattice(QW,  S:N,  gW)
			lattice(QSE, gS:N, W) = lattice(QNW, S:gN, gW)
			lattice(QNE, S:gN, W) = lattice(QSW, gS:N, gW)
		end if

		if (mask(places, E_WALL)) then
			lattice(QW, S:N, E) = lattice(QE, S:N, gE)
			lattice(QSW, gS:N, E) = lattice(QNE, S:gN, gE)
			lattice(QNW, S:gN, E) = lattice(QSE, gS:N, gE)
		end if

		if (mask(places, S_WALL)) then
			lattice(QN, S, W:E) = lattice(QS, gS, W:E)
			lattice(QNE, S, W:gE) = lattice(QSW, gS, gW:E)
			lattice(QNW, S, gW:E) = lattice(QSE, gS, W:gE)
		end if

		if (mask(places, N_WALL)) then
			lattice(QS, N, W:E) = lattice(QN, gN, W:E)
			lattice(QSE, N, W:gE) = lattice(QNW, gN, gW:E)
			lattice(QSW, N, gW:E) = lattice(QNE, gN, W:gE)
		end if

!        if (mask(places, NE_CORNER)) &
!			lattice(QSW, N, E) = lattice(QNE, gN, gE)
!
!        if (mask(places, NW_CORNER)) &
!			lattice(QSE, N, W) = lattice(QNW, gN, gW)
!
!		if (mask(places, SW_CORNER)) &
!	        lattice(QNE, S, W) = lattice(QSW, gS, gW)
!
!		if (mask(places, SE_CORNER)) &
!	        lattice(QNW, S, E) = lattice(QSE, gS, gE)

	end subroutine bounce_back_bc


	subroutine moving_wall_bc_places

	end subroutine


	subroutine flow_bc(places)
		integer, intent(in) :: places

		! TODO: communicate with neighbouring image
	end subroutine flow_bc

end module lbm_d2q9_boundry_conditions