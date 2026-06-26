module lbm_d2q9_boundry_conditions
	use lbm_d2q9_consts
	use lbm_d2q9_solver
	use lbm_d2q9_types

	implicit none
	private

	public :: periodic_bc
	public :: bounce_back_bc
	public :: moving_wall_bc
	public :: flow_bc
	public :: make_bc

contains

	function mask(xs, x) result(m)
        integer, intent(in) :: xs, x
        logical :: m
        m = iand(xs, x) /= 0
    end function


	function offset(rho_w, u_w, ci, q) result(off)
		integer, intent(in) :: ci, q
		real(dp), intent(in) :: rho_w, u_w
		real(dp) :: off

		off = 2.0_dp * weights(q) * rho_w * c(ci, q) * u_w / cs2
	end function


	function make_bc(func, places, speed) result(bc)
		procedure(bc_func), intent(in), pointer :: func
		integer, intent(in) :: places
		real(dp), intent(in) :: speed
		type(boundry_condition) :: bc

		bc%apply => func
		bc%places = places
		bc%speed = speed
	end function


    subroutine periodic_bc(places, speed)
		integer, intent(in) :: places
		real(dp), intent(in) :: speed

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

    end subroutine


	subroutine bounce_back_bc(places, speed)
		integer, intent(in) :: places
		real(dp), intent(in) :: speed

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

	end subroutine


	subroutine moving_wall_bc(places, speed)
		integer, intent(in) :: places
		real(dp), intent(in) :: speed
		real(dp) :: rho_w
		integer :: y, x

		if (mask(places, W_WALL)) then
			do y = S, N
				rho_w = measure_rho(y, W)
				lattice(QE, y, W) = lattice(QW, y, gW) - offset(rho_w, speed, CY, QW)
				lattice(QNE, y, W) = lattice(QSW, y-1, gW) - offset(rho_w, speed, CY, QSW)
				lattice(QSE, y, W) = lattice(QNW, y+1, gW) - offset(rho_w, speed, CY, QNW)
			end do
		end if

		if (mask(places, E_WALL)) then
			do y = S, N
				rho_w = measure_rho(y, E)
				lattice(QW, y, E) = lattice(QE, y, gE) - offset(rho_w, speed, CY, QE)
				lattice(QNW, y, E) = lattice(QSE, y-1, gE) - offset(rho_w, speed, CY, QSE)
				lattice(QSW, y, E) = lattice(QNE, y+1, gE) - offset(rho_w, speed, CY, QNE)
			end do
		end if

		if (mask(places, S_WALL)) then
			do x = W, E
				rho_w = measure_rho(S, x)
				lattice(QN, S, x) = lattice(QS, gS, x) - offset(rho_w, speed, CX, QS)
				lattice(QNE, S, x) = lattice(QSW, gS, x-1) - offset(rho_w, speed, CX, QSW)
				lattice(QNW, S, x) = lattice(QSE, gS, x+1) - offset(rho_w, speed, CX, QSE)
			end do
		end if

		if (mask(places, N_WALL)) then
			do x = W, E
				rho_w = measure_rho(N, x)
				lattice(QS, N, x) = lattice(QN, gN, x) - offset(rho_w, speed, CX, QN)
				lattice(QSE, N, x) = lattice(QNW, gN, x-1) - offset(rho_w, speed, CX, QNW)
				lattice(QSW, N, x) = lattice(QNE, gN, x+1) - offset(rho_w, speed, CX, QNE)
			end do
		end if

	end subroutine


	subroutine flow_bc(places)
		integer, intent(in) :: places

		! TODO: communicate with neighbouring image
	end subroutine

end module lbm_d2q9_boundry_conditions