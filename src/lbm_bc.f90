module lbm_d2q9_boundry_conditions
	use lbm_d2q9_consts
	use lbm_d2q9_solver
	use lbm_d2q9_types

	implicit none
	private

	public :: periodic_bc
	public :: image_sync_bc
	public :: bounce_back_bc
	public :: moving_wall_correction

contains

	function offset(rho, speed, ci, q) result(off)
		integer, intent(in) :: ci, q
		real(dp), intent(in) :: rho, speed
		real(dp) :: off

		off = 2.0_dp * weights(q) * rho * c(ci, q) * speed / cs2
	end function

	function opposite_im_idx(ci) result(op)
		integer, intent(in) :: ci
		integer :: op

		if (ci == CX) then
			if (Ix == 1) op = Nix
			if (Ix == Nix) op = 1
			return
		end if

		if (ci == CY) then
			if (Iy == 1) op = Niy
			if (Iy == Niy) op = 1
			return
		end if

		op = -1
	end function


    subroutine periodic_bc(places)
		integer, intent(in) :: places
		integer :: src_Ix, src_Iy

        if (mask(places, W_WALL)) then
	        src_Iy = Iy
	        src_Ix = opposite_im_idx(CX)
	        lattice(QE,  S:N,   W) = lattice(QE,  S:N,   gE)[src_Iy, src_Ix]
            lattice(QNE, S+1:N, W) = lattice(QNE, S+1:N, gE)[src_Iy, src_Ix]
            lattice(QSE, S:N-1, W) = lattice(QSE, S:N-1, gE)[src_Iy, src_Ix]
        end if

        if (mask(places, E_WALL)) then
	        src_Iy = Iy
	        src_Ix = opposite_im_idx(CX)
	        lattice(QW,  S:N,   E) = lattice(QW,  S:N,   gW)[src_Iy, src_Ix]
	        lattice(QNW, S+1:N, E) = lattice(QNW, S+1:N, gW)[src_Iy, src_Ix]
	        lattice(QSW, S:N-1, E) = lattice(QSW, S:N-1, gW)[src_Iy, src_Ix]
        end if

        if (mask(places, S_WALL)) then
	        src_Iy = opposite_im_idx(CY)
	        src_Ix = Ix
	        lattice(QN,  S, W:E  ) = lattice(QN,  gN, W:E  )[src_Iy, src_Ix]
            lattice(QNE, S, W+1:E) = lattice(QNE, gN, W+1:E)[src_Iy, src_Ix]
            lattice(QNW, S, W:E-1) = lattice(QNW, gN, W:E-1)[src_Iy, src_Ix]
        end if

        if (mask(places, N_WALL)) then
	        src_Iy = opposite_im_idx(CY)
	        src_Ix = Ix
	        lattice(QS,  N, W:E  ) = lattice(QS,  gS, W:E  )[src_Iy, src_Ix]
            lattice(QSW, N, W:E-1) = lattice(QSW, gS, W:E-1)[src_Iy, src_Ix]
            lattice(QSE, N, W+1:E) = lattice(QSE, gS, W+1:E)[src_Iy, src_Ix]
		end if

        if (mask(places, NE_CORNER)) then
	        src_Iy = opposite_im_idx(CY)
	        src_Ix = opposite_im_idx(CX)
			lattice(QSW, N, E) = lattice(QSW, gS, gW)[src_Iy, src_Ix]
		end if

        if (mask(places, NW_CORNER)) then
	        src_Iy = opposite_im_idx(CY)
	        src_Ix = opposite_im_idx(CX)
	        lattice(QSE, N, W) = lattice(QSE, gS, gE)[src_Iy, src_Ix]
		end if

        if (mask(places, SW_CORNER)) then
	        src_Iy = opposite_im_idx(CY)
	        src_Ix = opposite_im_idx(CX)
	        lattice(QNE, S, W) = lattice(QNE, gN, gE)[src_Iy, src_Ix]
		end if

        if (mask(places, SE_CORNER)) then
	        src_Iy = opposite_im_idx(CY)
	        src_Ix = opposite_im_idx(CX)
			lattice(QNW, S, E) = lattice(QNW, gN, gW)[src_Iy, src_Ix]
		end if

    end subroutine

	subroutine bounce_back_bc(places)
		integer, intent(in) :: places

		if (mask(places, W_WALL)) then
			lattice(QE,  S:N, W) = lattice(QW,  S:N,    gW)
			lattice(QSE, S:N, W) = lattice(QNW, S+1:gN, gW)
			lattice(QNE, S:N, W) = lattice(QSW, gS:N-1, gW)
		end if

		if (mask(places, E_WALL)) then
			lattice(QW,  S:N, E) = lattice(QE,  S:N,    gE)
			lattice(QSW, S:N, E) = lattice(QNE, S+1:gN, gE)
			lattice(QNW, S:N, E) = lattice(QSE, gS:N-1, gE)
		end if

		if (mask(places, S_WALL)) then
			lattice(QN,  S, W:E) = lattice(QS,  gS, W:E   )
			lattice(QNE, S, W:E) = lattice(QSW, gS, gW:E-1)
			lattice(QNW, S, W:E) = lattice(QSE, gS, W+1:gE)
		end if

		if (mask(places, N_WALL)) then
			lattice(QS,  N, W:E) = lattice(QN,  gN, W:E   )
			lattice(QSE, N, W:E) = lattice(QNW, gN, gW:E-1)
			lattice(QSW, N, W:E) = lattice(QNE, gN, W+1:gE)
		end if

	end subroutine

	subroutine moving_wall_correction(places, speed)
		integer, intent(in) :: places
		real(dp), intent(in) :: speed
		real(dp) :: rho_w
		integer :: y, x

		if (mask(places, N_WALL)) then
	        do x = W, E
	            rho_w = measure_rho(N, x)
				lattice(QS,  N, x) = lattice(QS,  N, x) - offset(rho_w, speed, CX, QN)
	            lattice(QSE, N, x) = lattice(QSE, N, x) - offset(rho_w, speed, CX, QNW)
				lattice(QSW, N, x) = lattice(QSW, N, x) - offset(rho_w, speed, CX, QNE)
	        end do
		end if

		if (mask(places, S_WALL)) then
			do x = W, E
	            rho_w = measure_rho(S, x)
				lattice(QN,  S, x) = lattice(QN,  S, x) - offset(rho_w, speed, CX, QS)
	            lattice(QNE, S, x) = lattice(QNE, S, x) - offset(rho_w, speed, CX, QSW)
				lattice(QNW, S, x) = lattice(QNW, S, x) - offset(rho_w, speed, CX, QSE)
	        end do
		end if

		if (mask(places, W_WALL)) then
			do y = S, N
	            rho_w = measure_rho(y, W)
				lattice(QE,  y, W) = lattice(QW,  y, W) - offset(rho_w, speed, CY, QW)
	            lattice(QNE, y, W) = lattice(QNE, y, W) - offset(rho_w, speed, CY, QSW)
				lattice(QSE, y, W) = lattice(QSE, y, W) - offset(rho_w, speed, CY, QNW)
	        end do
		end if

		if (mask(places, E_WALL)) then
			do y = S, N
	            rho_w = measure_rho(y, E)
				lattice(QW,  y, E) = lattice(QW,  y, E) - offset(rho_w, speed, CY, QE)
	            lattice(QNW, y, E) = lattice(QNW, y, E) - offset(rho_w, speed, CY, QSE)
				lattice(QSW, y, E) = lattice(QSW, y, E) - offset(rho_w, speed, CY, QNE)
	        end do
		end if

	end subroutine

	subroutine image_sync_bc(places)
		integer, intent(in) :: places
		integer :: src_Ix, src_Iy

		if (mask(places, W_WALL)) then
			src_Iy = Iy + c(CY, QW)
			src_Ix = Ix + c(CX, QW)
	        lattice(QE,  S:N,   W) = lattice(QE,  S:N,   gE)[src_Iy, src_Ix]
            lattice(QNE, S+1:N, W) = lattice(QNE, S+1:N, gE)[src_Iy, src_Ix]
            lattice(QSE, S:N-1, W) = lattice(QSE, S:N-1, gE)[src_Iy, src_Ix]
        end if

        if (mask(places, E_WALL)) then
	        src_Iy = Iy + c(CY, QE)
			src_Ix = Ix + c(CX, QE)
	        lattice(QW,  S:N,   E) = lattice(QW,  S:N,   gW)[src_Iy, src_Ix]
	        lattice(QNW, S+1:N, E) = lattice(QNW, S+1:N, gW)[src_Iy, src_Ix]
	        lattice(QSW, S:N-1, E) = lattice(QSW, S:N-1, gW)[src_Iy, src_Ix]
        end if

        if (mask(places, S_WALL)) then
	        src_Iy = Iy + c(CY, QS)
			src_Ix = Ix + c(CX, QS)
	        lattice(QN,  S, W:E  ) = lattice(QN,  gN, W:E  )[src_Iy, src_Ix]
            lattice(QNE, S, W+1:E) = lattice(QNE, gN, W+1:E)[src_Iy, src_Ix]
            lattice(QNW, S, W:E-1) = lattice(QNW, gN, W:E-1)[src_Iy, src_Ix]
        end if

        if (mask(places, N_WALL)) then
	        src_Iy = Iy + c(CY, QN)
			src_Ix = Ix + c(CX, QN)
	        lattice(QS,  N, W:E  ) = lattice(QS,  gS, W:E  )[src_Iy, src_Ix]
            lattice(QSW, N, W:E-1) = lattice(QSW, gS, W:E-1)[src_Iy, src_Ix]
            lattice(QSE, N, W+1:E) = lattice(QSE, gS, W+1:E)[src_Iy, src_Ix]
		end if

        if (mask(places, NE_CORNER)) then
            src_Iy = Iy + c(CY, QNE)
			src_Ix = Ix + c(CX, QNE)
			lattice(QSW, N, E) = lattice(QSW, gS, gW)[src_Iy, src_Ix]
        end if

        if (mask(places, NW_CORNER)) then
	        src_Iy = Iy + c(CY, QNW)
			src_Ix = Ix + c(CX, QNW)
	        lattice(QSE, N, W) = lattice(QSE, gS, gE)[src_Iy, src_Ix]
		end if

        if (mask(places, SW_CORNER)) then
	        src_Iy = Iy + c(CY, QSW)
			src_Ix = Ix + c(CX, QSW)
	        lattice(QNE, S, W) = lattice(QNE, gN, gE)[src_Iy, src_Ix]
		end if

        if (mask(places, SE_CORNER)) then
	        src_Iy = Iy + c(CY, QSE)
			src_Ix = Ix + c(CX, QSE)
			lattice(QNW, S, E) = lattice(QNW, gN, gW)[src_Iy, src_Ix]
		end if

	end subroutine

end module lbm_d2q9_boundry_conditions