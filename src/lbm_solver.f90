module lbm_d2q9_solver
	use lbm_d2q9_global
	implicit none

	! VARIABLES

    integer :: W, S, E, N                             !! border layers indices
    integer :: gW, gS, gE, gN                         !! ghost layers indices
    integer :: Ny, Nx                                 !! lengths of y and x domain
    integer :: Iy, Ix                                 !! image y and x indices
    integer :: Niy, Nix                               !! number of image y and x indices

    real(dp), allocatable :: lattice(:, :, :)[:, :]   !! distribution function
    real(dp), allocatable :: sbuffer(:, :, :)         !! streaming buffer


contains


	subroutine init_lbm()
		integer :: q

        Niy = int(sqrt(real(num_images())))
        Nix = num_images() / Niy

        Ny = LATTICE_SIZE / Niy
        Nx = LATTICE_SIZE / Nix

        W = 1
        S = 1
        E = Nx
        N = Ny

        gW = W - 1
        gS = S - 1
        gE = E + 1
        gN = N + 1

        sync all
        allocate(lattice(0:Nq-1, gS:gN, gW:gE)[Niy, *])
        allocate(sbuffer(0:Nq-1, gS:gN, gW:gE))
		sync all

        Iy = this_image(lattice, 1)
        Ix = this_image(lattice, 2)

        do q = 0, Nq-1
			lattice(q, :, :) = weights(q) * INIT_RHO
        end do
    end subroutine

	subroutine end_lbm()
		sync all
        if (allocated(lattice)) deallocate(lattice)
        if (allocated(sbuffer)) deallocate(sbuffer)
	end subroutine


	subroutine collision_step()
        integer :: x, y, q
        real(dp) :: rho, inv_rho, uy, ux
        real(dp) :: feq, au2, cu
        real(dp) :: f(0:Nq-1)

        do x = W, E
            do y = S, N
                rho = 0.0_dp
                uy = 0.0_dp
                ux = 0.0_dp

                do q = 0, Nq-1
                    f(q) = lattice(q, y, x)
                    rho = rho + f(q)
                    uy = uy + f(q) * c_dp(CY, q)
                    ux = ux + f(q) * c_dp(CX, q)
                end do

                inv_rho = 1.0_dp / rho
                uy = uy * inv_rho
                ux = ux * inv_rho

                au2 = -1.5_dp * (uy*uy + ux*ux)

                do q = 0, Nq-1
                    cu = c_dp(CY, q) * uy + c_dp(CX, q) * ux
                    feq = weights(q) * rho * (1.0_dp + 3.0_dp*cu + 4.5_dp*cu*cu + au2)
                    lattice(q, y, x) = ONE_MINUS_OMEGA * f(q) + OMEGA * feq
                end do
            end do
        end do
    end subroutine

    subroutine streaming_step()
        integer :: y, x, q

        call cpu_time(st2)
        do x = W, E
	        do y = S, N
		        do q = 0, Nq-1
		            ! sbuffer(q, y+c(CY, q), x+c(CX, q)) = lattice(q, y, x)
			        sbuffer(q, y, x) = lattice(q, y-c(CY, q), x-c(CX, q))
		        end do
	        end do
        end do
        call cpu_time(et2)
        total_slide_time = total_slide_time + (et2 - st2)

        call cpu_time(st2)
        lattice = sbuffer
	    call cpu_time(et2)
	    total_copy_time = total_copy_time + (et2 - st2)
    end subroutine

    subroutine communication_step()

	    if (num_images() == 1) then
		    call communication_ONLY()

	    else if (Iy == 1 .and. Ix /= 1 .and. Ix /= Nix) then
		    call communication_S()

	    else if (Iy == Niy .and. Ix /= 1 .and. Ix /= Nix) then
		    call communication_N()

	    else if (Ix == 1 .and. Iy /= 1 .and. Iy /= Niy) then
		    call communication_W()

	    else if (Ix == Nix .and. Iy /= 1 .and. Iy /= Niy) then
		    call communication_E()

	    else if (Iy == 1 .and. Ix == 1) then
		    call communication_SW()

	    else if (Iy == Niy .and. Ix == 1) then
		    call communication_NW()

	    else if (Iy == Niy .and. Ix == Nix) then
		    call communication_NE()

		else if (Iy == 1 .and. Ix == Nix) then
			call communication_SE()

	    else
		    call communication_IN()

	    end if

    end subroutine


	subroutine communication_ONLY()
		real(dp) :: rho_wall
		integer :: x

		! N wall
		lattice(QS,  gN, W  :E  ) = lattice(QN,  N, W:E)
		lattice(QSE, gN, W-1:E-1) = lattice(QNW, N, W:E)
		lattice(QSW, gN, W+1:E+1) = lattice(QNE, N, W:E)

		! E wall
		lattice(QW,  S  :N,   gE) = lattice(QE,  S:N, E)
		lattice(QSW, S+1:N+1, gE) = lattice(QNE, S:N, E)
		lattice(QNW, S-1:N-1, gE) = lattice(QSE, S:N, E)

		! S wall
		lattice(QN,  gS, W  :E  ) = lattice(QS,  S, W:E)
		lattice(QNE, gS, W-1:E-1) = lattice(QSW, S, W:E)
		lattice(QNW, gS, W+1:E+1) = lattice(QSE, S, W:E)

		! W wall
		lattice(QE,  S  :N,   gW) = lattice(QW,  S:N, W)
		lattice(QSE, S+1:N+1, gW) = lattice(QNW, S:N, W)
		lattice(QNE, S-1:N-1, gW) = lattice(QSW, S:N, W)

		! moving wall correction
		do x = W, E
	        rho_wall = sum(lattice(:, gN, x))
			lattice(QS,  gN, x) = lattice(QS,  gN, x) - OFFSET_QN  * rho_wall
	        lattice(QSE, gN, x) = lattice(QSE, gN, x) - OFFSET_QNW * rho_wall
			lattice(QSW, gN, x) = lattice(QSW, gN, x) - OFFSET_QNE * rho_wall
	    end do

	end subroutine

	subroutine communication_N()
		integer :: src_Ix, src_Iy, x
		real(dp) :: rho_wall

		! N wall
		lattice(QS,  gN, W  :E  ) = lattice(QN,  N, W:E)
		lattice(QSE, gN, W-1:E-1) = lattice(QNW, N, W:E)
		lattice(QSW, gN, W+1:E+1) = lattice(QNE, N, W:E)

		! moving wall correction
		do x = W, E
	        rho_wall = sum(lattice(:, gN, x))
			lattice(QS,  gN, x) = lattice(QS,  gN, x) - OFFSET_QN  * rho_wall
	        lattice(QSE, gN, x) = lattice(QSE, gN, x) - OFFSET_QNW * rho_wall
			lattice(QSW, gN, x) = lattice(QSW, gN, x) - OFFSET_QNE * rho_wall
	    end do

		! W wall internal communication
		src_Iy = Iy + c(CY, QW)
		src_Ix = Ix + c(CX, QW)
        lattice(QE,  S:N,   gW) = lattice(QE,  S:N,   E)[src_Iy, src_Ix]
        lattice(QNE, S:N-1, gW) = lattice(QNE, S:N-1, E)[src_Iy, src_Ix]
        lattice(QSE, S+1:N, gW) = lattice(QSE, S+1:N, E)[src_Iy, src_Ix]

		! E wall internal communication
		src_Iy = Iy + c(CY, QE)
		src_Ix = Ix + c(CX, QE)
        lattice(QW,  S:N,   gE) = lattice(QW,  S:N,   W)[src_Iy, src_Ix]
        lattice(QNW, S:N-1, gE) = lattice(QNW, S:N-1, W)[src_Iy, src_Ix]
        lattice(QSW, S+1:N, gE) = lattice(QSW, S+1:N, W)[src_Iy, src_Ix]

		! S wall internal communication
		src_Iy = Iy + c(CY, QS)
		src_Ix = Ix + c(CX, QS)
        lattice(QN,  gS, W:E  ) = lattice(QN,  N, W:E  )[src_Iy, src_Ix]
        lattice(QNE, gS, W:E-1) = lattice(QNE, N, W:E-1)[src_Iy, src_Ix]
        lattice(QNW, gS, W+1:E) = lattice(QNW, N, W+1:E)[src_Iy, src_Ix]

		! SW corner internal communication
		src_Iy = Iy + c(CY, QSW)
		src_Ix = Ix + c(CX, QSW)
        lattice(QNE, gS, gW) = lattice(QNE, N, E)[src_Iy, src_Ix]

		! SE corner internal communication
		src_Iy = Iy + c(CY, QSE)
		src_Ix = Ix + c(CX, QSE)
		lattice(QNW, gS, gE) = lattice(QNW, N, W)[src_Iy, src_Ix]

	end subroutine

	subroutine communication_E()
		integer :: src_Ix, src_Iy

		! E wall
		lattice(QW,  S  :N,   gE) = lattice(QE,  S:N, E)
		lattice(QSW, S+1:N+1, gE) = lattice(QNE, S:N, E)
		lattice(QNW, S-1:N-1, gE) = lattice(QSE, S:N, E)

		! W wall internal communication
		src_Iy = Iy + c(CY, QW)
		src_Ix = Ix + c(CX, QW)
        lattice(QE,  S:N,   gW) = lattice(QE,  S:N,   E)[src_Iy, src_Ix]
        lattice(QNE, S:N-1, gW) = lattice(QNE, S:N-1, E)[src_Iy, src_Ix]
        lattice(QSE, S+1:N, gW) = lattice(QSE, S+1:N, E)[src_Iy, src_Ix]

		! S wall internal communication
		src_Iy = Iy + c(CY, QS)
		src_Ix = Ix + c(CX, QS)
        lattice(QN,  gS, W:E  ) = lattice(QN,  N, W:E  )[src_Iy, src_Ix]
        lattice(QNE, gS, W:E-1) = lattice(QNE, N, W:E-1)[src_Iy, src_Ix]
        lattice(QNW, gS, W+1:E) = lattice(QNW, N, W+1:E)[src_Iy, src_Ix]

		! N wall internal communication
		src_Iy = Iy + c(CY, QN)
		src_Ix = Ix + c(CX, QN)
        lattice(QS,  gN, W:E  ) = lattice(QS,  S, W:E  )[src_Iy, src_Ix]
        lattice(QSW, gN, W+1:E) = lattice(QSW, S, W+1:E)[src_Iy, src_Ix]
        lattice(QSE, gN, W:E-1) = lattice(QSE, S, W:E-1)[src_Iy, src_Ix]

		! NW corner internal communication
		src_Iy = Iy + c(CY, QNW)
		src_Ix = Ix + c(CX, QNW)
        lattice(QSE, gN, gW) = lattice(QSE, S, E)[src_Iy, src_Ix]

		! SW corner internal communication
		src_Iy = Iy + c(CY, QSW)
		src_Ix = Ix + c(CX, QSW)
        lattice(QNE, gS, gW) = lattice(QNE, N, E)[src_Iy, src_Ix]

	end subroutine

	subroutine communication_S()
		integer :: src_Ix, src_Iy

		! S wall
		lattice(QN,  gS, W  :E  ) = lattice(QS,  S, W:E)
		lattice(QNE, gS, W-1:E-1) = lattice(QSW, S, W:E)
		lattice(QNW, gS, W+1:E+1) = lattice(QSE, S, W:E)

		! W wall internal communication
		src_Iy = Iy + c(CY, QW)
		src_Ix = Ix + c(CX, QW)
        lattice(QE,  S:N,   gW) = lattice(QE,  S:N,   E)[src_Iy, src_Ix]
        lattice(QNE, S:N-1, gW) = lattice(QNE, S:N-1, E)[src_Iy, src_Ix]
        lattice(QSE, S+1:N, gW) = lattice(QSE, S+1:N, E)[src_Iy, src_Ix]

		! E wall internal communication
		src_Iy = Iy + c(CY, QE)
		src_Ix = Ix + c(CX, QE)
        lattice(QW,  S:N,   gE) = lattice(QW,  S:N,   W)[src_Iy, src_Ix]
        lattice(QNW, S:N-1, gE) = lattice(QNW, S:N-1, W)[src_Iy, src_Ix]
        lattice(QSW, S+1:N, gE) = lattice(QSW, S+1:N, W)[src_Iy, src_Ix]

		! N wall internal communication
		src_Iy = Iy + c(CY, QN)
		src_Ix = Ix + c(CX, QN)
        lattice(QS,  gN, W:E  ) = lattice(QS,  S, W:E  )[src_Iy, src_Ix]
        lattice(QSW, gN, W+1:E) = lattice(QSW, S, W+1:E)[src_Iy, src_Ix]
        lattice(QSE, gN, W:E-1) = lattice(QSE, S, W:E-1)[src_Iy, src_Ix]

		! NE corner internal communication
		src_Iy = Iy + c(CY, QNE)
		src_Ix = Ix + c(CX, QNE)
		lattice(QSW, gN, gE) = lattice(QSW, S, W)[src_Iy, src_Ix]

		! NW corner internal communication
		src_Iy = Iy + c(CY, QNW)
		src_Ix = Ix + c(CX, QNW)
        lattice(QSE, gN, gW) = lattice(QSE, S, E)[src_Iy, src_Ix]

	end subroutine

	subroutine communication_W()
		integer :: src_Ix, src_Iy

		! W wall
		lattice(QE,  S  :N,   gW) = lattice(QW,  S:N, W)
		lattice(QSE, S+1:N+1, gW) = lattice(QNW, S:N, W)
		lattice(QNE, S-1:N-1, gW) = lattice(QSW, S:N, W)

		! E wall internal communication
		src_Iy = Iy + c(CY, QE)
		src_Ix = Ix + c(CX, QE)
        lattice(QW,  S:N,   gE) = lattice(QW,  S:N,   W)[src_Iy, src_Ix]
        lattice(QNW, S:N-1, gE) = lattice(QNW, S:N-1, W)[src_Iy, src_Ix]
        lattice(QSW, S+1:N, gE) = lattice(QSW, S+1:N, W)[src_Iy, src_Ix]

		! S wall internal communication
		src_Iy = Iy + c(CY, QS)
		src_Ix = Ix + c(CX, QS)
        lattice(QN,  gS, W:E  ) = lattice(QN,  N, W:E  )[src_Iy, src_Ix]
        lattice(QNE, gS, W:E-1) = lattice(QNE, N, W:E-1)[src_Iy, src_Ix]
        lattice(QNW, gS, W+1:E) = lattice(QNW, N, W+1:E)[src_Iy, src_Ix]

		! N wall internal communication
		src_Iy = Iy + c(CY, QN)
		src_Ix = Ix + c(CX, QN)
        lattice(QS,  gN, W:E  ) = lattice(QS,  S, W:E  )[src_Iy, src_Ix]
        lattice(QSW, gN, W+1:E) = lattice(QSW, S, W+1:E)[src_Iy, src_Ix]
        lattice(QSE, gN, W:E-1) = lattice(QSE, S, W:E-1)[src_Iy, src_Ix]

		! NE corner internal communication
		src_Iy = Iy + c(CY, QNE)
		src_Ix = Ix + c(CX, QNE)
		lattice(QSW, gN, gE) = lattice(QSW, S, W)[src_Iy, src_Ix]

		! SE corner internal communication
		src_Iy = Iy + c(CY, QSE)
		src_Ix = Ix + c(CX, QSE)
		lattice(QNW, gS, gE) = lattice(QNW, N, W)[src_Iy, src_Ix]

	end subroutine

	subroutine communication_NE()
		integer :: src_Ix, src_Iy, x
		real(dp) :: rho_wall

		! N wall
		lattice(QS,  gN, W  :E  ) = lattice(QN,  N, W:E)
		lattice(QSE, gN, W-1:E-1) = lattice(QNW, N, W:E)
		lattice(QSW, gN, W+1:E+1) = lattice(QNE, N, W:E)

		! E wall
		lattice(QW,  S  :N,   gE) = lattice(QE,  S:N, E)
		lattice(QSW, S+1:N+1, gE) = lattice(QNE, S:N, E)
		lattice(QNW, S-1:N-1, gE) = lattice(QSE, S:N, E)

		! moving wall correction
		do x = W, E
	        rho_wall = sum(lattice(:, gN, x))
			lattice(QS,  gN, x) = lattice(QS,  gN, x) - OFFSET_QN  * rho_wall
	        lattice(QSE, gN, x) = lattice(QSE, gN, x) - OFFSET_QNW * rho_wall
			lattice(QSW, gN, x) = lattice(QSW, gN, x) - OFFSET_QNE * rho_wall
	    end do

		! W wall internal communication
		src_Iy = Iy + c(CY, QW)
		src_Ix = Ix + c(CX, QW)
        lattice(QE,  S:N,   gW) = lattice(QE,  S:N,   E)[src_Iy, src_Ix]
        lattice(QNE, S:N-1, gW) = lattice(QNE, S:N-1, E)[src_Iy, src_Ix]
        lattice(QSE, S+1:N, gW) = lattice(QSE, S+1:N, E)[src_Iy, src_Ix]

		! S wall internal communication
		src_Iy = Iy + c(CY, QS)
		src_Ix = Ix + c(CX, QS)
        lattice(QN,  gS, W:E  ) = lattice(QN,  N, W:E  )[src_Iy, src_Ix]
        lattice(QNE, gS, W:E-1) = lattice(QNE, N, W:E-1)[src_Iy, src_Ix]
        lattice(QNW, gS, W+1:E) = lattice(QNW, N, W+1:E)[src_Iy, src_Ix]

		! SW corner internal communication
		src_Iy = Iy + c(CY, QSW)
		src_Ix = Ix + c(CX, QSW)
        lattice(QNE, gS, gW) = lattice(QNE, N, E)[src_Iy, src_Ix]

	end subroutine

	subroutine communication_NW()
		integer :: src_Ix, src_Iy, x
		real(dp) :: rho_wall

		! N wall
		lattice(QS,  gN, W  :E  ) = lattice(QN,  N, W:E)
		lattice(QSE, gN, W-1:E-1) = lattice(QNW, N, W:E)
		lattice(QSW, gN, W+1:E+1) = lattice(QNE, N, W:E)

		! W wall
		lattice(QE,  S  :N,   gW) = lattice(QW,  S:N, W)
		lattice(QSE, S+1:N+1, gW) = lattice(QNW, S:N, W)
		lattice(QNE, S-1:N-1, gW) = lattice(QSW, S:N, W)

		! moving wall correction
		do x = W, E
	        rho_wall = sum(lattice(:, gN, x))
			lattice(QS,  gN, x) = lattice(QS,  gN, x) - OFFSET_QN  * rho_wall
	        lattice(QSE, gN, x) = lattice(QSE, gN, x) - OFFSET_QNW * rho_wall
			lattice(QSW, gN, x) = lattice(QSW, gN, x) - OFFSET_QNE * rho_wall
	    end do

		! E wall internal communication
		src_Iy = Iy + c(CY, QE)
		src_Ix = Ix + c(CX, QE)
        lattice(QW,  S:N,   gE) = lattice(QW,  S:N,   W)[src_Iy, src_Ix]
        lattice(QNW, S:N-1, gE) = lattice(QNW, S:N-1, W)[src_Iy, src_Ix]
        lattice(QSW, S+1:N, gE) = lattice(QSW, S+1:N, W)[src_Iy, src_Ix]

		! S wall internal communication
		src_Iy = Iy + c(CY, QS)
		src_Ix = Ix + c(CX, QS)
        lattice(QN,  gS, W:E  ) = lattice(QN,  N, W:E  )[src_Iy, src_Ix]
        lattice(QNE, gS, W:E-1) = lattice(QNE, N, W:E-1)[src_Iy, src_Ix]
        lattice(QNW, gS, W+1:E) = lattice(QNW, N, W+1:E)[src_Iy, src_Ix]

		! SE corner internal communication
		src_Iy = Iy + c(CY, QSE)
		src_Ix = Ix + c(CX, QSE)
		lattice(QNW, gS, gE) = lattice(QNW, N, W)[src_Iy, src_Ix]

	end subroutine

	subroutine communication_SE()
		integer :: src_Ix, src_Iy

		! S wall
		lattice(QN,  gS, W  :E  ) = lattice(QS,  S, W:E)
		lattice(QNE, gS, W-1:E-1) = lattice(QSW, S, W:E)
		lattice(QNW, gS, W+1:E+1) = lattice(QSE, S, W:E)

		! E wall
		lattice(QW,  S  :N,   gE) = lattice(QE,  S:N, E)
		lattice(QSW, S+1:N+1, gE) = lattice(QNE, S:N, E)
		lattice(QNW, S-1:N-1, gE) = lattice(QSE, S:N, E)

		! W wall internal communication
		src_Iy = Iy + c(CY, QW)
		src_Ix = Ix + c(CX, QW)
        lattice(QE,  S:N,   gW) = lattice(QE,  S:N,   E)[src_Iy, src_Ix]
        lattice(QNE, S:N-1, gW) = lattice(QNE, S:N-1, E)[src_Iy, src_Ix]
        lattice(QSE, S+1:N, gW) = lattice(QSE, S+1:N, E)[src_Iy, src_Ix]

		! N wall internal communication
		src_Iy = Iy + c(CY, QN)
		src_Ix = Ix + c(CX, QN)
        lattice(QS,  gN, W:E  ) = lattice(QS,  S, W:E  )[src_Iy, src_Ix]
        lattice(QSW, gN, W+1:E) = lattice(QSW, S, W+1:E)[src_Iy, src_Ix]
        lattice(QSE, gN, W:E-1) = lattice(QSE, S, W:E-1)[src_Iy, src_Ix]

		! NW corner internal communication
		src_Iy = Iy + c(CY, QNW)
		src_Ix = Ix + c(CX, QNW)
        lattice(QSE, gN, gW) = lattice(QSE, S, E)[src_Iy, src_Ix]

	end subroutine

	subroutine communication_SW()
		integer :: src_Ix, src_Iy

		! S wall
		lattice(QN,  gS, W  :E  ) = lattice(QS,  S, W:E)
		lattice(QNE, gS, W-1:E-1) = lattice(QSW, S, W:E)
		lattice(QNW, gS, W+1:E+1) = lattice(QSE, S, W:E)

		! W wall
		lattice(QE,  S  :N,   gW) = lattice(QW,  S:N, W)
		lattice(QSE, S+1:N+1, gW) = lattice(QNW, S:N, W)
		lattice(QNE, S-1:N-1, gW) = lattice(QSW, S:N, W)

		! E wall internal communication
		src_Iy = Iy + c(CY, QE)
		src_Ix = Ix + c(CX, QE)
        lattice(QW,  S:N,   gE) = lattice(QW,  S:N,   W)[src_Iy, src_Ix]
        lattice(QNW, S:N-1, gE) = lattice(QNW, S:N-1, W)[src_Iy, src_Ix]
        lattice(QSW, S+1:N, gE) = lattice(QSW, S+1:N, W)[src_Iy, src_Ix]

		! N wall internal communication
		src_Iy = Iy + c(CY, QN)
		src_Ix = Ix + c(CX, QN)
        lattice(QS,  gN, W:E  ) = lattice(QS,  S, W:E  )[src_Iy, src_Ix]
        lattice(QSW, gN, W+1:E) = lattice(QSW, S, W+1:E)[src_Iy, src_Ix]
        lattice(QSE, gN, W:E-1) = lattice(QSE, S, W:E-1)[src_Iy, src_Ix]

		! NE corner internal communication
		src_Iy = Iy + c(CY, QNE)
		src_Ix = Ix + c(CX, QNE)
		lattice(QSW, gN, gE) = lattice(QSW, S, W)[src_Iy, src_Ix]

	end subroutine

	subroutine communication_IN()
		integer :: src_Ix, src_Iy

		! W wall internal communication
		src_Iy = Iy + c(CY, QW)
		src_Ix = Ix + c(CX, QW)
        lattice(QE,  S:N,   gW) = lattice(QE,  S:N,   E)[src_Iy, src_Ix]
        lattice(QNE, S:N-1, gW) = lattice(QNE, S:N-1, E)[src_Iy, src_Ix]
        lattice(QSE, S+1:N, gW) = lattice(QSE, S+1:N, E)[src_Iy, src_Ix]

		! E wall internal communication
		src_Iy = Iy + c(CY, QE)
		src_Ix = Ix + c(CX, QE)
        lattice(QW,  S:N,   gE) = lattice(QW,  S:N,   W)[src_Iy, src_Ix]
        lattice(QNW, S:N-1, gE) = lattice(QNW, S:N-1, W)[src_Iy, src_Ix]
        lattice(QSW, S+1:N, gE) = lattice(QSW, S+1:N, W)[src_Iy, src_Ix]

		! S wall internal communication
		src_Iy = Iy + c(CY, QS)
		src_Ix = Ix + c(CX, QS)
        lattice(QN,  gS, W:E  ) = lattice(QN,  N, W:E  )[src_Iy, src_Ix]
        lattice(QNE, gS, W:E-1) = lattice(QNE, N, W:E-1)[src_Iy, src_Ix]
        lattice(QNW, gS, W+1:E) = lattice(QNW, N, W+1:E)[src_Iy, src_Ix]

		! N wall internal communication
		src_Iy = Iy + c(CY, QN)
		src_Ix = Ix + c(CX, QN)
        lattice(QS,  gN, W:E  ) = lattice(QS,  S, W:E  )[src_Iy, src_Ix]
        lattice(QSW, gN, W+1:E) = lattice(QSW, S, W+1:E)[src_Iy, src_Ix]
        lattice(QSE, gN, W:E-1) = lattice(QSE, S, W:E-1)[src_Iy, src_Ix]

		! NE corner internal communication
		src_Iy = Iy + c(CY, QNE)
		src_Ix = Ix + c(CX, QNE)
		lattice(QSW, gN, gE) = lattice(QSW, S, W)[src_Iy, src_Ix]

		! NW corner internal communication
		src_Iy = Iy + c(CY, QNW)
		src_Ix = Ix + c(CX, QNW)
        lattice(QSE, gN, gW) = lattice(QSE, S, E)[src_Iy, src_Ix]

		! SW corner internal communication
		src_Iy = Iy + c(CY, QSW)
		src_Ix = Ix + c(CX, QSW)
        lattice(QNE, gS, gW) = lattice(QNE, N, E)[src_Iy, src_Ix]

		! SE corner internal communication
		src_Iy = Iy + c(CY, QSE)
		src_Ix = Ix + c(CX, QSE)
		lattice(QNW, gS, gE) = lattice(QNW, N, W)[src_Iy, src_Ix]

	end subroutine

end module lbm_d2q9_solver