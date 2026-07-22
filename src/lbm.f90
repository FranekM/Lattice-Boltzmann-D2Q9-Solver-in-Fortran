program lbm_d2q9
	use iso_fortran_env, only: dp => real64
	implicit none

	! HYPER-PARAMETERS

    integer, parameter :: Nq = 9            !! number of velocity channels

    integer, parameter :: QR = 0            !! rest channel index
    integer, parameter :: QE = 1            !! east channel index
    integer, parameter :: QN = 2            !! north channel index
    integer, parameter :: QW = 3            !! west channel index
    integer, parameter :: QS = 4            !! sounth chnnel index
    integer, parameter :: QNE = 5           !! northeast channel index
    integer, parameter :: QNW = 6           !! northwest channel index
    integer, parameter :: QSW = 7           !! southwest channel index
    integer, parameter :: QSE = 8           !! southeast channel index

    integer, parameter :: CY = 0            !! y component index
    integer, parameter :: CX = 1            !! x component index

	real(dp), parameter :: cs2 = 1.0_dp / 3.0_dp    !! speed of sound on the lattice

	integer, parameter :: c(0:1, 0:8) = reshape([ &
         0,  0,     &! q=0: no movment
         0,  1,     &! q=1: move east
         1,  0,     &! q=2: move north
         0, -1,     &! q=3: move west
        -1,  0,     &! q=4: move south
         1,  1,     &! q=5: move NE
         1, -1,     &! q=6: move NW
        -1, -1,     &! q=7: move SW
        -1,  1      &! q=8: move SE
    ], [2,9])        !! 2DQ9 direction vectors (in column-major order -> y, x)

	real(dp), parameter :: c_dp(0:1, 0:8) = reshape([ &
         0.0_dp,  0.0_dp,     &! q=0: no movment
         0.0_dp,  1.0_dp,     &! q=1: move east
         1.0_dp,  0.0_dp,     &! q=2: move north
         0.0_dp, -1.0_dp,     &! q=3: move west
        -1.0_dp,  0.0_dp,     &! q=4: move south
         1.0_dp,  1.0_dp,     &! q=5: move NE
         1.0_dp, -1.0_dp,     &! q=6: move NW
        -1.0_dp, -1.0_dp,     &! q=7: move SW
        -1.0_dp,  1.0_dp      &! q=8: move SE
    ], [2,9])        !! 2DQ9 direction vectors (in column-major order -> y, x)

    real(dp), parameter :: weights(0:8) = [                     &! weights:
        4.0d0/9.0d0,                                            &! for the rest channel
        1.0d0/9.0d0, 1.0d0/9.0d0, 1.0d0/9.0d0, 1.0d0/9.0d0,     &! for the straight channels
        1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0  &! for the diagonal channels
    ]                                                            !! D2Q9 weights for equilibirum function


	! TIME MEASUREMENT VARIABLES

	real :: total_start_time            = 0
	real :: total_end_time              = 0

	real :: total_time                  = 0
	real :: total_collision_time        = 0
	real :: total_communication_time    = 0
	real :: total_copy_time             = 0
	real :: total_sync_time             = 0

	real :: st  = 0
	real :: et  = 0
	real :: st2 = 0
	real :: et2 = 0


	! SIMULATION PARAMETERS

	real(dp), parameter :: OMEGA            = 1.041_dp
	real(dp), parameter :: ONE_MINUS_OMEGA  = 1.0_dp - 1.041_dp
	real(dp), parameter :: INIT_RHO         = 1.0_dp
	real(dp), parameter :: SPEED            = 0.1_dp

	integer, parameter :: LATTICE_SIZE  = 1536
	integer, parameter :: TIME_STEPS    = 1000

	real(dp), parameter :: OFFSET_QN  = (2.0_dp * weights(QN ) * c_dp(CX, QN ) * SPEED / cs2)
	real(dp), parameter :: OFFSET_QNW = (2.0_dp * weights(QNW) * c_dp(CX, QNW) * SPEED / cs2)
	real(dp), parameter :: OFFSET_QNE = (2.0_dp * weights(QNE) * c_dp(CX, QNE) * SPEED / cs2)

	character(32), parameter :: LOG_PARENT   = "out/scaling_periodic"
	character(32), parameter :: LOG_FILENAME = "re1000-"
	character(32), parameter :: LOG_EXT      = ".tsv"

	character(32), parameter :: STAT_FILENAME = "out/scaling_periodic/stats"
	character(32), parameter :: STAT_EXT      = ".tsv"


	! VARIABLES

    integer :: W, S, E, N                             !! border layers indices
    integer :: gW, gS, gE, gN                         !! ghost layers indices
    integer :: Ny, Nx                                 !! lengths of y and x domain
    integer :: Iy, Ix                                 !! image y and x indices
    integer :: Niy, Nix                               !! number of image y and x indices

    real(dp), allocatable :: lattice(:, :, :)[:, :]   !! distribution function
    real(dp), allocatable :: sbuffer(:, :, :)         !! streaming buffer

	integer :: t, q, y, x

	! collision
    integer :: coll_x, coll_y, coll_q
    real(dp) :: rho, inv_rho, uy, ux
    real(dp) :: feq, au2, cu
    real(dp) :: fR, fE, fN, fW, fS, fNE, fNW, fSW, fSE

	! communication
	integer :: src_Ix, src_Iy, comm_x
    real(dp) :: rho_wall

	! logging
	integer :: log_unit, img
	integer :: log_y, log_x, gy, gx
    real(dp) :: mlups
	character(256) :: log_file

	! initial condition
	real(dp) :: eps =  0.01_dp
	real(dp) :: sigma = 10.0_dp
	real(dp) :: gaussian_rho = 0.0_dp


	! EXECUTION

	if (this_image() == 1)                                           &
		print '(A,A,I0,A)',  'LBM D2Q9 LID DRIVEN CAVITY ',          &
							 'starting on: ', num_images(), ' images'

	! lbm initialization
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

!	do x = 0, Nx
!		do y = 0, Ny
!			do q = 0, Nq-1
!				gaussian_rho = INIT_RHO + eps * exp(-((x - (Nx/2))**2 + (y - (Ny/2))**2) / (2.0_dp * sigma**2))
!				lattice(q, y, x) = weights(q) * gaussian_rho * (1.0_dp + 3.0_dp*cu + 4.5_dp*cu*cu)
!			end do
!		end do
!	end do

	! start simulation
	call cpu_time(total_start_time)
	do t = 0, TIME_STEPS

		! start collision
		call cpu_time(st)

		do coll_x = W, E
            do coll_y = S, N
                rho = 0.0_dp
                uy = 0.0_dp
                ux = 0.0_dp


                fR = lattice(QR, coll_y-c(CY, QR), coll_x-c(CX, QR))
                rho = rho + fR
                uy = uy + fR * c_dp(CY, QR)
                ux = ux + fR * c_dp(CX, QR)

                fE = lattice(QE, coll_y-c(CY, QE), coll_x-c(CX, QE))
                rho = rho + fE
                uy = uy + fE * c_dp(CY, QE)
                ux = ux + fE * c_dp(CX, QE)

                fN = lattice(QN, coll_y-c(CY, QN), coll_x-c(CX, QN))
                rho = rho + fN
                uy = uy + fN * c_dp(CY, QN)
                ux = ux + fN * c_dp(CX, QN)

                fW = lattice(QW, coll_y-c(CY, QW), coll_x-c(CX, QW))
                rho = rho + fW
                uy = uy + fW * c_dp(CY, QW)
                ux = ux + fW * c_dp(CX, QW)

                fS = lattice(QS, coll_y-c(CY, QS), coll_x-c(CX, QS))
                rho = rho + fS
                uy = uy + fS * c_dp(CY, QS)
                ux = ux + fS * c_dp(CX, QS)

                fNE = lattice(QNE, coll_y-c(CY, QNE), coll_x-c(CX, QNE))
                rho = rho + fNE
                uy = uy + fNE * c_dp(CY, QNE)
                ux = ux + fNE * c_dp(CX, QNE)

                fNW = lattice(QNW, coll_y-c(CY, QNW), coll_x-c(CX, QNW))
                rho = rho + fNW
                uy = uy + fNW * c_dp(CY, QNW)
                ux = ux + fNW * c_dp(CX, QNW)

                fSW = lattice(QSW, coll_y-c(CY, QSW), coll_x-c(CX, QSW))
                rho = rho + fSW
                uy = uy + fSW * c_dp(CY, QSW)
                ux = ux + fSW * c_dp(CX, QSW)

                fSE = lattice(QSE, coll_y-c(CY, QSE), coll_x-c(CX, QSE))
                rho = rho + fSE
                uy = uy + fSE * c_dp(CY, QSE)
                ux = ux + fSE * c_dp(CX, QSE)


                inv_rho = 1.0_dp / rho
                uy = uy * inv_rho
                ux = ux * inv_rho
                au2 = -1.5_dp * (uy*uy + ux*ux)


                cu = c_dp(CY, QR) * uy + c_dp(CX, QR) * ux
                feq = weights(QR) * rho * (1.0_dp + 3.0_dp*cu + 4.5_dp*cu*cu + au2)
                sbuffer(QR, coll_y, coll_x) = ONE_MINUS_OMEGA * fR + OMEGA * feq

                cu = c_dp(CY, QE) * uy + c_dp(CX, QE) * ux
                feq = weights(QE) * rho * (1.0_dp + 3.0_dp*cu + 4.5_dp*cu*cu + au2)
                sbuffer(QE, coll_y, coll_x) = ONE_MINUS_OMEGA * fE + OMEGA * feq

                cu = c_dp(CY, QN) * uy + c_dp(CX, QN) * ux
                feq = weights(QN) * rho * (1.0_dp + 3.0_dp*cu + 4.5_dp*cu*cu + au2)
                sbuffer(QN, coll_y, coll_x) = ONE_MINUS_OMEGA * fN + OMEGA * feq

                cu = c_dp(CY, QW) * uy + c_dp(CX, QW) * ux
                feq = weights(QW) * rho * (1.0_dp + 3.0_dp*cu + 4.5_dp*cu*cu + au2)
                sbuffer(QW, coll_y, coll_x) = ONE_MINUS_OMEGA * fW + OMEGA * feq

                cu = c_dp(CY, QS) * uy + c_dp(CX, QS) * ux
                feq = weights(QS) * rho * (1.0_dp + 3.0_dp*cu + 4.5_dp*cu*cu + au2)
                sbuffer(QS, coll_y, coll_x) = ONE_MINUS_OMEGA * fS + OMEGA * feq

                cu = c_dp(CY, QNE) * uy + c_dp(CX, QNE) * ux
                feq = weights(QNE) * rho * (1.0_dp + 3.0_dp*cu + 4.5_dp*cu*cu + au2)
                sbuffer(QNE, coll_y, coll_x) = ONE_MINUS_OMEGA * fNE + OMEGA * feq

                cu = c_dp(CY, QNW) * uy + c_dp(CX, QNW) * ux
                feq = weights(QNW) * rho * (1.0_dp + 3.0_dp*cu + 4.5_dp*cu*cu + au2)
                sbuffer(QNW, coll_y, coll_x) = ONE_MINUS_OMEGA * fNW + OMEGA * feq

                cu = c_dp(CY, QSW) * uy + c_dp(CX, QSW) * ux
                feq = weights(QSW) * rho * (1.0_dp + 3.0_dp*cu + 4.5_dp*cu*cu + au2)
                sbuffer(QSW, coll_y, coll_x) = ONE_MINUS_OMEGA * fSW + OMEGA * feq

                cu = c_dp(CY, QSE) * uy + c_dp(CX, QSE) * ux
                feq = weights(QSE) * rho * (1.0_dp + 3.0_dp*cu + 4.5_dp*cu*cu + au2)
                sbuffer(QSE, coll_y, coll_x) = ONE_MINUS_OMEGA * fSE + OMEGA * feq

            end do
        end do

		call cpu_time(st2)

        lattice = sbuffer

	    call cpu_time(et2)
	    total_copy_time = total_copy_time + (et2 - st2)

		call cpu_time(et)
		total_collision_time = total_collision_time + (et - st)

		! end collision
		! start synchronization

		call cpu_time(st)

		sync all

		call cpu_time(et)
		total_sync_time = total_sync_time + (et - st)

		! end synchronization
		! start communication

		call cpu_time(st)

	    ! W wall communication
		src_Iy = modulo(Iy + c(CY, QW) - 1, Niy) + 1
		src_Ix = modulo(Ix + c(CX, QW) - 1, Nix) + 1
        lattice(QE,  S:N,   gW) = lattice(QE,  S:N,   E)[src_Iy, src_Ix]
        lattice(QNE, S:N-1, gW) = lattice(QNE, S:N-1, E)[src_Iy, src_Ix]
        lattice(QSE, S+1:N, gW) = lattice(QSE, S+1:N, E)[src_Iy, src_Ix]

		! E wall communication
		src_Iy = modulo(Iy + c(CY, QE) - 1, Niy) + 1
		src_Ix = modulo(Ix + c(CX, QE) - 1, Nix) + 1
        lattice(QW,  S:N,   gE) = lattice(QW,  S:N,   W)[src_Iy, src_Ix]
        lattice(QNW, S:N-1, gE) = lattice(QNW, S:N-1, W)[src_Iy, src_Ix]
        lattice(QSW, S+1:N, gE) = lattice(QSW, S+1:N, W)[src_Iy, src_Ix]

		! S wall communication
		src_Iy = modulo(Iy + c(CY, QS) - 1, Niy) + 1
		src_Ix = modulo(Ix + c(CX, QS) - 1, Nix) + 1
        lattice(QN,  gS, W:E  ) = lattice(QN,  N, W:E  )[src_Iy, src_Ix]
        lattice(QNE, gS, W:E-1) = lattice(QNE, N, W:E-1)[src_Iy, src_Ix]
        lattice(QNW, gS, W+1:E) = lattice(QNW, N, W+1:E)[src_Iy, src_Ix]

		! N wall communication
		src_Iy = modulo(Iy + c(CY, QN) - 1, Niy) + 1
		src_Ix = modulo(Ix + c(CX, QN) - 1, Nix) + 1
        lattice(QS,  gN, W:E  ) = lattice(QS,  S, W:E  )[src_Iy, src_Ix]
        lattice(QSW, gN, W+1:E) = lattice(QSW, S, W+1:E)[src_Iy, src_Ix]
        lattice(QSE, gN, W:E-1) = lattice(QSE, S, W:E-1)[src_Iy, src_Ix]

		! NE corner communication
		src_Iy = modulo(Iy + c(CY, QNE) - 1, Niy) + 1
		src_Ix = modulo(Ix + c(CX, QNE) - 1, Nix) + 1
		lattice(QSW, gN, gE) = lattice(QSW, S, W)[src_Iy, src_Ix]

		! NW corner communication
		src_Iy = modulo(Iy + c(CY, QNW) - 1, Niy) + 1
		src_Ix = modulo(Ix + c(CX, QNW) - 1, Nix) + 1
        lattice(QSE, gN, gW) = lattice(QSE, S, E)[src_Iy, src_Ix]

		! SW corner communication
		src_Iy = modulo(Iy + c(CY, QSW) - 1, Niy) + 1
		src_Ix = modulo(Ix + c(CX, QSW) - 1, Nix) + 1
        lattice(QNE, gS, gW) = lattice(QNE, N, E)[src_Iy, src_Ix]

		! SE corner communication
		src_Iy = modulo(Iy + c(CY, QSE) - 1, Niy) + 1
		src_Ix = modulo(Ix + c(CX, QSE) - 1, Nix) + 1
		lattice(QNW, gS, gE) = lattice(QNW, N, W)[src_Iy, src_Ix]

		call cpu_time(et)
		total_communication_time = total_communication_time + (et - st)

		! end communication

	end do
	call cpu_time(total_end_time)
	total_time = (total_end_time - total_start_time)

	! end simulation
	! start log statistics

	write (log_file, '(A,A)'), trim(STAT_FILENAME), trim(STAT_EXT)

	do img = 1, num_images()
		sync all

		if (this_image() == img) then

			open(                        &
			    newunit  = log_unit,     &
			    file     = log_file,     &
				position = "append",     &
			    status   = "unknown",    &
			    action   = "write"       &
			)

			mlups = real(LATTICE_SIZE, dp) * real(LATTICE_SIZE, dp)
			mlups = mlups * real(TIME_STEPS, dp)
			mlups = mlups / 1000000_dp
			mlups = mlups / total_time

			write (log_unit, '(I0,A,I0,A,G,A,G,A,G,A,G,A,G,A,G,A,G,A,G)'),  &
				num_images(),               char(9),                        &
				this_image(),               char(9),                        &
				mlups,                      char(9),                        &
				total_time,                 char(9),                        &
				total_copy_time,            char(9),                        &
				total_sync_time,            char(9),                        &
				total_communication_time,   char(9),                        &
				total_collision_time

			close(log_unit)

		end if
	end do

	! end log statistics
	! start log results

	write (log_file, '(A,A,I0,A,A,I0,A)'),  &
	trim(LOG_PARENT),   "/",            &
	num_images(),       "/",            &
	trim(LOG_FILENAME), this_image(),   &
	trim(LOG_EXT)

	open(                       &
	    newunit = log_unit,     &
	    file    = log_file,     &
	    status  = "replace",    &
	    action  = "write"       &
	)

	do log_y = S, N
		do log_x = W, E
            rho = sum(lattice(:, log_y, log_x))
			uy = sum(lattice(:, log_y, log_x) * c(CY, :)) / rho
			ux = sum(lattice(:, log_y, log_x) * c(CX, :)) / rho

			gy = (Iy - 1) * Ny + log_y
			gx = (Ix - 1) * Nx + log_x
			write(log_unit, '(ES12.5,A,I0,A,I0,A,I0,A,ES20.10,A,ES20.10,A,ES20.10)') &
                            OMEGA,      char(9),        &
							TIME_STEPS, char(9),        &
							gy,         char(9),        &
							gx,         char(9),        &
							rho,        char(9),        &
							uy,         char(9),        &
							ux
		end do
	end do

	close(log_unit)

	! end log results

	sync all
    if (allocated(lattice)) deallocate(lattice)
    if (allocated(sbuffer)) deallocate(sbuffer)

	! end simulation

	sync all
	if (this_image() == 1)                                           &
		print '(A,A,I0,A)',  'LBM D2Q9 LID DRIVEN CAVITY ',          &
					         'finished on: ', num_images(), ' images'

end program lbm_d2q9