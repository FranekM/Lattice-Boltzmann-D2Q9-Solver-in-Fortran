module lbm_d2q9_fast_global
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
	real :: total_streaming_time        = 0
	real :: total_communication_time    = 0
	real :: total_slide_time            = 0
	real :: total_copy_time             = 0

	real :: st  = 0
	real :: et  = 0
	real :: st2 = 0
	real :: et2 = 0


	! SIMULATION PARAMETERS

	real(dp), parameter :: OMEGA = 1.781_dp
	real(dp), parameter :: ONE_MINUS_OMEGA = 1.0_dp - 1.781_dp
	real(dp), parameter :: INIT_RHO = 1.0_dp
	real(dp), parameter :: SPEED = 0.1_dp

	integer, parameter :: LATTICE_SIZE = 1024
	integer, parameter :: TIME_STEPS = 1000000

	real(dp), parameter :: OFFSET_QN  = (2.0_dp * weights(QN ) * c_dp(CX, QN ) * SPEED / cs2)
	real(dp), parameter :: OFFSET_QNW = (2.0_dp * weights(QNW) * c_dp(CX, QNW) * SPEED / cs2)
	real(dp), parameter :: OFFSET_QNE = (2.0_dp * weights(QNE) * c_dp(CX, QNE) * SPEED / cs2)

	character(32), parameter :: LOG_PARENT   = "out/fast"
	character(32), parameter :: LOG_FILENAME = "re5000-"
	character(32), parameter :: LOG_EXT      = ".tsv"

	character(32), parameter :: STAT_FILENAME = "out/fast/stats"
	character(32), parameter :: STAT_EXT      = ".tsv"

end module lbm_d2q9_fast_global