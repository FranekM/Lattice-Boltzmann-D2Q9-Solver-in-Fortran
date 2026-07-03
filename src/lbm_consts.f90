!> Module for constants of D2Q9 Lattice Boltzman Method
!!
module lbm_d2q9_consts
    use iso_fortran_env, only: dp => real64

    implicit none

    integer, parameter :: Nd = 2            !! number of spatial dimensions
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


    integer, parameter :: W_WALL = 1        !! west wall mask for cummunication step
    integer, parameter :: S_WALL = 2        !! south wall mask for cummunication step
    integer, parameter :: E_WALL = 4        !! east wall mask for cummunication step
    integer, parameter :: N_WALL = 8        !! north wall mask for cummunication step
    integer, parameter :: ALL_WALLS = 15    !! all walls mask

    ! TODO: consider droping _corners. Do e.g. SW_CORNER = S_WALL + W_WALL instead
    integer, parameter :: SW_CORNER = 16    !! southwest corner mask for cummunication step
    integer, parameter :: SE_CORNER = 32    !! southeast corner mask for cummunication step
    integer, parameter :: NW_CORNER = 64    !! northwest corner mask for cummunication step
    integer, parameter :: NE_CORNER = 128   !! northeast corner mask for cummunication step
    integer, parameter :: ALL_CORNERS = 240 !! all corners mask

    integer, parameter :: ALL_PLACES = 255  !! all places (both walls and walls) mask

    integer, parameter :: MAX_BC_COUNT = 8  !! max number of boundry conditions (max 4 walls and 4 corners)

    real(dp), parameter :: PI = acos(-1.0d0)    !! pi appproximation for calcuations
    real(dp), parameter :: INIT_VAL = 999.999   !! initial lattice value (should be overriden)

    real(dp), parameter :: cs2 = 1.0_dp / 3.0_dp    !! speed of sound on the lattice

    real(dp), parameter :: weights(0:8) = [                     &! weights:
        4.0d0/9.0d0,                                            &! for the rest channel
        1.0d0/9.0d0, 1.0d0/9.0d0, 1.0d0/9.0d0, 1.0d0/9.0d0,     &! for the straight channels
        1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0, 1.0d0/36.0d0  &! for the diagonal channels
    ]                                                            !! D2Q9 weights for equilibirum function

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

contains

    function mask(xs, x) result(m)
        integer, intent(in) :: xs, x
        logical :: m
        m = iand(xs, x) /= 0
    end function

end module lbm_d2q9_consts