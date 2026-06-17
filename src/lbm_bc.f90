module lbm_d2q9_boundry_conditions
	use lbm_d2q9_consts
	use lbm_d2q9_solver
	implicit none
	private

	public :: periodic_bc
	public :: bounce_back_bc
	public :: flow_bc

contains

	logical function mask(xs, x)
        integer, intent(in) :: xs, x
        mask = iand(xs, x) /= 0
    end function mask


    subroutine periodic_bc(walls, corners)
        integer, intent(in) :: walls, corners

        ! Copying straight channels (ghost walls -> border walls)
        ! Example (4th channel, copy wall gS -> N):
        ! g g g g g -> g g g g g
        ! g n n n g -> g * * * g
        ! g n n n g -> g n n n g
        ! g * * * g -> g g g g g
        ! g = ghost, n = normal, * = copied value
        if (mask(walls, W_WALL)) lattice(QE, S:N, W) = lattice(QE, S:N, gE) ! 1st q (E): gohst w. E -> border w. W
        if (mask(walls, S_WALL)) lattice(QN, S, W:E) = lattice(QN, gN, W:E) ! 2nd q (N): gohst w. N -> border w. S
        if (mask(walls, E_WALL)) lattice(QW, S:N, E) = lattice(QW, S:N, gW) ! 3rd q (W): gohst w. W -> border w. E
        if (mask(walls, N_WALL)) lattice(QS, N, W:E) = lattice(QS, gS, W:E) ! 4th q (S): gohst w. S -> border w. N

        ! Copying diagonal channels (ghost walls -> border walls)
        ! Example (5th channel, copy wall gE -> W and then gN -> S):
        ! g g g g g -> g g g g g -> g * * * g -> g g g g g
        ! g n n n * -> g * n n g -> g n n n g -> g n n n g
        ! g n n n * -> g * n n g -> g n n n g -> g * * * g
        ! g g g g g -> g g g g g -> g g g g g -> g g g g g
        ! g = ghost, n = normal, * = copied value
        if (mask(walls, W_WALL)) lattice(QNE, S:N, W) = lattice(QNE, S:N, gE) ! 5th q (NE): ghost w. E -> border w. W
        if (mask(walls, S_WALL)) lattice(QNE, S, W:E) = lattice(QNE, gN, W:E) ! 5th q (NE): ghost w. N -> border w. S
        if (mask(walls, E_WALL)) lattice(QNW, S:N, E) = lattice(QNW, S:N, gW) ! 6th q (NW): ghost w. W -> border w. E
        if (mask(walls, S_WALL)) lattice(QNW, S, W:E) = lattice(QNW, gN, W:E) ! 6th q (NW): ghost w. N -> border w. S
        if (mask(walls, E_WALL)) lattice(QSW, S:N, E) = lattice(QSW, S:N, gW) ! 7th q (SW): ghost w. W -> border w. E
        if (mask(walls, N_WALL)) lattice(QSW, N, W:E) = lattice(QSW, gS, W:E) ! 7th q (SW): ghost w. S -> border w. N
        if (mask(walls, W_WALL)) lattice(QSE, S:N, W) = lattice(QSE, S:N, gE) ! 8th q (SE): ghost w. E -> border w. W
        if (mask(walls, N_WALL)) lattice(QSE, N, W:E) = lattice(QSE, gS, W:E) ! 8th q (SE): ghost w. S -> border w. N

        ! communicating diagnal channels on corners (ghost corners -> border corners)
        ! Diagonal channel corner exmaple (8th channel, copy corner gSE -> NW):
        ! g g g g g -> g g g g g
        ! g n n n g -> g * n n g
        ! g n n n g -> g n n n g
        ! g g g g * -> g g g g g
        ! g = ghost, n = normal, * = copied value
        if (mask(corners, NE_CORNER)) lattice(QNE, S, W) = lattice(QNE, gN, gE) ! 5th q (NE): ghost c. NE -> border c. SW
        if (mask(corners, NW_CORNER)) lattice(QNW, S, E) = lattice(QNW, gN, gW) ! 6th q (NW): ghost c. NW -> border c. SE
        if (mask(corners, SW_CORNER)) lattice(QSW, N, E) = lattice(QSW, gS, gW) ! 7th q (SW): ghost c. SW -> border c. NE
        if (mask(corners, SE_CORNER)) lattice(QSE, N, W) = lattice(QSE, gS, gE) ! 8th q (SE): ghost c. SE -> border c. NW
    end subroutine periodic_bc


	subroutine bounce_back_bc(walls, corners)
		integer, intent(in) :: walls, corners

		if (mask(walls, W_WALL)) then
			lattice(QE, S:N, W) = lattice(QW, S:N, gW)
			lattice(QSE, S:N, W) = lattice(QNW, S:N, gW)
			lattice(QNE, S:N, W) = lattice(QSW, S:N, gW)
		end if

		if (mask(walls, E_WALL)) then
			lattice(QW, S:N, E) = lattice(QE, S:N, gE)
			lattice(QSW, S:N, E) = lattice(QNE, S:N, gE)
			lattice(QNW, S:N, E) = lattice(QSE, S:N, gE)
		end if

		if (mask(walls, S_WALL)) then
			lattice(QN, S, W:E) = lattice(QS, gS, W:E)
			lattice(QNE, S, W:E) = lattice(QSW, gS, W:E)
			lattice(QNW, S, W:E) = lattice(QSE, gS, W:E)
		end if

		if (mask(walls, N_WALL)) then
			lattice(QS, N, W:E) = lattice(QN, gN, W:E)
			lattice(QSE, N, W:E) = lattice(QNW, gN, W:E)
			lattice(QSW, N, W:E) = lattice(QNE, gN, W:E)
		end if

        if (mask(corners, NE_CORNER)) &
			lattice(QSW, N, E) = lattice(QNE, gN, gE)

        if (mask(corners, NW_CORNER)) &
			lattice(QSE, N, W) = lattice(QNW, gN, gW)

		if (mask(corners, SW_CORNER)) &
	        lattice(QNE, S, W) = lattice(QSW, gS, gW)

		if (mask(corners, SE_CORNER)) &
	        lattice(QNW, S, E) = lattice(QSE, gS, gE)
	end subroutine bounce_back_bc


	subroutine flow_bc(walls, corners)
		integer, intent(in) :: walls, corners

		! TODO: communicate with neighbouring image
	end subroutine flow_bc

end module lbm_d2q9_boundry_conditions