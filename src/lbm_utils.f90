module lbm_utils
    use iso_fortran_env, only: output_unit
    use lbm_d2q9_boundry_conditions
    use lbm_d2q9_consts
    use lbm_d2q9_solver
    use lbm_d2q9_types

    implicit none
    private

    public :: debug_print
    public :: prepare_output_file
    public :: test_mass_conservation
    public :: test_streaming_pbc

    real(dp), parameter :: EPS = 1.0e-10

    integer, parameter :: q_map(0:2, 0:2) = reshape([ &
            QNW, QN, QNE,   &! array that   0  1  2      6  2  5
            QW,  QR, QE,    &! helps to     3  4  5  ->  3  0  1
            QSW, QS, QSE    &! rearange q   6  7  8      7  4  8
    ], [3,3])

contains

    logical function prepare_output_file(filename)
		character(*), intent(in) :: filename

		logical :: exists
		character(10) :: answer
		integer :: unit

		inquire(file=filename, exist=exists)

		if (.not. exists) then
			prepare_output_file = .true.
			return
		end if

		print *, ' File already exists: ', trim(filename)
		print *, ' Override? [Y/n]'
		read(*, *) answer

		if (answer == 'Y') then
			open(newunit=unit, file=filename, status="old")
			close(unit, status="delete")
			prepare_output_file = .true.
		else
			prepare_output_file = .false.
		end if

	end function prepare_output_file


    subroutine debug_print()
        real(dp) :: output(1:(N-S+1)*3, 1:(E-W+1)*3)
        integer :: lat_x, lat_y     ! lattice x, y indices
        integer :: out_x, out_y     ! output  x, y indices
        integer :: qx, qy           ! for mapping q indices

        ! arrange values from lattice into output array
        do out_x = lbound(output, 2), ubound(output, 2)
            qx = mod(out_x-1, 3)
            lat_x = W + (out_x - 1) / 3
            do out_y = lbound(output, 1), ubound(output, 1)
                qy = mod(out_y-1, 3)
                lat_y = S + (out_y - 1) / 3
                output(out_y, out_x) = lattice(q_map(qy, qx), lat_y, lat_x)
            end do
        end do

        ! print output array row-wise
        do out_y = ubound(output, 1), lbound(output, 1), -1
            print '(*(F7.2,1X))', output(out_y, :)
        end do
    end subroutine debug_print


    logical function test_streaming_pbc(Nt) result(r)
        integer, intent(in) :: Nt   ! test max time steps
        real(dp) old_mass           ! mass before streaming step
        integer :: x, y, t, q       ! loop indicies

        ! requires that the lbm is already initilized
        ! initial condition: zero whole lattice
        ! set unique indetifiers only for channels at position (S, W)
        lattice = 0.0_dp
        sbuffer = 0.0_dp
        do q = 0, Nq-1
            lattice(q, S, W) = dble(q+1)
        end do

        ! main time loop
        do t = 1, Nt
            ! measure mass before streaming
            old_mass = sum(lattice(:, S:N, W:E))

            ! perform next streaming step
            call streaming_step()
            call communication_step()

            ! test mass conservation
            if (abs(old_mass - sum(lattice(:, S:N, W:E))) > EPS) then
                r = .false. ! test failed
                return
            end if

            ! chanels should be stream along c directions
            ! and periodic boundry condition should apply
            do q = 0, Nq-1
                y = S + modulo(t * c(CY, q), Ny)
                x = W + modulo(t * c(CX, q), Nx)
                if (abs(lattice(q, y, x) - dble(q+1)) > EPS) then
                    r = .false. ! test failed
                    return
                end if
            end do
        end do

        r = .true.  ! test passsed
    end function test_streaming_pbc


    logical function test_mass_conservation(Nt) result(r)
        integer, intent(in) :: Nt       !! number of time steps
        real(dp) :: old_mass            !! mass before the streaming step
        integer :: t                    !! loop time index

        do t = 1, Nt
            old_mass = sum(lattice(:, S:N, W:E))    ! mass bofore streaming
            call streaming_step()                   ! streaming step
            call communication_step()               ! communication after streaming

            ! mass before and after streaming step should not change
            if (abs(old_mass - sum(lattice(:, S:N, W:E))) > EPS) then ! note: ghost layers
                r = .false. ! test failed
                return
            end if
        end do

        r = .true. ! test passed
    end function test_mass_conservation

end module lbm_utils