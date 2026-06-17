module lbm_utils
    use iso_fortran_env, only: output_unit
    use lbm_d2q9_consts
    use lbm_d2q9_solver
    use lbm_d2q9_boundry_conditions

    implicit none
    private

    public :: debug_print
    public :: write_vti
    public :: log_write_vti_progress
    public :: test_streaming_with_pbc
    public :: test_mass_conservation

    real(dp), parameter :: EPS = 1.0e-10

    integer, parameter :: q_map(0:2, 0:2) = reshape([ &
            QNW, QN, QNE,   &! array that   0  1  2      6  2  5
            QW,  QR, QE,    &! helps to     3  4  5  ->  3  0  1
            QSW, QS, QSE    &! rearange q   6  7  8      7  4  8
    ], [3,3])

contains

    subroutine debug_print()
        double precision :: output(1:(N-S+1)*3, 1:(E-W+1)*3)
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

	
    subroutine write_vti(t_step)
        integer, intent(in) :: t_step

        integer :: dim_x, dim_y ! lattice dimensions
        integer :: x, y         ! loop indicies
        integer :: unit         ! output file stream

        double precision :: rho, ux, uy ! measured desity and volocity components
        character(len=100) :: filename  ! filename (one file for each time step)

        dim_x = (E-W)
        dim_y = (N-S)

        write(filename, '("res/lbm_",I6.6,".vti")') t_step

        open(newunit=unit, file=filename, status='replace')

        write(unit, '(A)')			'<?xml version="1.0"?>'
        write(unit, '(A)')			'<VTKFile type="ImageData" version="0.1" byte_order="LittleEndian">'
        write(unit, '(A,I0,A,I0,A)')	'<ImageData WholeExtent="0 ',dim_x,' 0 ',dim_y,' 0 0" Origin="0 0 0" Spacing="1 1 1">'
        write(unit, '(A,I0,A,I0,A)')		'<Piece Extent="0 ',dim_x,' 0 ',dim_y,' 0 0">'
        write(unit, '(A)')						'<PointData Scalars="rho" Vectors="velocity">'

        write(unit, '(A)')							'<DataArray type="Float64" Name="rho" format="ascii">'
        do y = S, N
            do x = W, E
                rho = measure_rho(y, x)
                write(unit, *)							rho
            enddo
        enddo
        write(unit, '(A)')							'</DataArray>'
        write(unit, '(A)')							'<DataArray type="Float64" Name="velocity" NumberOfComponents="3" format="ascii">'
        do y = S, N
            do x = W, E
                rho = measure_rho(y, x)
                ux  = measure_u(CX, y, x, rho)
                uy  = measure_u(CY, y, x, rho)
                write(unit, *)           				ux, uy, 0.0d0
            enddo
        enddo
        write(unit, '(A)')							'</DataArray>'
        write(unit, '(A)')						'</PointData>'
        write(unit, '(A)')					'</Piece>'
        write(unit, '(A)')				'</ImageData>'
        write(unit, '(A)')			'</VTKFile>'

        close(unit)

    end subroutine


    subroutine log_write_vti_progress(t, Nt)
        integer, intent(in) :: t, Nt

        write(output_unit, '(A,I0,A,I0,A,F5.1,A)', advance='no') &
            achar(13), t, '/', Nt, ' frames ready (', real(t)/real(Nt) * 100.0, '%)'
        flush(output_unit)
    end subroutine log_write_vti_progress


    logical function test_streaming_with_pbc(Nt) result(r)
        integer, intent(in) :: Nt   ! test max time steps
        real(dp) old_mass           ! mass before streaming step
        integer :: x, y, t, q       ! loop indicies
        integer :: dim_x, dim_y     ! lattice dimensions

        dim_x = (E-W) + 1
        dim_y = (N-S) + 1

        ! initial condition: zero whole lattice
        ! set unique indetifiers only for channels at position (S, W)
        lattice = 0.0_dp
        sbuffer = 0.0_dp
        do q = 0, Nq-1
            lattice(q, S, W) = dble(q+1)
        end do
        call clear_boundry_conditions()
        call add_boundry_condition(periodic_bc, ALL_WALLS, ALL_CORNERS)

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
                y = S + modulo(t * c(CY, q), dim_y)
                x = W + modulo(t * c(CX, q), dim_x)
                if (abs(lattice(q, y, x) - dble(q+1)) > EPS) then
                    r = .false. ! test failed
                    return
                end if
            end do
        end do

        r = .true.  ! test passsed
    end function test_streaming_with_pbc


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