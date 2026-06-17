module lbm_d2q9_solver
    use lbm_d2q9_consts

    implicit none

    integer :: W, S, E, N           !! border layers indices
    integer :: gW, gS, gE, gN       !! ghost layers indices

    real(dp), allocatable :: lattice(:, :, :)   !! distribution function
    real(dp), allocatable :: sbuffer(:, :, :)   !! streaming buffer

    type(boundry_condition) :: bc_list(Nq)
    integer :: bc_count = 0

contains

    real(dp) function measure_rho(y, x)
        integer, intent(in) :: x, y
        
        measure_rho = sum(lattice(:, y, x))
    end function measure_rho

    real(dp) function measure_u(component, y, x, rho)
        integer, intent(in) :: component, y, x
        real(dp), intent(in) :: rho

!        if (rho == 0.0_dp) print *, 'WARNING! rho == ', rho

        measure_u = sum(lattice(:, y, x) * c(component, :)) / rho
    end function measure_u

    function feq(uy, ux, rho)
        real(dp), intent(in) :: rho, uy, ux
        real(dp) :: feq(0:Nq-1)
        real(dp) :: cu(0:Nq-1)
        real(dp) :: u2
        
        cu = c(CY, :) * uy + c(CX, :) * ux
        u2 = uy*uy + ux*ux
        feq = weights * rho * (1.0d0 + 3.0d0*cu + 4.5d0*cu*cu - 1.5d0*u2)
    end function feq


    subroutine init_lbm(orig_y, orig_x, size_y, size_x, init_val)
        integer, intent(in) :: orig_y, orig_x, size_y, size_x
        real(dp), intent(in) :: init_val

        W = orig_x
        S = orig_y
        E = orig_x + size_x - 1
        N = orig_y + size_y - 1

        gW = W - 1
        gS = S - 1
        gE = E + 1
        gN = N + 1

        if (allocated(lattice)) deallocate(lattice)
        if (allocated(sbuffer)) deallocate(sbuffer)

        allocate(lattice(0:Nq-1, gS:gN, gW:gE))
        allocate(sbuffer(0:Nq-1, gS:gN, gW:gE))

        lattice = init_val
        sbuffer = init_val
    end subroutine init_lbm

    subroutine set_initial_condition(rho_func, uy_func, ux_func)
        procedure(scalar_field) :: rho_func, uy_func, ux_func
        integer :: x, y

        do x = W, E
            do y = S, N
                lattice(:, y, x) = feq(uy_func(y, x), ux_func(y, x), rho_func(y, x))
            end do
        end do
    end subroutine set_initial_condition

    subroutine add_boundry_condition(func, walls, corners)
        procedure(bc_func) :: func
        integer, intent(in) :: walls, corners

        bc_count = bc_count + 1
        bc_list(bc_count)%apply => func
        bc_list(bc_count)%walls = walls
        bc_list(bc_count)%corners = corners
    end subroutine add_boundry_condition

    subroutine clear_boundry_conditions()
        bc_count = 0
    end subroutine clear_boundry_conditions


    subroutine collision_step(omega)
        real(dp), intent(in) :: omega   !
        real(dp) :: rho, uy, ux         ! macro properties
        integer :: x, y                 ! loop indicies

        do x = W, E
            do y = S, N
                rho = measure_rho(y, x)
                uy = measure_u(CY, y, x, rho)
                ux = measure_u(CX, y, x, rho)
                lattice(:, y, x) = (1.0d0 - omega) * lattice(:, y, x) + omega * feq(uy, ux, rho)
            end do
        end do
    end subroutine collision_step

    subroutine streaming_step()
        integer :: q ! loop index

        do q = 0, Nq-1
            sbuffer(q, S+c(CY, q):N+c(CY,q), W+c(CX,q):E+c(CX,q)) = lattice(q, S:N, W:E)
        end do

        lattice = sbuffer
    end subroutine streaming_step

    subroutine communication_step()
        integer :: i ! loop index

        do i = 1, bc_count
            call bc_list(i)%apply(   &
                 bc_list(i)%walls,   &
                 bc_list(i)% corners &
            )
        end do
    end subroutine communication_step

end module lbm_d2q9_solver


!    subroutine streaming_step_pbc()
!        integer :: q ! loop indicies
!
!        do q = 0, Nq-1
!            lattice(q, :, :) = cshift(lattice(q, :, :), -c(CX, q), 2)
!            lattice(q, :, :) = cshift(lattice(q, :, :), -c(CY, q), 1)
!        end do
!    end subroutine streaming_step_pbc