module lbm_d2q9_solver
    use lbm_d2q9_consts
    use lbm_d2q9_types

    implicit none

    integer :: W, S, E, N           !! border layers indices
    integer :: gW, gS, gE, gN       !! ghost layers indices
    integer :: Ny, Nx               !! lengths of y and x domain

    real(dp), private :: omega, dpy, dpx     !! simulation parameters

    real(dp), allocatable :: lattice(:, :, :)   !! distribution function
    real(dp), allocatable :: sbuffer(:, :, :)   !! streaming buffer

    type(boundry_condition), private :: bc_list(16)
    integer, private :: bc_count = 0

contains

    real(dp) function measure_rho(y, x)
        integer, intent(in) :: x, y
        
        measure_rho = sum(lattice(:, y, x))
    end function

    real(dp) function measure_u(component, y, x, rho)
        integer, intent(in) :: component, y, x
        real(dp), intent(in) :: rho

        measure_u = sum(lattice(:, y, x) * c(component, :)) / rho
    end function

    real(dp) function pressure_acc(dpc, length, rho)
        real(dp), intent(in) :: dpc, rho
        integer, intent(in) :: length
        real(dp) :: pressure_force

        pressure_force = dpc / length

        pressure_acc = pressure_force / rho ! or ... / (2.0_dp * rho)
    end function

    function feq(uy, ux, rho)
        real(dp), intent(in) :: rho, uy, ux
        real(dp) :: feq(0:Nq-1)
        real(dp) :: cu(0:Nq-1)
        real(dp) :: u2
        
        cu = c(CY, :) * uy + c(CX, :) * ux
        u2 = uy*uy + ux*ux
        feq = weights * rho * (1.0d0 + 3.0d0*cu + 4.5d0*cu*cu - 1.5d0*u2)
    end function


    subroutine init_lbm(&
        size_y,         &
        size_x,         &
        set_omega,      &
        set_dpy,        &
        set_dpx,        &
        set_bc_list,    &
        set_bc_count    &
    )
        integer, intent(in) :: size_y, size_x
        real(dp), intent(in) :: set_omega, set_dpy, set_dpx
        type(boundry_condition), intent(in) :: set_bc_list(Nq)
        integer, intent(in) :: set_bc_count

        Ny = size_y
        Nx = size_x

        W = 1
        S = 1
        E = size_x
        N = size_y

        gW = W - 1
        gS = S - 1
        gE = E + 1
        gN = N + 1

        omega = set_omega
        dpy = set_dpy
        dpx = set_dpx

        bc_count = set_bc_count
        bc_list(1:bc_count) = set_bc_list(1:bc_count)

        if (allocated(lattice)) deallocate(lattice)
        if (allocated(sbuffer)) deallocate(sbuffer)

        allocate(lattice(0:Nq-1, gS:gN, gW:gE))
        allocate(sbuffer(0:Nq-1, gS:gN, gW:gE))

        lattice = INIT_VAL
        sbuffer = INIT_VAL
    end subroutine init_lbm

    subroutine set_initial_condition(rho, uy, ux)
        class(scalar_field), intent(in) :: rho, uy, ux
        integer :: x, y

        do x = W, E
            do y = S, N
                lattice(:, y, x) = feq(uy%eval(y, x), ux%eval(y, x), rho%eval(y, x))
            end do
        end do
    end subroutine


    subroutine collision_step()
        real(dp) :: rho, uy, ux !! macro properties
        integer :: x, y         !! loop indicies

        do x = W, E
            do y = S, N
                rho = measure_rho(y, x)
                uy = measure_u(CY, y, x, rho) + pressure_acc(dpy, Ny, rho)
                ux = measure_u(CX, y, x, rho) + pressure_acc(dpx, Nx, rho)
                lattice(:, y, x) = (1.0d0 - omega) * lattice(:, y, x) + omega * feq(uy, ux, rho)
            end do
        end do
    end subroutine

    subroutine streaming_step()
        integer :: q ! loop index

        do q = 0, Nq-1
            sbuffer(q, S+c(CY, q):N+c(CY,q), W+c(CX,q):E+c(CX,q)) = lattice(q, S:N, W:E)
        end do

        lattice = sbuffer
    end subroutine

    subroutine communication_step()
        integer :: i ! loop index

        do i = 1, bc_count
            call bc_list(i)%apply(bc_list(i)%places, bc_list(i)%speed)
        end do
    end subroutine

end module lbm_d2q9_solver


!    subroutine streaming_step_pbc()
!        integer :: q ! loop indicies
!
!        do q = 0, Nq-1
!            lattice(q, :, :) = cshift(lattice(q, :, :), -c(CX, q), 2)
!            lattice(q, :, :) = cshift(lattice(q, :, :), -c(CY, q), 1)
!        end do
!    end subroutine streaming_step_pbc