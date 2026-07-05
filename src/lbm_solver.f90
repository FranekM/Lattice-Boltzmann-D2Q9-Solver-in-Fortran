module lbm_d2q9_solver
    use lbm_d2q9_consts
    use lbm_d2q9_types

    implicit none

    integer :: W, S, E, N                             !! border layers indices
    integer :: gW, gS, gE, gN                         !! ghost layers indices
    integer :: Ny, Nx                                 !! lengths of y and x domain
    integer :: Iy, Ix                                 !! image y and x indices
    integer :: Niy, Nix                               !! number of image y and x indices

    real(dp), private :: omega, dpy, dpx              !! simulation parameters

    real(dp), allocatable :: lattice(:, :, :)[:, :]   !! distribution function
    real(dp), allocatable :: sbuffer(:, :, :)         !! streaming buffer

    type(boundry_condition), private :: bc_list(MAX_BC_COUNT)
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

        pressure_acc = pressure_force / (2.0_dp * rho) ! coś mnie trafi !!!
    end function

    ! TODO: avoid creation of temporary arrays
    function feq(uy, ux, rho)
        real(dp), intent(in) :: rho, uy, ux
        real(dp) :: feq(0:Nq-1)
        real(dp) :: cu(0:Nq-1)
        real(dp) :: u2
        
        cu = c(CY, :) * uy + c(CX, :) * ux
        u2 = uy*uy + ux*ux
        feq = weights * rho * (1.0d0 + 3.0d0*cu + 4.5d0*cu*cu - 1.5d0*u2)
    end function


    ! TODO: simplify access helpers and move them to another module
    integer function local_to_global_x(x, imx)
        integer, intent(in) :: x
        integer, optional, intent(in) :: imx
        integer :: im_x

        if (.not. present(imx)) then
            im_x = Ix
        else
            im_x = imx
        end if
        local_to_global_x = (im_x - 1) * Nx + x
    end function

    integer function local_to_global_y(y, imy)
        integer, intent(in) :: y
        integer, optional, intent(in) :: imy
        integer :: im_y

        if (.not. present(imy)) then
            im_y = Iy
        else
            im_y = imy
        end if
        local_to_global_y = (im_y - 1) * Ny + y
    end function

    integer function global_to_local_x(x, imx)
        integer, intent(in) :: x
        integer, optional, intent(in) :: imx
        integer :: im_x

        if (.not. present(imx)) then
            im_x = Ix
        else
            im_x = imx
        end if
        global_to_local_x = x - (im_x - 1) * Nx
    end function

    integer function global_to_local_y(y, imy)
        integer, intent(in) :: y
        integer, optional, intent(in) :: imy
        integer :: im_y

        if (.not. present(imy)) then
            im_y = Iy
        else
            im_y = imy
        end if
        global_to_local_y = y - (im_y - 1) * Ny
    end function

    logical function is_global_x_in(x, imx)
        integer, intent(in) :: x
        integer, optional, intent(in) :: imx
        integer :: im_x
        integer :: x_begin, x_end

        if (.not. present(imx)) then
            im_x = Ix
        else
            im_x = imx
        end if
        x_begin = (im_x - 1) * Nx + 1
		x_end = im_x * Nx
		is_global_x_in = x >= x_begin .and. x <= x_end
    end function

    logical function is_global_y_in(y, imy)
        integer, intent(in) :: y
        integer, optional, intent(in) :: imy
        integer :: im_y
        integer :: y_begin, y_end

        if (.not. present(imy)) then
            im_y = Iy
        else
            im_y = imy
        end if
        y_begin = (im_y - 1) * Ny + 1
		y_end = im_y * Ny
		is_global_y_in = y >= y_begin .and. y <= y_end
    end function


    ! TODO: move boundy setup helpers to another module
    function is_global_boundry(place) result(is)
		integer, intent(in) :: place
        logical :: is

        select case(place)
            case (S_WALL)
                is = (Iy == 1)
            case (N_WALL)
                is = (Iy == Niy)
            case (W_WALL)
                is = (Ix == 1)
            case (E_WALL)
                is = (Ix == Nix)
            case (NW_CORNER)
                is = (Iy == Niy .or. Ix == 1)
            case (SW_CORNER)
                is = (Iy == 1 .or. Ix == 1)
            case (SE_CORNER)
                is = (Iy == 1 .or. Ix == Nix)
            case (NE_CORNER)
                is = (Iy == Niy .or. Ix == Nix)
            case default
                is = .false.
        end select
	end function

    function is_local_corner(place) result(is)
        integer, intent(in) :: place
        logical :: is

        select case(place)
            case (NW_CORNER)
                is = .not. (Iy == Niy .and. Ix == 1)
            case (SW_CORNER)
                is = .not. (Iy == 1 .and. Ix == 1)
            case (SE_CORNER)
                is = .not. (Iy == 1 .and. Ix == Nix)
            case (NE_CORNER)
                is = .not. (Iy == Niy .and. Ix == Nix)
            case default
                is = .false.
        end select
    end function

    subroutine copy_bc(bc_in, bc_out, places)
        type(boundry_condition), intent(in) :: bc_in
        type(boundry_condition), intent(out) :: bc_out
        integer, optional, intent(in) :: places

        bc_out%postproc => bc_in%postproc
        bc_out%apply    => bc_in%apply
        bc_out%speed    =  bc_in%speed
        bc_out%places   =  bc_in%places

        if (present(places)) bc_out%places = places
    end subroutine

    subroutine copy_bc_from(bc_list_in, bc_count_in, place_from, bc_out, places)
        type(boundry_condition), intent(in) :: bc_list_in(MAX_BC_COUNT)
        type(boundry_condition), intent(out) :: bc_out
        integer, intent(in) :: bc_count_in, place_from, places
        integer :: bc

        do bc = 1, bc_count_in
            if (mask(bc_list_in(bc)%places, place_from)) then
                call copy_bc(bc_list_in(bc), bc_out, places)
                exit
            end if
        end do
    end subroutine


    subroutine init_lbm(size_y, size_x, set_omega, set_dpy, set_dpx)
        integer, intent(in) :: size_y, size_x
        real(dp), intent(in) :: set_omega, set_dpy, set_dpx

        Niy = int(sqrt(real(num_images())))
        Nix = num_images() / Niy

        Ny = size_y / Niy
        Nx = size_x / Niy

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

        Iy  = this_image(lattice, 1)
        Ix  = this_image(lattice, 2)

        omega = set_omega
        dpy = set_dpy
        dpx = set_dpx

    end subroutine

    subroutine set_initial_condition(rho, uy, ux)
        class(scalar_field), intent(in) :: rho, uy, ux
        integer :: x, y     !! loop indices
        integer :: gx, gy   !! global x, y (across images)

        do x = W, E
            gx = (Ix - 1) * Nx + x
            do y = S, N
                gy = (Iy - 1) * Ny + y
                lattice(:, y, x) = feq(uy%eval(gy, gx), ux%eval(gy, gx), rho%eval(gy, gx))
            end do
        end do
    end subroutine

    subroutine set_boundry_condition(g_bc_list, g_bc_count, internal_bc)
        type(boundry_condition), intent(in) :: g_bc_list(MAX_BC_COUNT)
        type(boundry_condition), intent(in) :: internal_bc
        integer, intent(in) :: g_bc_count
        integer :: place, shift, bc

		place = 1
        do shift = 0, 7
            place = ishft(1, shift)
            bc_count = bc_count + 1

            if (is_global_boundry(place)) then
                if (is_local_corner(place)) then
                    select case(place)
                        case (NW_CORNER)
                            if (Ix == 1) &
                                call copy_bc_from(g_bc_list, g_bc_count, W_WALL, bc_list(bc_count), place)
                            if (Iy == Niy) &
                                call copy_bc_from(g_bc_list, g_bc_count, N_WALL, bc_list(bc_count), place)
                        case (SW_CORNER)
                            if (Ix == 1) &
                                call copy_bc_from(g_bc_list, g_bc_count, W_WALL, bc_list(bc_count), place)
                            if (Iy == 1) &
                                call copy_bc_from(g_bc_list, g_bc_count, S_WALL, bc_list(bc_count), place)
                        case (SE_CORNER)
                            if (Ix == Nix) &
                                call copy_bc_from(g_bc_list, g_bc_count, E_WALL, bc_list(bc_count), place)
                            if (Iy == 1) &
                                call copy_bc_from(g_bc_list, g_bc_count, S_WALL, bc_list(bc_count), place)
                        case (NE_CORNER)
                            if (Ix == Nix) &
                                call copy_bc_from(g_bc_list, g_bc_count, E_WALL, bc_list(bc_count), place)
                            if (Iy == Niy) &
                                call copy_bc_from(g_bc_list, g_bc_count, N_WALL, bc_list(bc_count), place)
                    end select
                else
                    do bc = 1, g_bc_count
                        if (mask(g_bc_list(bc)%places, place)) then
                            call copy_bc(g_bc_list(bc), bc_list(bc_count), place)
                            exit
                        end if
                    end do
                end if

            else
                call copy_bc(internal_bc, bc_list(bc_count), place)
            end if
        end do
    end subroutine

    subroutine end_lbm()
        integer :: bc

        W = -1
        S = -1
        E = -1
        N = -1
        gW = -1
        gS = -1
        gE = -1
        gN = -1
        Ny = -1
        Nx = -1
        Iy = -1
        Ix = -1
        Niy = -1
        Nix = -1
        omega = -1.0
        dpy = -1.0
        dpx = -1.0

        bc_count = 0
        do bc = 1, MAX_BC_COUNT
            nullify(bc_list(bc)%apply)
            nullify(bc_list(bc)%postproc)
            bc_list(bc)%speed = 0.0_dp
            bc_list(bc)%places = 0
        end do

        sync all

        if (allocated(lattice)) deallocate(lattice)
        if (allocated(sbuffer)) deallocate(sbuffer)
    end subroutine


    subroutine collision_step()
        real(dp) :: rho, uy, ux !! macro properties
        integer :: x, y         !! loop indicies

        do x = W, E
            do y = S, N
                rho = measure_rho(y, x)
                uy = measure_u(CY, y, x, rho) + pressure_acc(dpy, Ny, rho)
                ux = measure_u(CX, y, x, rho) + pressure_acc(dpx, Nx, rho)
                lattice(:, y, x) = (1.0_dp - omega) * lattice(:, y, x) + omega * feq(uy, ux, rho)
            end do
        end do
    end subroutine

    subroutine streaming_step()
        integer :: y, x, q ! loop indices

        do y = S, N
            do x = W, E
                do q = 0, Nq-1
                    sbuffer(q, y+c(CY, q), x+c(CX,q)) = lattice(q, y, x)
                end do
            end do
        end do

        lattice = sbuffer
    end subroutine

    subroutine communication_step()
        integer :: i ! loop index

        do i = 1, bc_count
            if (associated(bc_list(i)%apply)) &
                call bc_list(i)%apply(bc_list(i)%places)
            if (associated(bc_list(i)%postproc)) &
                call bc_list(i)%postproc(bc_list(i)%places, bc_list(i)%speed)
        end do
    end subroutine

end module lbm_d2q9_solver