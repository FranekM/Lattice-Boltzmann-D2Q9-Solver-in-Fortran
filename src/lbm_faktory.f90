module lbm_faktory
	use lbm_d2q9_boundry_conditions
	use lbm_d2q9_consts
	use lbm_d2q9_types
	use lbm_loggers
	use lbm_fields
	use lbm_config

	implicit none
	private

	public :: build_logger
	public :: build_field
	public :: build_bc_list

contains

	subroutine build_logger(log, logger_type, omega, i)
		class(logger), allocatable, intent(out) :: log
		character(*), intent(in) :: logger_type
		real(dp), intent(in) :: omega
		integer, intent(in) :: i

		select case(logger_type)
            case("point")
                allocate(log, source = point_logger(    &
                    file  = configs(i)%log%file,        &
                    ext   = configs(i)%log%ext,         &
                    omega = omega,                      &
                    y     = configs(i)%log%y,           &
                    x     = configs(i)%log%x            &
                ))

            case("profile_x")
                allocate(log, source = profile_x_logger(&
                    file = configs(i)%log%file,         &
                    ext   = configs(i)%log%ext,         &
                    omega = omega,                      &
                    y    = configs(i)%log%y             &
                ))

            case("profile_y")
                allocate(log, source = profile_y_logger(&
                    file = configs(i)%log%file,         &
                    ext   = configs(i)%log%ext,         &
                    omega = omega,                      &
                    x    = configs(i)%log%x             &
                ))

            case("field")
                allocate(log, source = field_logger(    &
                    file = configs(i)%log%file,         &
                    ext   = configs(i)%log%ext,         &
                    omega = omega                       &
                ))

        end select
	end subroutine

	subroutine build_field(field, field_cfg, i)
		class(scalar_field), allocatable, intent(out) :: field
		type(field_config), intent(in) :: field_cfg
		integer, intent(in) :: i

		select case(field_cfg%type)
            case("constant")
                allocate(field, source = constant_field(&
                    value = field_cfg%const%value       &
                ))

            case("sine")
                allocate(field, source = sine_field(    &
                    val_0   = field_cfg%sine%val_0,     &
                    eps_y   = field_cfg%sine%eps_y,     &
                    eps_x   = field_cfg%sine%eps_x,     &
                    Ly      = configs(i)%Ny,            &
                    Lx      = configs(i)%Nx,            &
                    ny      = field_cfg%sine%ny,        &
                    nx      = field_cfg%sine%nx,        &
                    mode_y  = field_cfg%sine%mode_y,    &
                    mode_x  = field_cfg%sine%mode_x     &
                ))

            case("gaussian")
                allocate(field, source = gaussian_field(&
                    val_0   = field_cfg%gauss%val_0,    &
                    sigma   = field_cfg%gauss%sigma,    &
                    eps     = field_cfg%gauss%eps,      &
                    xc      = field_cfg%gauss%xc,       &
                    yc      = field_cfg%gauss%yc        &
                ))

        end select
	end subroutine

	subroutine build_bc_list(bc_list, i)
		type(boundry_condition), intent(out) :: bc_list(:)
		integer, intent(in) :: i
		integer :: j

		do j = 1, configs(i)%bc_count
            select case (configs(i)%bc_list(j))
                case ("bounce_back")
                    bc_list(j)%apply => bounce_back_bc

                case ("periodic")
                    bc_list(j)%apply => periodic_bc

            end select
            select case(configs(i)%bc_post(j))
	            case ("moving_wall")
	                bc_list(j)%postproc => moving_wall_correction

            end select
            bc_list(j)%places = configs(i)%bc_places(j)
			bc_list(j)%speed  = configs(i)%bc_speed(j)
        end do
	end subroutine

end module lbm_faktory