module lbm_fields
	use lbm_d2q9_consts
	use lbm_d2q9_types

	implicit none
	private


	public :: constant_field
	public :: sine_field
	public :: gaussian_field


	type, extends(scalar_field) :: constant_field
	    real(dp) :: value
	contains
	    procedure :: eval => constant_eval
	end type

	type, extends(scalar_field) :: sine_field
		real(dp) :: val_0, eps_y, eps_x
		integer :: Ly, Lx, ny, nx, mode_y, mode_x
	contains
		procedure :: eval => sine_eval
	end type sine_field

	type, extends(scalar_field) :: gaussian_field
		real(dp) :: val_0, sigma, eps
		integer :: xc, yc
	contains
		procedure :: eval => gaussian_eval
	end type gaussian_field


contains

	real(dp) function constant_eval(this, y, x)
		class(constant_field), intent(in) :: this
	    integer, intent(in) :: y,x

		constant_eval = this%value
	end function

	real(dp) function sine_eval(this, y, x)
	    class(sine_field), intent(in) :: this
	    integer, intent(in) :: y,x
	    real(dp) :: ky, kx, sin_y, sin_x

	    ky = 2.0_dp * PI / this%Ly
	    kx = 2.0_dp * PI / this%Lx
	    sin_y = this%eps_y * sin(this%ny * ky * y)
	    sin_x = this%eps_x * sin(this%nx * kx * x)

	    sine_eval = this%val_0 + this%mode_x * sin_x + this%mode_y * sin_y
	end function

	real(dp) function gaussian_eval(this, y, x)
		class(gaussian_field), intent(in) :: this
		integer, intent(in) :: y, x
        real(dp) :: r2, den

		r2 = (x - this%xc)**2 + (y - this%yc)**2
		den = 2.0_dp * this%sigma**2

		gaussian_eval = this%val_0 + this%eps * exp(-r2 / den)
	end function

end module