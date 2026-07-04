module lbm_d2q9_fast_logger
	use lbm_d2q9_fast_global
	use lbm_d2q9_fast_solver
	implicit none

contains

	subroutine field_log()
		integer :: log_unit
        real(dp) :: rho, uy, ux
		integer :: y, x, gy, gx
		character(256) :: log_file

		write (log_file, '(A,A,I0,A,A,I0,A)'),  &
			trim(LOG_PARENT),   "/",            &
			num_images(),       "/",            &
			trim(LOG_FILENAME), this_image(),   &
			trim(LOG_EXT)

		open(                       &
		    newunit = log_unit,     &
		    file    = log_file,     &
		    status  = "replace",    &
		    action  = "write"       &
		)

		do y = S, N
			do x = W, E
 				rho = sum(lattice(:, y, x))
				uy = sum(lattice(:, y, x) * c(CY, :)) / rho
				ux = sum(lattice(:, y, x) * c(CX, :)) / rho

				gy = (Iy - 1) * Ny + y
				gx = (Ix - 1) * Nx + x
				write(log_unit, '(ES12.5,A,I0,A,I0,A,I0,A,ES20.10,A,ES20.10,A,ES20.10)') &
                                OMEGA,      char(9),        &
								TIME_STEPS, char(9),        &
								gy,         char(9),        &
								gx,         char(9),        &
								rho,        char(9),        &
								uy,         char(9),        &
								ux
			end do
		end do

		close(log_unit)

	end subroutine

	subroutine stats_log()
		integer :: log_unit, img
		character(256) :: log_file
		real(dp) :: mlups

		write (log_file, '(A,A)'), trim(STAT_FILENAME), trim(STAT_EXT)

		do img = 1, num_images()
			sync all

			if (this_image() == img) then

				open(                        &
				    newunit  = log_unit,     &
				    file     = log_file,     &
					position = "append",     &
				    status   = "unknown",    &
				    action   = "write"       &
				)

				mlups = real(LATTICE_SIZE, dp) * real(LATTICE_SIZE, dp)
				mlups = mlups * real(TIME_STEPS, dp)
				mlups = mlups / 1000000_dp
				mlups = mlups / total_time

				write (log_unit, '(I0,A,I0,A,G,A,G,A,G,A,G,A,G,A,G,A,G)'),  &
					num_images(),               char(9),                    &
					this_image(),               char(9),                    &
					mlups,                      char(9),                    &
					total_time,                 char(9),                    &
					total_streaming_time,       char(9),                    &
					total_slide_time,           char(9),                    &
					total_copy_time,            char(9),                    &
					total_communication_time,   char(9),                    &
					total_collision_time

				close(log_unit)

			end if
		end do

	end subroutine

end module lbm_d2q9_fast_logger