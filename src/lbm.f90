program lbm_d2q9
	use lbm_d2q9_global
	use lbm_d2q9_solver
	use lbm_d2q9_logger
	implicit none

	integer :: t


	if (this_image() == 1)                                           &
		print '(A,A,I0,A)',  'LBM D2Q9 LID DRIVEN CAVITY ',          &
							 'starting on: ', num_images(), ' images'

	call init_lbm()

	call cpu_time(total_start_time)
	do t = 0, TIME_STEPS

		call cpu_time(st)
		call collision_step()
		call cpu_time(et)
		total_collision_time = total_collision_time + (et - st)

		call cpu_time(st)
		call communication_step()
		call cpu_time(et)
		total_communication_time = total_communication_time + (et - st)

		call cpu_time(st)
		call streaming_step()
		sync all
		call cpu_time(et)
		total_streaming_time = total_streaming_time + (et - st)

	end do
	call cpu_time(total_end_time)
	total_time = (total_end_time - total_start_time)

	call stats_log()
	call field_log()

	call end_lbm()

	sync all
	if (this_image() == 1)                                           &
		print '(A,A,I0,A)',  'LBM D2Q9 LID DRIVEN CAVITY ',          &
					         'finished on: ', num_images(), ' images'

end program lbm_d2q9