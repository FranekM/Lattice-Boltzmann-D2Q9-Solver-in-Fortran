module lbm_io
	implicit none
	private

	public :: prepare_output_file

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

		print *, ' !!! File already exists: ', trim(filename)
		print *, '     Override? [Y/n]'
		read(*, *) answer

		if (answer == 'Y') then
			open(newunit=unit, file=filename, status="old")
			close(unit, status="delete")
			prepare_output_file = .true.
		else
			print *, "     Skipping experiment"
			prepare_output_file = .false.
		end if

	end function prepare_output_file

end module lbm_io