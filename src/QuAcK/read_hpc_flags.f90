subroutine read_hpc_flags(working_dir, switch_hpc, use_gpu)

  implicit none

  character(len=256), intent(in) :: working_dir

  logical, intent(out)           :: switch_hpc
  logical, intent(out)           :: use_gpu

  character(len=1)              :: ans
  integer                       :: ios
  character(len=256)            :: file_path

  file_path = trim(working_dir) // '/input/hpc_flags'
  open(unit=1, file=file_path, status='old', action='read', iostat=ios)

    if(ios /= 0) then

      switch_hpc = .False.
      use_gpu = .False.

    else

      switch_hpc = .False.
      use_gpu = .False.

      read(1,*) 
      read(1,*) ans
      if(ans == 'T') switch_hpc = .true.

      read(1,*) 
      read(1,*) ans
      if(ans == 'T') use_gpu = .true.

    endif

  close(unit=1)

end subroutine 
