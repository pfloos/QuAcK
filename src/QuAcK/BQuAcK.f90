subroutine BQuAcK(working_dir,dotest,doHFB)

! Restricted branch of QuAcK

  implicit none
  include 'parameters.h'

  character(len=256),intent(in) :: working_dir

  logical,intent(in)            :: dotest

  logical,intent(in)            :: doHFB

! Local variables

  double precision              :: start_HF     ,end_HF       ,t_HF

  write(*,*)
  write(*,*) '******************************'
  write(*,*) '* Bogoliubov Branch of QuAcK *'
  write(*,*) '******************************'
  write(*,*)

!-------------------!
! Memory allocation !
!-------------------!

!---------------------!
! Hartree-Fock module !
!---------------------!

  if(doHFB) then

    call wall_time(start_HF)
!   call HFB(dotest)
    call wall_time(end_HF)

    t_HF = end_HF - start_HF
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for HFB = ',t_HF,' seconds'
    write(*,*)

  end if

end subroutine
