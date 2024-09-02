subroutine RMP(dotest,doMP2,doMP3,regularize,nOrb,nC,nO,nV,nR,ERI,ENuc,ERHF,eHF)

! Moller-Plesset module

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dotest

  logical,intent(in)            :: doMP2
  logical,intent(in)            :: doMP3

  logical,intent(in)            :: regularize
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: ERHF
  double precision,intent(in)   :: eHF(nOrb)
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)

! Local variables

  double precision              :: start_MP     ,end_MP       ,t_MP
  double precision              :: Ec

! Output variables

!------------------------------------------------------------------------
! Compute MP3 energy
!------------------------------------------------------------------------                               

  if(doMP2) then    
       
    call wall_time(start_MP)
    call RMP2(dotest,regularize,nOrb,nC,nO,nV,nR,ERI,ENuc,ERHF,eHF,Ec)
    call wall_time(end_MP)

    t_MP = end_MP - start_MP
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for MP2 = ',t_MP,' seconds'
    write(*,*)
                    
  end if

!------------------------------------------------------------------------
! Compute MP3 energy
!------------------------------------------------------------------------                               
                    
  if(doMP3) then    

    call wall_time(start_MP)
    call RMP3(nOrb,nC,nO,nV,nR,ERI,ENuc,ERHF,eHF)
    call wall_time(end_MP)

    t_MP = end_MP - start_MP
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for MP2 = ',t_MP,' seconds'
    write(*,*)

  end if

end subroutine 
