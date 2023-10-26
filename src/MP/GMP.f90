subroutine GMP(doMP2,regularize,nBas,nC,nO,nV,nR,ERI,ENuc,EHF,epsHF)

! Moller-Plesset module

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: doMP2

  logical,intent(in)            :: regularize
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: EHF
  double precision,intent(in)   :: epsHF(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  double precision              :: start_MP     ,end_MP       ,t_MP

! Output variables

!------------------------------------------------------------------------
! Compute MP2 energy
!------------------------------------------------------------------------                               

  if(doMP2) then    
       
    call wall_time(start_MP)
    call GMP2(regularize,nBas,nC,nO,nV,nR,ERI,ENuc,EHF,epsHF)
    call wall_time(end_MP)

    t_MP = end_MP - start_MP
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for MP2 = ',t_MP,' seconds'
    write(*,*)
                    
  end if

end subroutine 
