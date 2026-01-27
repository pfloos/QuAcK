subroutine UMP(dotest,doMP2,doMP3,regularize,MOM,nCVS,occupations,nBas,nC,nO,nV,nR,ERI_aaaa,ERI_aabb,ERI_bbbb,ENuc,EHF,epsHF)

! Moller-Plesset module

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dotest

  logical,intent(in)            :: doMP2
  logical,intent(in)            :: doMP3

  logical,intent(in)            :: regularize
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC(nspin)
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nV(nspin)
  integer,intent(in)            :: nR(nspin)
  logical,intent(in)            :: MOM
  integer,intent(in)            :: nCVS(nspin)
  integer,intent(in)            :: occupations(maxval(nO),nspin)
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: EHF
  double precision,intent(in)   :: epsHF(nBas,nspin)
  double precision,intent(in)   :: ERI_aaaa(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_aabb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_bbbb(nBas,nBas,nBas,nBas)

! Local variables

  double precision              :: start_MP     ,end_MP       ,t_MP
  double precision              :: Ec(nsp)

! Output variables

!------------------------------------------------------------------------
! Compute MP2 energy
!------------------------------------------------------------------------                               

  if(doMP2) then    
    if(.not. MOM) then   
       call wall_time(start_MP)
       call UMP2(dotest,nBas,nC,nO,nV,nR,ERI_aaaa,ERI_aabb,ERI_bbbb,ENuc,EHF,epsHF,Ec)
       call wall_time(end_MP)
    else
       call wall_time(start_MP)
       call MOM_UMP2(dotest,nBas,nC,nO,nV,nR,nCVS,occupations,ERI_aaaa,ERI_aabb,ERI_bbbb,ENuc,EHF,epsHF,Ec)
       call wall_time(end_MP)
    end if
    t_MP = end_MP - start_MP
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for MP2 = ',t_MP,' seconds'
    write(*,*)
                    
  end if

!------------------------------------------------------------------------
! Compute MP3 energy
!------------------------------------------------------------------------                               
                    
  if(doMP3) then    

    call wall_time(start_MP)
    write(*,*) 'MP3 NYI for UHF reference'
    call wall_time(end_MP)

    t_MP = end_MP - start_MP
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for MP2 = ',t_MP,' seconds'
    write(*,*)

  end if

end subroutine 
