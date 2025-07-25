subroutine G_AO_RHF(nBas,nOrb,nO,eta,cHF,eHF,wcoord,G_AO)

! G(i w)

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nO

  double precision,intent(in)   :: eta
  double precision,intent(in)   :: cHF(nBas,nOrb)

  complex*16,intent(in)         :: wcoord
  complex*16,allocatable        :: Gtmp(:,:)

! Local variables

  integer                       :: iorb

  double precision              :: chem_pot

! Output variables
  double precision,intent(inout):: eHF(nOrb)
  complex*16,intent(out)        :: G_AO(nBas,nBas)
  
!--------------------------
! Build G(i w) in AO basis
!--------------------------

! write(*,*)     
! write(*,*)'**************'
! write(*,*)'* RHF G(i w) *'
! write(*,*)'**************'
! write(*,*)

 allocate(Gtmp(nOrb,nOrb))
 Gtmp(:,:) = czero
  
 chem_pot = 0.5d0*(eHF(nO)+eHF(nO+1))
 eHF(:) = eHF(:)-chem_pot

 do iorb=1,nOrb
  if(iorb<nO+1) then
   Gtmp(iorb,iorb) = 1d0/(wcoord-eHF(iorb)-im*eta)
  else
   Gtmp(iorb,iorb) = 1d0/(wcoord-eHF(iorb)+im*eta)
  endif
 enddo

 G_AO=matmul(matmul(cHF,Gtmp),transpose(cHF))

 ! Restore values and deallocate dyn arrays
 deallocate(Gtmp)
 eHF(:) = eHF(:)+chem_pot
  
end subroutine

