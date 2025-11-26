subroutine G_MO_RHF_w(nOrb,nO,eta,eHF,wcoord,G_MO)

! G(i w)

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nO

  double precision,intent(in)   :: eta

  complex*16,intent(in)         :: wcoord

! Local variables

  integer                       :: iorb

  double precision              :: chem_pot

! Output variables
  double precision,intent(inout):: eHF(nOrb)
  complex*16,intent(out)        :: G_MO(nOrb,nOrb)
  
!--------------------------
! Build G(i w) in MO basis
!--------------------------

! write(*,*)     
! write(*,*)'**************'
! write(*,*)'* RHF G(i w) *'
! write(*,*)'**************'
! write(*,*)

 G_MO(:,:) = czero
  
 chem_pot = 0.5d0*(eHF(nO)+eHF(nO+1))
 eHF(:) = eHF(:)-chem_pot

 do iorb=1,nOrb
  if(iorb<nO+1) then
   G_MO(iorb,iorb) = 1d0/(wcoord-eHF(iorb)-im*eta)
  else
   G_MO(iorb,iorb) = 1d0/(wcoord-eHF(iorb)+im*eta)
  endif
 enddo

 ! Restore values and deallocate dyn arrays
 eHF(:) = eHF(:)+chem_pot
  
end subroutine

