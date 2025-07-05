subroutine Xoiw_mo_RHF(nOrb,nO,eta,eHF,weval,Chi0_mo_iw)

! Restricted Xo(i w) in MO basis

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nO

  double precision,intent(in)   :: eta
  double precision,intent(inout):: eHF(nOrb)

  complex *16,intent(in)        :: weval

! Local variables

  integer                       :: iorb,aorb

  double precision              :: chem_pot

! Ouput variables

  complex *16,intent(out)       :: Chi0_mo_iw(nOrb*nOrb,nOrb*nOrb)

!------------------------------------------------------------------------
! Build Xo(i w) in MO basis
!------------------------------------------------------------------------

  chem_pot = 0.5d0*(eHF(nO)+eHF(nO+1))
  eHF(:) = eHF(:)-chem_pot

  Chi0_mo_iw=czero

  do iorb=1,nO
   do aorb=nO+1,nOrb
     Chi0_mo_iw(1+(aorb-1)+(iorb-1)*nOrb,1+(aorb-1)+(iorb-1)*nOrb)=1d0/(weval-(eHF(aorb)-eHF(iorb))+im*eta) &
                                                                  -1d0/(weval+(eHF(aorb)-eHF(iorb))-im*eta)
   enddo
  enddo
  
  Chi0_mo_iw=2d0*Chi0_mo_iw

  ! Recover eHF
  eHF(:) = eHF(:)+chem_pot

end subroutine

