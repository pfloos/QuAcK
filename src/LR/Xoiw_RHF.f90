subroutine Xoiw_RHF(nOrb,nO,eta,eHF,weval,Chi0_mo_iw)

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

  double precision              :: start_Xoiw   ,end_Xoiw     ,t_Xoiw
  double precision              :: chem_pot

! Ouput variables

  complex *16,intent(out)       :: Chi0_mo_iw(nOrb*nOrb,nOrb*nOrb)

!

  call wall_time(start_Xoiw)

  chem_pot = 0.5d0*(eHF(nO)+eHF(nO+1))
  eHF(:) = eHF(:)-chem_pot

  Chi0_mo_iw=czero

!------------------------------------------------------------------------
! Build Xo(i w) in MO basis
!------------------------------------------------------------------------

!  write(*,*)
!  write(*,*)'************************************'
!  write(*,*)'* Build RHF Xo(i w) in MO basis    *'
!  write(*,*)'************************************'
!  write(*,*)

  do iorb=1,nO
   do aorb=nO+1,nOrb
     Chi0_mo_iw(1+(aorb-1)+(iorb-1)*nOrb,1+(aorb-1)+(iorb-1)*nOrb)=1d0/(weval-(eHF(aorb)-eHF(iorb))+im*eta) &
                                                                  -1d0/(weval+(eHF(aorb)-eHF(iorb))-im*eta)
   enddo
  enddo
  
  Chi0_mo_iw=2d0*Chi0_mo_iw

  ! Deallocate arrays
  eHF(:) = eHF(:)+chem_pot

  call wall_time(end_Xoiw)
  t_Xoiw = end_Xoiw - start_Xoiw
 ! write(*,'(A65,1X,F9.3,A8)') 'Total wall time for Xo(i w) AO = ',t_Xoiw,' seconds'
 ! write(*,*)

end subroutine

