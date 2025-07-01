subroutine Xoiw_HFB(nOrb,nOrb_twice,eta,eHFB,weval,U_QP,Chi0_mo_iw)

! Restricted Xo(i w) in MO basis

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nOrb_twice

  double precision,intent(in)   :: eta
  double precision,intent(inout):: eHFB(nOrb_twice)
  double precision,intent(in)   :: U_QP(nOrb_twice,nOrb_twice)

  complex *16,intent(in)        :: weval

! Local variables

  integer                       :: iorb,jorb,korb,lorb

  double precision              :: start_Xoiw   ,end_Xoiw     ,t_Xoiw

! Ouput variables

  complex *16,intent(out)       :: Chi0_mo_iw(nOrb*nOrb,nOrb*nOrb)

!

  call wall_time(start_Xoiw)

  Chi0_mo_iw=czero

!------------------------------------------------------------------------
! Build Xo(i w) in MO basis
!------------------------------------------------------------------------

!  write(*,*)
!  write(*,*)'************************************'
!  write(*,*)'* Build HFB Xo(i w) in MO basis    *'
!  write(*,*)'************************************'
!  write(*,*)

  do iorb=1,nOrb
   do jorb=1,nOrb
     Chi0_mo_iw(1+(jorb-1)+(iorb-1)*nOrb,1+(jorb-1)+(iorb-1)*nOrb)=1d0/(weval-(eHFB(jorb)-eHFB(iorb))+im*eta) &
                                                                  -1d0/(weval+(eHFB(jorb)-eHFB(iorb))-im*eta)
   enddo
  enddo
  
  Chi0_mo_iw=2d0*Chi0_mo_iw

  ! Deallocate arrays

  call wall_time(end_Xoiw)
  t_Xoiw = end_Xoiw - start_Xoiw
 ! write(*,'(A65,1X,F9.3,A8)') 'Total wall time for Xo(i w) AO = ',t_Xoiw,' seconds'
 ! write(*,*)

end subroutine

