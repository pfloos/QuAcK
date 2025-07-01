subroutine Xoiw_HFB(nOrb,nOrb_twice,eta,eHFB,weval,Mat1,Mat2,Chi0_mo_iw)

! Restricted Xo(i w) in MO basis

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nOrb_twice

  double precision,intent(in)   :: eta
  double precision,intent(inout):: eHFB(nOrb_twice)
  double precision,intent(in)   :: Mat1(nOrb,nOrb)
  double precision,intent(in)   :: Mat2(nOrb,nOrb)

  complex *16,intent(in)        :: weval

! Local variables

  integer                       :: porb,qorb,rorb,sorb
  integer                       :: Istate,Jstate

  double precision              :: start_Xoiw   ,end_Xoiw     ,t_Xoiw
  double precision              :: factor1,factor2

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

  do porb=1,nOrb
   do qorb=1,nOrb
    do rorb=qorb,nOrb
     do sorb=porb,nOrb
      do Istate=1,nOrb
       do Jstate=1,nOrb
        factor1=Mat2(porb,Istate)*Mat2(qorb,Istate)*Mat1(rorb,Jstate)*Mat1(sorb,Jstate) &
               -Mat2(porb,Istate)*Mat1(qorb,Istate)*Mat2(rorb,Jstate)*Mat1(sorb,Jstate)
        factor2=Mat2(rorb,Jstate)*Mat2(sorb,Jstate)*Mat1(porb,Istate)*Mat1(qorb,Istate) &
               -Mat1(rorb,Jstate)*Mat2(sorb,Jstate)*Mat1(porb,Istate)*Mat2(qorb,Istate)
     Chi0_mo_iw(1+(sorb-1)+(porb-1)*nOrb,1+(rorb-1)+(qorb-1)*nOrb)=    &
       +Chi0_mo_iw(1+(sorb-1)+(porb-1)*nOrb,1+(rorb-1)+(qorb-1)*nOrb)  &
       +factor1/(weval-(-eHFB(Istate)-eHFB(Jstate))+im*eta)            &
       -factor2/(weval+(-eHFB(Istate)-eHFB(Jstate))-im*eta)
       enddo
      enddo
     enddo
    enddo
   enddo
  enddo
  
  Chi0_mo_iw=4d0*Chi0_mo_iw

  ! Deallocate arrays

  call wall_time(end_Xoiw)
  t_Xoiw = end_Xoiw - start_Xoiw
 ! write(*,'(A65,1X,F9.3,A8)') 'Total wall time for Xo(i w) AO = ',t_Xoiw,' seconds'
 ! write(*,*)

end subroutine

