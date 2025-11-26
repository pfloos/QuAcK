subroutine Xo_MO_RHFB_w(nOrb,nOrb_twice,eta,eHFB,weval,Mat1,Mat2,Chi0_mo_iw)

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

  integer                       :: index_ps,index_rq
  integer                       :: porb,qorb,rorb,sorb
  integer                       :: Istate,Jstate,nState

  complex *16                   :: factor1,factor2

  complex *16,allocatable       :: Mat_rsI_1(:,:,:)
  complex *16,allocatable       :: Mat_rsI_2(:,:,:)
  complex *16,allocatable       :: Mat_rsI_3(:,:,:)
  complex *16,allocatable       :: Mat_rsI_4(:,:,:)

! Ouput variables

  complex *16,intent(out)       :: Chi0_mo_iw(nOrb*nOrb,nOrb*nOrb)

!------------------------------------------------------------------------
! Build Xo(i w) in MO basis (M^5 scaling)
!------------------------------------------------------------------------

  nState=nOrb
  Chi0_mo_iw=czero
  allocate(Mat_rsI_1(nOrb,nOrb,nState))
  allocate(Mat_rsI_2(nOrb,nOrb,nState))
  allocate(Mat_rsI_3(nOrb,nOrb,nState))
  allocate(Mat_rsI_4(nOrb,nOrb,nState))

  ! M^6
!  do porb=1,nOrb
!   do qorb=1,nOrb
!    do rorb=1,nOrb
!     do sorb=1,nOrb
!      do Istate=1,nState
!       do Jstate=1,nState
!        factor1=Mat2(porb,Istate)*Mat2(qorb,Istate)*Mat1(rorb,Jstate)*Mat1(sorb,Jstate) & ! G_he G_he
!               -Mat2(porb,Istate)*Mat1(qorb,Istate)*Mat2(rorb,Jstate)*Mat1(sorb,Jstate)   ! G_hh G_ee
!        factor2=Mat2(rorb,Jstate)*Mat2(sorb,Jstate)*Mat1(porb,Istate)*Mat1(qorb,Istate) & ! G_he G_he
!               -Mat1(rorb,Jstate)*Mat2(sorb,Jstate)*Mat1(porb,Istate)*Mat2(qorb,Istate)   ! G_hh G_ee
!        index_ps=1+(sorb-1)+(porb-1)*nOrb
!        index_rq=1+(rorb-1)+(qorb-1)*nOrb
!        Chi0_mo_iw(index_ps,index_rq)=Chi0_mo_iw(index_ps,index_rq)                     &
!                                  +factor1/(weval-(-eHFB(Istate)-eHFB(Jstate))+im*eta)  &
!                                  -factor2/(weval+(-eHFB(Istate)-eHFB(Jstate))-im*eta)
!       enddo
!      enddo
!     enddo
!    enddo
!   enddo
!  enddo

  ! M^5
  Mat_rsI_1=czero
  Mat_rsI_2=czero
  Mat_rsI_3=czero
  Mat_rsI_4=czero
   ! Build the intermediates (M^4)
  do rorb=1,nOrb
   do sorb=1,nOrb
    do Istate=1,nState
     do Jstate=1,nState
      Mat_rsI_1(rorb,sorb,Istate)=Mat_rsI_1(rorb,sorb,Istate)+Mat1(rorb,Jstate)*Mat1(sorb,Jstate) &
                                                             /(weval-(-eHFB(Istate)-eHFB(Jstate))+im*eta)
      Mat_rsI_2(rorb,sorb,Istate)=Mat_rsI_2(rorb,sorb,Istate)+Mat2(rorb,Jstate)*Mat1(sorb,Jstate) &
                                                             /(weval-(-eHFB(Istate)-eHFB(Jstate))+im*eta)
      Mat_rsI_3(rorb,sorb,Istate)=Mat_rsI_3(rorb,sorb,Istate)+Mat2(rorb,Jstate)*Mat2(sorb,Jstate) &
                                                             /(weval+(-eHFB(Istate)-eHFB(Jstate))-im*eta)
      Mat_rsI_4(rorb,sorb,Istate)=Mat_rsI_4(rorb,sorb,Istate)+Mat1(rorb,Jstate)*Mat2(sorb,Jstate) &
                                                             /(weval+(-eHFB(Istate)-eHFB(Jstate))-im*eta)
     enddo
    enddo
   enddo
  enddo
   ! Complete the construction of Xo (M^5)
  do porb=1,nOrb
   do qorb=1,nOrb
    do rorb=1,nOrb
     do sorb=1,nOrb
      do Istate=1,nState
       factor1=Mat2(porb,Istate)*Mat2(qorb,Istate)*Mat_rsI_1(rorb,sorb,Istate) & ! G_he G_he
              -Mat2(porb,Istate)*Mat1(qorb,Istate)*Mat_rsI_2(rorb,sorb,Istate)   ! G_hh G_ee
       factor2=Mat1(porb,Istate)*Mat1(qorb,Istate)*Mat_rsI_3(rorb,sorb,Istate) & ! G_he G_he
              -Mat1(porb,Istate)*Mat2(qorb,Istate)*Mat_rsI_4(rorb,sorb,Istate)   ! G_hh G_ee
       index_ps=1+(sorb-1)+(porb-1)*nOrb
       index_rq=1+(rorb-1)+(qorb-1)*nOrb
       Chi0_mo_iw(index_ps,index_rq)=Chi0_mo_iw(index_ps,index_rq) + factor1 - factor2
      enddo
     enddo
    enddo
   enddo
  enddo
  
  Chi0_mo_iw=2d0*Chi0_mo_iw  ! Times 2 because of the spin contributions
                             ! [ i.e., for Ghe Ghe take (up,up,up,up) and (down,down,down,down) 
                             !                                 & 
                             !       for Ghh Gee (up,down,down,up) and (down,up,up,down) ]
  deallocate(Mat_rsI_1)
  deallocate(Mat_rsI_2)
  deallocate(Mat_rsI_3)
  deallocate(Mat_rsI_4)

end subroutine

