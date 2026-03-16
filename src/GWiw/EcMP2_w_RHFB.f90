subroutine EcMP2_w_RHFB(nOrb,nOrb_twice,verbose,eHFB,nfreqs,ntimes,wweight,wcoord,vMAT,&
                        U_QP,EHFB_tot,EcMP2)

! Restricted Sigma_c(E)

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: verbose
  integer,intent(in)            :: nfreqs
  integer,intent(in)            :: ntimes
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nOrb_twice

  double precision,intent(in)   :: EHFB_tot
  double precision,intent(inout):: eHFB(nOrb_twice)
  double precision,intent(in)   :: wweight(nfreqs)
  double precision,intent(in)   :: wcoord(nfreqs)
  double precision,intent(in)   :: U_QP(nOrb_twice,nOrb_twice)
  double precision,intent(in)   :: vMAT(nOrb*nOrb,nOrb*nOrb)

! Local variables

  integer                       :: ifreq
  integer                       :: iorb
  integer                       :: nOrb2

  double precision              :: eta
  double precision              :: trace1
  double precision,allocatable  :: Mat1(:,:)
  double precision,allocatable  :: Mat2(:,:)
  double precision,allocatable  :: Tmp_mo_w(:,:)
  double precision,allocatable  :: eigval_Xov(:)

  complex *16                   :: weval
  complex *16,allocatable       :: Chi0_mo_w(:,:)

! Ouput variables

  double precision,intent(out)  :: EcMP2

!

  nOrb2=nOrb*nOrb
  EcMP2=0d0
  allocate(eigval_Xov(nOrb2))
  allocate(Chi0_mo_w(nOrb2,nOrb2),Tmp_mo_w(nOrb2,nOrb2))
  allocate(Mat1(nOrb,nOrb),Mat2(nOrb,nOrb))
  Mat1(1:nOrb,1:nOrb)=U_QP(1:nOrb,1:nOrb)
  Mat2(1:nOrb,1:nOrb)=U_QP(nOrb+1:nOrb_twice,1:nOrb)

! Imaginary freqs contribution

  do ifreq=1,nfreqs
   
   ! Initialize
   trace1=0d0;
   Tmp_mo_w=0d0

   ! Xo (iw)
   if(ntimes>0) then
    stop
   else
    eta=0d0
    call Xo_MO_RHFB_w(nOrb,nOrb_twice,eta,eHFB,im*wcoord(ifreq),Mat1,Mat2,Chi0_mo_w)
   endif

   ! Tr [ ( Xo v )^2 ]
   Tmp_mo_w(:,:)=matmul(Real(Chi0_mo_w(:,:)),vMAT(:,:))
   Tmp_mo_w(:,:)=matmul(Tmp_mo_w(:,:),Tmp_mo_w(:,:))
   do iorb=1,nOrb2
    trace1=trace1+Tmp_mo_w(iorb,iorb)
   enddo

   ! Compute EcMP2
   EcMP2=EcMP2-wweight(ifreq)*trace1/(8d0*pi)

  enddo

! Print results
 
  if(verbose/=0) then 
   write(*,*)
   write(*,*) '*********************************************'
   write(*,*) '* EcMP2 computed with imaginary frequencies *'
   write(*,*) '*********************************************'
   write(*,*)
   write(*,*)'-------------------------------------------------------------------------------'
   write(*,'(2X,A60,F15.6,A3)') '           MP2 correlation energy = ',EcMP2,' au'
   write(*,'(2X,A60,F15.6,A3)') '           MP2 total energy       = ',EHFB_tot+EcMP2,' au'
   write(*,*)'-------------------------------------------------------------------------------'
   write(*,*)
  endif

  ! Deallocate arrays
  deallocate(Chi0_mo_w)
  deallocate(Tmp_mo_w)
  deallocate(eigval_Xov)
  deallocate(Mat1,Mat2)

end subroutine

