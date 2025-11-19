subroutine EcRPA_EcGM_w_RHFB(nOrb,nOrb_twice,verbose,eHFB,nfreqs,ntimes,wweight,wcoord,vMAT,&
                             U_QP,EHFB_tot,EcRPA,EcGM)

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
  double precision              :: trace1,trace2,trace3
  double precision,allocatable  :: Mat1(:,:)
  double precision,allocatable  :: Mat2(:,:)
  double precision,allocatable  :: Tmp_mo_w(:,:)
  double precision,allocatable  :: eigval_Xov(:)

  complex *16                   :: weval
  complex *16,allocatable       :: Chi0_mo_w(:,:)

! Ouput variables

  double precision,intent(out)  :: EcRPA,EcGM

!

  nOrb2=nOrb*nOrb
  EcRPA=0d0; EcGM=0d0;
  allocate(eigval_Xov(nOrb2))
  allocate(Chi0_mo_w(nOrb2,nOrb2),Tmp_mo_w(nOrb2,nOrb2))
  allocate(Mat1(nOrb,nOrb),Mat2(nOrb,nOrb))
  Mat1(1:nOrb,1:nOrb)=U_QP(1:nOrb,1:nOrb)
  Mat2(1:nOrb,1:nOrb)=U_QP(nOrb+1:nOrb_twice,1:nOrb)

! Imaginary freqs contribution

  do ifreq=1,nfreqs
   
   ! Initialize
   trace1=0d0; trace2=0d0; trace3=0d0;
   Tmp_mo_w=0d0

   ! Xo (iw)
   if(ntimes>0) then
    stop
   else
    eta=0d0
    call Xoiw_mo_RHFB(nOrb,nOrb_twice,eta,eHFB,im*wcoord(ifreq),Mat1,Mat2,Chi0_mo_w)
   endif

   ! Tr [ Xo v ] and define eps
   Tmp_mo_w(:,:)=matmul(Real(Chi0_mo_w(:,:)),vMAT(:,:))
   do iorb=1,nOrb2
    trace1=trace1+Tmp_mo_w(iorb,iorb)
   enddo

   ! Diagonalize Xo v
   call diagonalize_general_matrix(nOrb2,Tmp_mo_w,eigval_Xov,Tmp_mo_w)

   ! GM 
   ! Tr [ X v ] = Tr [ (1 - Xo v)^-1 Xo v ] = Tr [ (I - U^-1 Xo v U)^-1 (U^-1 Xo v U) ] = Tr [ ( I - Eigval_mat)^-1 Eigval_mat ]
   !            = sum_p eigval_p / (1-eigval_p)
   ! RPA 
   ! Tr [ Log(I - Xo v) ] = Tr [ U^-1 Log(I - U^-1 Xo v U) U ] = Tr [ U U^-1 Log(I- Eigval_mat) ] = Tr [ Log(I - Eigval_mat) ]
   !            = sum_p Log(1-eigval_p)
   do iorb=1,nOrb2
    trace2=trace2+Log(abs(1d0-eigval_Xov(iorb)))
    trace3=trace3+eigval_Xov(iorb)/(1d0-eigval_Xov(iorb))
   enddo

   ! Compute EcRPA and EcGM from traces 
   EcRPA=EcRPA+wweight(ifreq)*(trace2+trace1)/(2d0*pi)
   EcGM =EcGM -wweight(ifreq)*(trace3-trace1)/(2d0*pi)

  enddo

! Print results
 
  if(verbose/=0) then 
   write(*,*)
   write(*,*) '******************************************************'
   write(*,*) '* EcRPA and EcGM computed with imaginary frequencies *'
   write(*,*) '******************************************************'
   write(*,*)
   write(*,*)'-------------------------------------------------------------------------------'
   write(*,'(2X,A60,F15.6,A3)') '         phRPA correlation energy = ',EcRPA,' au'
   write(*,'(2X,A60,F15.6,A3)') '         phRPA total energy       = ',EHFB_tot+EcRPA,' au'
   write(*,'(2X,A60,F15.6,A3)') '            GM correlation energy = ',EcGM,' au'
   write(*,'(2X,A60,F15.6,A3)') '            GM total energy       = ',EHFB_tot+EcGM,' au'
   write(*,*)'-------------------------------------------------------------------------------'
   write(*,*)
  endif

  ! Deallocate arrays
  deallocate(Chi0_mo_w)
  deallocate(Tmp_mo_w)
  deallocate(eigval_Xov)
  deallocate(Mat1,Mat2)

end subroutine

