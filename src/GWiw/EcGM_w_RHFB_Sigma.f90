subroutine EcGM_w_RHFB_Sigma(nOrb,nOrb_twice,verbose,eQP_state,nfreqs,wweight,wcoord,vMAT,U_QP,&
                            ERHFB,EcGM)

! Use the restricted Sigma_c(E) to compute the R-Bogoliubov Galitskii-Migdal correlation energy

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: verbose
  integer,intent(in)            :: nfreqs
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nOrb_twice

  double precision,intent(in)   :: ERHFB
  double precision,intent(inout):: eQP_state(nOrb_twice)
  double precision,intent(in)   :: wweight(nfreqs)
  double precision,intent(in)   :: wcoord(nfreqs)
  double precision,intent(in)   :: U_QP(nOrb_twice,nOrb_twice)
  double precision,intent(in)   :: vMAT(nOrb*nOrb,nOrb*nOrb)

! Local variables

  integer                       :: kind_int
  integer                       :: ifreq
  integer                       :: iorb,jorb
  integer                       :: nfreqs2

  double precision              :: eta
  double precision              :: lim_inf,lim_sup
  double precision              :: alpha,beta

  double precision,allocatable  :: wweight2(:)
  double precision,allocatable  :: wcoord2(:)
  double precision,allocatable  :: Mat1(:,:)
  double precision,allocatable  :: Mat2(:,:)

  complex *16                   :: trace
  complex *16,allocatable       :: Sigma_c_he(:,:,:)
  complex *16,allocatable       :: Sigma_c_hh(:,:,:)
  complex *16,allocatable       :: Sigma_c_eh(:,:,:)
  complex *16,allocatable       :: Sigma_c_ee(:,:,:)
  complex *16,allocatable       :: Tmp_mo(:,:)
  complex *16,allocatable       :: wcoord2_cpx(:)

! Ouput variables

  double precision,intent(out)  :: EcGM

! Allocate and initialize arrays and variables

  eta=0d0
  nfreqs2=10*nfreqs
  allocate(Sigma_c_he(nfreqs2,nOrb,nOrb))
  allocate(Sigma_c_hh(nfreqs2,nOrb,nOrb))
  allocate(Sigma_c_eh(nfreqs2,nOrb,nOrb))
  allocate(Sigma_c_ee(nfreqs2,nOrb,nOrb))
  allocate(Tmp_mo(nOrb,nOrb))
  allocate(Mat1(nOrb,nOrb),Mat2(nOrb,nOrb))
  Mat1(1:nOrb,1:nOrb)=U_QP(1:nOrb,1:nOrb)
  Mat2(1:nOrb,1:nOrb)=U_QP(nOrb+1:nOrb_twice,1:nOrb)

! Prepare second quadrature

  kind_int = 1
  lim_inf = 0d0; lim_sup = 1d0;
  alpha = 0d0;   beta  = 0d0;
  allocate(wweight2(nfreqs2),wcoord2(nfreqs2),wcoord2_cpx(nfreqs2))
  call cgqf(nfreqs2,kind_int,alpha,beta,lim_inf,lim_sup,wcoord2,wweight2)
  wweight2(:)=wweight2(:)/((1d0-wcoord2(:))**2d0)
  wcoord2(:)=wcoord2(:)/(1d0-wcoord2(:))
  wcoord2_cpx(:)=wcoord2(:)*im

! Build Sigma_c(iw)

  call build_Sigmac_w_RHFB(nOrb,nOrb_twice,nfreqs2,eta,0,wcoord2_cpx,eQP_state,nfreqs,0,wweight,wcoord, & 
                           vMAT,U_QP,Sigma_c_he,Sigma_c_hh,Sigma_c_eh,Sigma_c_ee,.false.,.false.)

! Integration along imag. freqs contributions

  EcGM=0d0
  do ifreq=1,nfreqs2
   
   call G_MO_RHFB(nOrb,nOrb_twice,eta,eQP_state,wcoord2_cpx(ifreq),Mat1,Mat1,Mat2, Mat2,Tmp_mo) ! G_he(iw2)
   trace=czero 
   Tmp_mo(:,:)=matmul(Sigma_c_he(ifreq,:,:),Tmp_mo(:,:))  ! This is Sigma_c_he(iw2) G_he(iw2)
   do iorb=1,nOrb
    trace=trace+Tmp_mo(iorb,iorb)                  !  Compute Tr [ Sigma_c_he(iw2) G_he(iw2) ]
   enddo
   call G_MO_RHFB(nOrb,nOrb_twice,eta,eQP_state,wcoord2_cpx(ifreq),Mat2,Mat1,Mat1,-Mat2,Tmp_mo) ! G_ee(iw2)
   Tmp_mo(:,:)=matmul(Sigma_c_hh(ifreq,:,:),Tmp_mo(:,:))  ! This is Sigma_c_hh(iw2) G_ee(iw2)
   do iorb=1,nOrb
    trace=trace-Tmp_mo(iorb,iorb)                  !  Substract Tr [ Sigma_c_hh(iw2) G_ee(iw2) ]
   enddo

   EcGM=EcGM+real(trace*wweight2(ifreq))

  enddo

  EcGM=EcGM/pi

! Print results
 
  if(verbose/=0) then 
   write(*,*)
   write(*,*) '**********************************************************'
   write(*,*) '* EcGM computed using Sigma_c with imaginary frequencies *'
   write(*,*) '**********************************************************'
   write(*,*)
   write(*,*)'-------------------------------------------------------------------------------'
   write(*,'(2X,A60,F15.6,A3)') '            GM correlation energy = ',EcGM,' au'
   write(*,'(2X,A60,F15.6,A3)') '            GM total energy       = ',ERHFB+EcGM,' au'
   write(*,*)'-------------------------------------------------------------------------------'
   write(*,*)
  endif

  ! Deallocate arrays

  deallocate(Sigma_c_he)
  deallocate(Sigma_c_hh)
  deallocate(Sigma_c_eh)
  deallocate(Sigma_c_ee)
  deallocate(Mat1,Mat2)
  deallocate(Tmp_mo,wweight2,wcoord2,wcoord2_cpx)

end subroutine

