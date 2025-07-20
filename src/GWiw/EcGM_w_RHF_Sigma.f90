subroutine EcGM_w_RHF_Sigma(nOrb,nO,verbose,eHF,nfreqs,wweight,wcoord,vMAT,&
                            ERHF,EcGM)

! Restricted Sigma_c(E)

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: verbose
  integer,intent(in)            :: nfreqs
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nO

  double precision,intent(in)   :: ERHF
  double precision,intent(inout):: eHF(nOrb)
  double precision,intent(in)   :: wweight(nfreqs)
  double precision,intent(in)   :: wcoord(nfreqs)
  double precision,intent(in)   :: vMAT(nOrb*nOrb,nOrb*nOrb)

! Local variables

  integer                       :: kind_int
  integer                       :: ifreq
  integer                       :: iorb,jorb
  integer                       :: nOrb2
  integer                       :: nfreqs2

  double precision              :: eta
  double precision              :: lim_inf,lim_sup
  double precision              :: alpha,beta
  double precision,allocatable  :: wweight2(:)
  double precision,allocatable  :: wcoord2(:)

  complex *16                   :: trace
  complex *16,allocatable       :: wcoord2_cpx(:)
  complex *16,allocatable       :: Sigma_c(:,:,:)
  complex *16,allocatable       :: Tmp_mo(:,:)

! Ouput variables

  double precision,intent(out)  :: EcGM

!
  eta=0d0
  nfreqs2=10*nfreqs
  allocate(Sigma_c(nfreqs2,nOrb,nOrb))
  allocate(Tmp_mo(nOrb,nOrb))
  EcGM=0d0
  trace=czero

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

   call build_Sigmac_w_RHF(nOrb,nO,nfreqs2,eta,0,wcoord2_cpx,eHF,nfreqs,0,&
                           wweight,wcoord,vMAT,Sigma_c)

! Imaginary freqs contribution

  do ifreq=1,nfreqs2
   
   call G_MO_RHF(nOrb,nO,0d0,eHF,wcoord2_cpx(ifreq),Tmp_mo) ! This is G(iw2)
   Tmp_mo(:,:)=matmul(Sigma_c(ifreq,:,:),Tmp_mo(:,:))       ! This is Sigma_c x G
 
   trace=czero 
   do iorb=1,nOrb
    trace=trace+Tmp_mo(iorb,iorb)                           !  Tr [ Sigma_c x G ]
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
   write(*,'(2X,A60,F15.6,A3)') '            GM total energy       = ',ERHF+EcGM,' au'
   write(*,*)'-------------------------------------------------------------------------------'
   write(*,*)
  endif

  ! Deallocate arrays
  deallocate(Sigma_c,Tmp_mo,wweight2,wcoord2,wcoord2_cpx)

end subroutine

