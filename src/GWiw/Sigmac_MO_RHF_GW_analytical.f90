subroutine Sigmac_MO_RHF_GW_analytical(nBas,nOrb,nO,verbose,c,eHF,nfreqs,wcoord,ERI_AO,Sigma_c,Ec)


! Use the restricted Sigma_c(E) to compute the linnearized approximation to G

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: verbose
  integer,intent(in)            :: nfreqs
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nO

  double precision,intent(inout):: eHF(nOrb)
  double precision,intent(in)   :: wcoord(nfreqs)
  double precision,intent(in)   :: c(nBas,nOrb)
  double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: nV
  integer                       :: nS
  integer                       :: ifreq
  integer                       :: iorb
  integer                       :: ispin
  integer                       :: i,a,m

  double precision              :: EcRPA,EcGM
  double precision,allocatable  :: Aph(:,:)
  double precision,allocatable  :: Bph(:,:)
  double precision,allocatable  :: Om(:)
  double precision,allocatable  :: XpY(:,:)
  double precision,allocatable  :: XmY(:,:)
  double precision,allocatable  :: rho(:,:,:)
  double precision,allocatable  :: ERI_MO(:,:,:,:)

  complex *16,allocatable       :: wcoord_cpx(:)

! Ouput variables

  double precision,intent(out)  :: Ec
  complex *16,intent(out)       :: Sigma_c(nfreqs,nOrb,nOrb)

! Allocate and initialize arrays and variables
  allocate(ERI_MO(nOrb,nOrb,nOrb,nOrb))
  call AOtoMO_ERI_RHF(nBas,nOrb,c,ERI_AO,ERI_MO)

  Sigma_c=czero

! Prepare second quadrature

  allocate(wcoord_cpx(nfreqs))
  wcoord_cpx(:)=wcoord(:)*im

! Compute linear response

  nV=nOrb-nO
  nS=nO*nV

  allocate(Aph(nS,nS))
  allocate(Bph(nS,nS))
  allocate(Om(nS))
  allocate(XpY(nS,nS))
  allocate(XmY(nS,nS))
  allocate(rho(nOrb,nOrb,nS))

  ispin  = 1
  call phRLR_A(ispin,.true.,nOrb,0,nO,nV,0,nS,1d0,eHF,ERI_MO,Aph)
  call phRLR_B(ispin,.true.,nOrb,0,nO,nV,0,nS,1d0,ERI_MO,Bph)
  call phRLR(.false.,nS,Aph,Bph,EcRPA,Om,XpY,XmY)
  call RGW_excitation_density(nOrb,0,nO,0,nS,ERI_MO,XpY,rho)

! Integration along imag. freqs contributions

  do ifreq=1,nfreqs
   
   call RGW_self_energy_iomega(0d0,wcoord_cpx(ifreq),nBas,nOrb,0,nO,nV,0,nS,eHF,Om, &       ! This is Sigma_c(iw)
   rho,EcGM,Sigma_c(ifreq,:,:))

   if(verbose/=0) then
    write(*,'(a,*(f20.8))') ' Analytic ',wcoord_cpx(ifreq)
    do iorb=1,nOrb
     write(*,'(*(f20.8))') Sigma_c(ifreq,iorb,:)
    enddo
   endif
  
  enddo

  Ec = 0d0
  do m=1,nS
   do a=nO+1,nOrb
    do i=1,nO
     Ec=Ec-4d0*rho(a,i,m)*rho(a,i,m)/(eHF(a) - eHF(i) + Om(m))
    enddo
   enddo
  enddo

! Compute new total energy and Occ numbers

  deallocate(wcoord_cpx,Aph,Bph,Om,XpY,XmY,rho)
  deallocate(ERI_MO)

end subroutine

