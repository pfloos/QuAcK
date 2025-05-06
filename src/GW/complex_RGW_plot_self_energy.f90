subroutine complex_RGW_plot_self_energy(nBas,eta,nC,nO,nV,nR,nS,Re_eHF,Im_eHF,Re_eGW,Im_eGW,Om,rho)

! Dump several GW quantities for external plotting

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: Re_eHF(nBas),Im_eHF(nBas)
  double precision,intent(in)   :: Re_eGW(nBas),Im_eGW(nBas)
  complex*16,intent(in)         :: Om(nS)
  complex*16,intent(in)         :: rho(nBas,nBas,nS)

! Local variables

  integer                       :: p,g
  integer                       :: nP
  character(len=256)            :: fmtP
  integer                       :: nGrid
  double precision              :: wmin,wmax,dw
  double precision,external     :: RGW_Re_SigC,RGW_Im_SigC,RGW_Re_dSigC
  double precision,allocatable  :: w(:)
  double precision,allocatable  :: ReSigC(:,:),ImSigC(:,:)
  double precision,allocatable  :: Z(:,:)
  double precision,allocatable  :: ReDS(:,:),ImDS(:,:)
  double precision,allocatable  :: A(:,:)

! Construct grid

  nGrid = 5000
  allocate(w(nGrid),ReSigC(nBas,nGrid),ImSigC(nBas,nGrid),Z(nBas,nGrid),&
          ReDS(nBas,nGrid),ImDS(nBas,nGrid),A(nBas,nGrid))

! Initialize 

  ReSigC(:,:) = 0d0
  ImSigC(:,:) = 0d0
  Z(:,:)      = 0d0
  ReDS(:,:)      = 0d0
  ImDS(:,:)      = 0d0

! Minimum and maximum frequency values

  wmin = 0d0
  wmax = +1d0
  dw = (wmax - wmin)/dble(ngrid)

  do g=1,nGrid
    w(g) = wmin + dble(g)*dw
  end do

! Occupied part of the self-energy and renormalization factor
 
  do g=1,nGrid
    do p=nC+1,nBas-nR

      call complex_RGW_SigC_dSigC(p,eta,nBas,nC,nO,nV,nR,nS,w(g),0d0,&
              Re_eGW,Im_eGW,Om,rho,ReSigC(p,g),ImSigC(p,g),ReDS(p,g),ImDS(p,g))
    end do
  end do
  
  Z(:,:) = (1d0-ReDS(:,:))/((1d0 - ReDS(:,:))**2 + ImDS(:,:)**2)

! Compute spectral function

  do g=1,nGrid
    do p=nC+1,nBas-nR
      A(p,g) = abs(0d0 - Im_eHF(p) - ImSigC(p,g))&
              /((w(g) - Re_eHF(p) - ReSigC(p,g))**2 + (0d0 - Im_eHF(p) - ImSigC(p,g))**2)
    end do
  end do

  A(:,:) = A(:,:)/pi

! Dump quantities in files as a function of w

  open(unit=8 ,file='RGW_SigC.dat')
  open(unit=9 ,file='RGW_freq.dat')
  open(unit=10 ,file='RGW_Z.dat')
  open(unit=11 ,file='RGW_A.dat')

  nP = nBas - nR - nC 
  write(fmtP, '(A,I0,A)') '(F12.6,1X,', nP, '(F12.6,1X))'
  do g=1,nGrid
    write(8 ,fmtP) w(g)*HaToeV,(ReSigC(p,g)*HaToeV,p=nC+1,nBas-nR)
    write(9 ,fmtP) w(g)*HaToeV,((w(g)-Re_eHF(p))*HaToeV,p=nC+1,nBas-nR)
    write(10,fmtP) w(g)*HaToeV,(Z(p,g),p=nC+1,nBas-nR)
    write(11,fmtP) w(g)*HaToeV,(A(p,g),p=nC+1,nBas-nR)
  end do

! Closing files and deallocation

  close(unit=8)
  close(unit=9)
  close(unit=10)
  close(unit=11)
  deallocate(w,ReSigC,ImSigC,Z,&
          ReDS,ImDS,A)

end subroutine 
