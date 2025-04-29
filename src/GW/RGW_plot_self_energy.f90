subroutine RGW_plot_self_energy(nBas,eta,nC,nO,nV,nR,nS,eHF,eGW,Om,rho)

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
  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: eGW(nBas)
  double precision,intent(in)   :: Om(nS)
  double precision,intent(in)   :: rho(nBas,nBas,nS)

! Local variables

  integer                       :: p,g
  integer                       :: nGrid
  double precision              :: wmin,wmax,dw
  double precision,external     :: RGW_Re_SigC,RGW_Im_SigC,RGW_Re_dSigC
  double precision,allocatable  :: w(:)
  double precision,allocatable  :: ReSigC(:,:),ImSigC(:,:)
  double precision,allocatable  :: Z(:,:)
  double precision,allocatable  :: A(:,:)

! Construct grid

  nGrid = 5000
  allocate(w(nGrid),ReSigC(nBas,nGrid),ImSigC(nBas,nGrid),Z(nBas,nGrid),A(nBas,nGrid))

! Initialize 

  ReSigC(:,:) = 0d0
  ImSigC(:,:) = 0d0
  Z(:,:)      = 0d0

! Minimum and maximum frequency values

  wmin = -5d0
  wmax = +5d0
  dw = (wmax - wmin)/dble(ngrid)

  do g=1,nGrid
    w(g) = wmin + dble(g)*dw
  end do

! Occupied part of the self-energy and renormalization factor
 
  do g=1,nGrid
    do p=nC+1,nBas-nR

      ReSigC(p,g) = RGW_Re_SigC(p,w(g),eta,nBas,nC,nO,nV,nR,nS,eGW,Om,rho)
      ImSigC(p,g) = RGW_Im_SigC(p,w(g),eta,nBas,nC,nO,nV,nR,nS,eGW,Om,rho)
      Z(p,g)      = RGW_Re_dSigC(p,w(g),eta,nBas,nC,nO,nV,nR,nS,eGW,Om,rho)

    end do
  end do
 
  Z(:,:) = 1d0/(1d0 + Z(:,:))

! Compute spectral function

  do g=1,nGrid
    do p=nC+1,nBas-nR
      A(p,g) = abs(ImSigC(p,g))/((w(g) - eHF(p) - ReSigC(p,g))**2 + ImSigC(p,g)**2)
    end do
  end do

  A(:,:) = A(:,:)/pi

! Dump quantities in files as a function of w

  open(unit=8 ,file='RGW_SigC.dat')
  open(unit=9 ,file='RGW_freq.dat')
  open(unit=10 ,file='RGW_Z.dat')
  open(unit=11 ,file='RGW_A.dat')
  
  do g=1,nGrid
    write(8 ,'(F12.6,1X,150(F12.6,1X))') w(g)*HaToeV,(ReSigC(p,g)*HaToeV,p=nC+1,nBas-nR)
    write(9 ,'(F12.6,1X,150(F12.6,1X))') w(g)*HaToeV,((w(g)-eHF(p))*HaToeV,p=nC+1,nBas-nR)
    write(10,'(F12.6,1X,150(F12.6,1X))') w(g)*HaToeV,(Z(p,g),p=nC+1,nBas-nR)
    write(11,'(F12.6,1X,150(F12.6,1X))') w(g)*HaToeV,(A(p,g),p=nC+1,nBas-nR)
  end do

! Closing files

  close(unit=8)
  close(unit=9)
  close(unit=10)
  close(unit=11)

end subroutine 
