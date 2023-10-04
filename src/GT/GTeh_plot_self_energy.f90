subroutine GTeh_plot_self_energy(nBas,nC,nO,nV,nR,nS,eHF,eGT,Om,rhoL,rhoR)

! Dump several GTeh quantities for external plotting

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: eGT(nBas)
  double precision,intent(in)   :: Om(nS)
  double precision,intent(in)   :: rhoL(nBas,nBas,nS)
  double precision,intent(in)   :: rhoR(nBas,nBas,nS)

! Local variables

  double precision              :: eta
  integer                       :: p,g
  integer                       :: nGrid
  double precision              :: wmin,wmax,dw
  double precision,external     :: GTeh_SigC,GTeh_dSigC
  double precision,allocatable  :: w(:)
  double precision,allocatable  :: SigC(:,:)
  double precision,allocatable  :: Z(:,:)
  double precision,allocatable  :: S(:,:)

! Broadening parameter

  eta = 0.1d0

! Construct grid

  nGrid = 5000
  allocate(w(nGrid),SigC(nBas,nGrid),Z(nBas,nGrid),S(nBas,nGrid))

! Initialize 

  SigC(:,:) = 0d0
  Z(:,:)    = 0d0

! Minimum and maximum frequency values

  wmin = -5d0
  wmax = +5d0
  dw = (wmax - wmin)/dble(ngrid)

  do g=1,nGrid
    w(g) = wmin + dble(g)*dw
  enddo

! Occupied part of the self-energy and renormalization factor
 
  do g=1,nGrid
    do p=nC+1,nBas-nR

      SigC(p,g) = GTeh_SigC(p,w(g),eta,nBas,nC,nO,nV,nR,nS,eGT,Om,rhoL,rhoR)
      Z(p,g)    = GTeh_dSigC(p,w(g),eta,nBas,nC,nO,nV,nR,nS,eGT,Om,rhoL,rhoR)

    end do
  end do
 
  Z(:,:) = 1d0/(1d0 + Z(:,:))

! Compute spectral function

  do g=1,nGrid
    do p=nC+1,nBas-nR
      S(p,g) = eta/((w(g) - eHF(p) - SigC(p,g))**2 + eta**2)
    enddo
  enddo

  S(:,:) = S(:,:)/pi

! Dump quantities in files as a function of w

  open(unit=8 ,file='GTeh_SigC.dat')
  open(unit=9 ,file='GTeh_freq.dat')
  open(unit=10 ,file='GTeh_Z.dat')
  open(unit=11 ,file='GTeh_A.dat')

  do g=1,nGrid
    write(8 ,*) w(g)*HaToeV,(SigC(p,g)*HaToeV,p=nC+1,nBas-nR)
    write(9 ,*) w(g)*HaToeV,((w(g)-eHF(p))*HaToeV,p=nC+1,nBas-nR)
    write(10,*) w(g)*HaToeV,(Z(p,g),p=nC+1,nBas-nR)
    write(11,*) w(g)*HaToeV,(S(p,g),p=nC+1,nBas-nR)
  enddo

! Closing files

  close(unit=8)
  close(unit=9)
  close(unit=10)
  close(unit=11)

end subroutine 
