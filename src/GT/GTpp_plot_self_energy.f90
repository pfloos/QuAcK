subroutine GTpp_plot_self_energy(nBas,nC,nO,nV,nR,nOOs,nVVs,nOOt,nVVt,eHF,eGT,Om1s,rho1s,Om2s,rho2s, &
                                 Om1t,rho1t,Om2t,rho2t)

! Dump several GTpp quantities for external plotting

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nOOs,nOOt
  integer,intent(in)            :: nVVs,nVVt

  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: eGT(nBas)
  double precision,intent(in)   :: Om1s(nVVs),Om1t(nVVt)
  double precision,intent(in)   :: rho1s(nBas,nBas,nVVs),rho1t(nBas,nBas,nVVt)
  double precision,intent(in)   :: Om2s(nOOs),Om2t(nOOt)
  double precision,intent(in)   :: rho2s(nBas,nBas,nOOs),rho2t(nBas,nBas,nOOt)

! Local variables

  double precision              :: eta
  integer                       :: p,g
  integer                       :: nGrid
  double precision              :: wmin,wmax,dw
  double precision,external     :: GTpp_SigC,GTpp_dSigC
  double precision,allocatable  :: w(:)
  double precision,allocatable  :: SigC(:,:)
  double precision,allocatable  :: Z(:,:)
  double precision,allocatable  :: S(:,:)

! Broadening parameter

  eta = 0.01d0

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

      SigC(p,g) = GTpp_SigC(p,w(g),eta,nBas,nC,nO,nV,nR,nOOs,nVVs,nOOt,nVVt,eGT,Om1s,rho1s,Om2s,rho2s,Om1t,rho1t,Om2t,rho2t)
      Z(p,g)    = GTpp_dSigC(p,w(g),eta,nBas,nC,nO,nV,nR,nOOs,nVVs,nOOt,nVVt,eGT,Om1s,rho1s,Om2s,rho2s,Om1t,rho1t,Om2t,rho2t)

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

  open(unit=8 ,file='GTpp_SigC.dat')
  open(unit=9 ,file='GTpp_freq.dat')
  open(unit=10 ,file='GTpp_Z.dat')
  open(unit=11 ,file='GTpp_A.dat')

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
