subroutine GWC_spectral_function(nBas,nC,nO,nV,nR,nS,eHF,eGW,Om,rho)

! Plot the spectral function at the GW+C level

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
  double precision,intent(in)   :: eGW(nBas)
  double precision,intent(in)   :: Om(nS)
  double precision,intent(in)   :: rho(nBas,nBas,nS)

! Local variables

  double precision              :: eta
  integer                       :: p,g
  integer                       :: nGrid
  double precision              :: wmin,wmax,dw
  double precision,external     :: GW_ReSigC,GW_ImSigC,GW_RedSigC,GW_ImdSigC
  double precision,allocatable  :: w(:)
  double precision,allocatable  :: ReSigC(:,:),ImSigC(:,:)
  double precision,allocatable  :: RedSigC(:,:),ImdSigC(:,:)
  double precision,allocatable  :: A(:,:)

! Broadening parameter

  eta = 0.01d0

! Construct grid

  nGrid = 5000
  allocate(w(nGrid),A(nBas,nGrid))

! Minimum and maximum frequency values

  wmin = -5d0
  wmax = +5d0
  dw = (wmax - wmin)/dble(ngrid)

  do g=1,nGrid
    w(g) = wmin + dble(g)*dw
  end do

! Compute QP part of the spectral function

  allocate(ReSigC(nBas,nGrid),ImSigC(nBas,nGrid))

  do g=1,nGrid
    do p=nC+1,nBas-nR

      ReSigC(p,g) = GW_ReSigC(p,eGW(p),eta,nBas,nC,nO,nV,nR,nS,eHF,Om,rho)
      ImSigC(p,g) = GW_ImSigC(p,eGW(p),eta,nBas,nC,nO,nV,nR,nS,eHF,Om,rho)

    end do
  end do

  do g=1,nGrid
    do p=nC+1,nBas-nR
      A(p,g) = abs(ImSigC(p,g))/((w(g) - eHF(p) - ReSigC(p,g))**2 + ImSigC(p,g)**2)
    end do
  end do

  A(:,:) = A(:,:)/pi

  deallocate(ReSigC,ImSigC)

! Dump quantities in files as a function of w

  open(unit=11 ,file='GWC_AQP.dat')

  do g=1,nGrid
    write(11,*) w(g)*HaToeV,(A(p,g),p=nC+1,nBas-nR)
  end do

! Closing files

  close(unit=11)

! Compute cumulant part of the spectral function 

  allocate(RedSigC(nBas,nGrid),ImdSigC(nBas,nGrid))

  do g=1,nGrid
    do p=nC+1,nBas-nR

      RedSigC(p,g) = GW_RedSigC(p,eHF(p),eta,nBas,nC,nO,nV,nR,nS,eHF,Om,rho)
      ImdSigC(p,g) = GW_ImdSigC(p,eHF(p),eta,nBas,nC,nO,nV,nR,nS,eHF,Om,rho)

    end do
  end do

  do g=1,nGrid
    do p=nC+1,nBas-nR
      A(p,g) = RedSigC(p,g) + (w(g) - eHF(p))*ImdSigC(p,g)
    end do
  end do

  do g=1,nGrid
    do p=nC+1,nBas-nR

      RedSigC(p,g) = GW_RedSigC(p,eHF(p)+w(g),eta,nBas,nC,nO,nV,nR,nS,eHF,Om,rho)

    end do
  end do

  do g=1,nGrid
    do p=nC+1,nBas-nR
      A(p,g) = (RedSigC(p,g) - A(p,g))/(w(g) - eHF(p))**2
    end do
  end do

  A(:,:) = A(:,:)/pi

  deallocate(RedSigC,ImdSigC)

! Dump quantities in files as a function of w

  open(unit=12 ,file='GWC_AC.dat')

  do g=1,nGrid
    write(12,*) w(g)*HaToeV,(A(p,g),p=nC+1,nBas-nR)
  end do
 
! Closing files

  close(unit=12)

end subroutine 
