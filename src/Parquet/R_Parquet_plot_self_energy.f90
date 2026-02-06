subroutine R_Parquet_plot_self_energy(nOrb,nC,nO,nV,nR,nS,nOOs,nVVs,nOOt,nVVt,ERI, &
                              eh_sing_rho,eh_sing_Om,eh_trip_rho,eh_trip_Om,   &
                              ee_sing_rho,ee_sing_Om,ee_trip_rho,ee_trip_Om,   &
                              hh_sing_rho,hh_sing_Om,hh_trip_rho,hh_trip_Om,   &
                              eHF,eQP)

! Dump several GW quantities for external plotting

  implicit none
  include 'parameters.h'
  
! Input variables
  integer,intent(in)            :: nOrb,nC,nO,nV,nR
  integer,intent(in)            :: nS,nOOs,nVVs,nOOt,nVVt
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: eh_sing_rho(nOrb,nOrb,nS)
  double precision,intent(in)   :: eh_sing_Om(nS)
  double precision,intent(in)   :: eh_trip_rho(nOrb,nOrb,nS)
  double precision,intent(in)   :: eh_trip_Om(nS)
  double precision,intent(in)   :: ee_sing_rho(nOrb,nOrb,nVVs)
  double precision,intent(in)   :: ee_sing_Om(nVVs)
  double precision,intent(in)   :: ee_trip_rho(nOrb,nOrb,nVVt)
  double precision,intent(in)   :: ee_trip_Om(nVVt)
  double precision,intent(in)   :: hh_sing_rho(nOrb,nOrb,nOOs)
  double precision,intent(in)   :: hh_sing_Om(nOOs)
  double precision,intent(in)   :: hh_trip_rho(nOrb,nOrb,nOOt)
  double precision,intent(in)   :: hh_trip_Om(nOOt)
  double precision,intent(in)   :: eHF(nOrb)
  double precision,intent(in)   :: eQP(nOrb)

! Local variables

  integer                       :: p,g
  integer                       :: nP
  character(len=256)            :: fmtP
  integer                       :: nGrid
  double precision              :: eta
  double precision              :: wmin,wmax,dw
  double precision              :: tmp_Re_SigC,tmp_Im_SigC,tmp_Re_dSigC
  double precision,allocatable  :: w(:)
  double precision,allocatable  :: ReSigC(:,:),ImSigC(:,:)
  double precision,allocatable  :: Z(:,:)
  double precision,allocatable  :: A(:,:)

  eta = 0.1d0
  
! Construct grid

  nGrid = 5000
  allocate(w(nGrid),ReSigC(nOrb,nGrid),ImSigC(nOrb,nGrid),Z(nOrb,nGrid),A(nOrb,nGrid))

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
    do p=nC+1,nOrb-nR

      call R_Parquet_self_energy_iieta(p,w(g),eta,nOrb,nC,nO,nV,nR,nS,nOOs,nVVs,nOOt,nVVt,ERI, &
                                       eh_sing_rho,eh_sing_Om,eh_trip_rho,eh_trip_Om,       &
                                       ee_sing_rho,ee_sing_Om,ee_trip_rho,ee_trip_Om,       &
                                       hh_sing_rho,hh_sing_Om,hh_trip_rho,hh_trip_Om,       &
                                       eQP,tmp_Re_SigC,tmp_Im_SigC,tmp_Re_dSigC)
      
      ReSigC(p,g) = tmp_Re_SigC
      ImSigC(p,g) = tmp_Im_SigC
      Z(p,g)      = tmp_Re_dSigC

    end do !(p,w(g),eta,nOrb,nC,nO,nV,nR,nS,eGW,Om,rho)
  end do
 
  Z(:,:) = 1d0/(1d0 - Z(:,:))

! Compute spectral function

  do g=1,nGrid
    do p=nC+1,nOrb-nR
      A(p,g) = abs(ImSigC(p,g))/((w(g) - eHF(p) - ReSigC(p,g))**2 + ImSigC(p,g)**2)
    end do
  end do

  A(:,:) = A(:,:)/pi

! Dump quantities in files as a function of w

  open(unit=8 ,file='R_Parquet_ReSigC.dat')
  open(unit=9 ,file='R_Parquet_freq.dat')
  open(unit=10 ,file='R_Parquet_Z.dat')
  open(unit=11 ,file='R_Parquet_A.dat')
  open(unit=12 ,file='R_Parquet_ImSigC.dat')
  
  nP = nOrb - nR - nC
  write(fmtP, '(A,I0,A)') '(F12.6,1X,', nP, '(F12.6,1X))'
  do g=1,nGrid
    write(8 ,fmtP) w(g)*HaToeV,(ReSigC(p,g)*HaToeV,p=nC+1,nOrb-nR)
    write(9 ,fmtP) w(g)*HaToeV,((w(g)-eHF(p))*HaToeV,p=nC+1,nOrb-nR)
    write(10,fmtP) w(g)*HaToeV,(Z(p,g),p=nC+1,nOrb-nR)
    write(11,fmtP) w(g)*HaToeV,(A(p,g),p=nC+1,nOrb-nR)
    write(12,fmtP) w(g)*HaToeV,(ImSigC(p,g)*HaToeV,p=nC+1,nOrb-nR)
  end do

! Closing files

  close(unit=8)
  close(unit=9)
  close(unit=10)
  close(unit=11)
  close(unit=12)

end subroutine 
