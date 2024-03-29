subroutine RGWC(dotest,eta,nBas,nC,nO,nV,nR,nS,Om,rho,eHF,e,eGW,Z)

! Perform GW+C calculation

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dotest

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS

  double precision,intent(in)   :: Om(nS)
  double precision,intent(in)   :: rho(nBas,nBas,nS)
  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: e(nBas)
  double precision,intent(in)   :: eGW(nBas)
  double precision,intent(in)   :: Z(nBas)

! Local variables

  integer                       :: p,q,i,a,m
  integer                       :: iSat
  double precision              :: num,eps
  double precision,parameter    :: cutoff = 0d-3

  logical,parameter             :: do_hole_branch = .true.
  logical,parameter             :: do_electron_branch = .false.

  double precision,allocatable  :: de(:,:,:)
  double precision,allocatable  :: Re_eQP(:),Im_eQP(:)
  double precision,allocatable  :: Re_ZQP(:),Im_ZQP(:)
  double precision,allocatable  :: Re_ZSat(:,:,:),Im_ZSat(:,:,:)
  double precision,allocatable  :: Re_eSat(:,:,:),Im_eSat(:,:,:)
  double precision              :: Re_dSig,Im_dSig
  double precision              :: Re_zeta,Im_zeta

  integer                       :: g
  integer                       :: nGrid
  double precision              :: wmin,wmax,dw
  double precision,allocatable  :: w(:)
  double precision,allocatable  :: AGWC(:,:)

! Output variables

! Hello world

  write(*,*)
  write(*,*)'*******************************'
  write(*,*)'* Restricted GW+C Calculation *'
  write(*,*)'*******************************'
  write(*,*)

! Memory allocation

  allocate(Re_eQP(nBas),Im_eQP(nBas),Re_ZQP(nBas),Im_ZQP(nBas), & 
           Re_eSat(nBas,nBas,nS),Im_eSat(nBas,nBas,nS),         & 
           Re_ZSat(nBas,nBas,nS),Im_ZSat(nBas,nBas,nS))

! Useful quantities

  allocate(de(nBas,nBas,nS))

  do p=nC+1,nBas-nR
    do i=nC+1,nO  
      do m=1,nS
        de(p,i,m) = e(i) - eHF(p) - Om(m)
      end do
    end do
  end do

  do p=nC+1,nBas-nR
    do a=nO+1,nBas-nR
      do m=1,nS
        de(p,a,m) = e(a) - eHF(p) + Om(m)
      end do
    end do
  end do

! GW+C quasiparticle energies and weights

  Re_eQP(:) = eHF(:)
  Im_eQP(:) = 0d0

  Re_ZQP(:) = 0d0
  Im_ZQP(:) = 0d0

  do p=nC+1,nBas-nR

    Re_dSig = 0d0
    Im_dSig = 0d0

    do q=nC+1,nBas-nR
      do m=1,nS
        num = 2d0*rho(p,q,m)**2
        eps = de(p,q,m)

        Re_eQP(p) = Re_eQP(p) - eps*num/(eps**2 + eta**2)
        Im_eQP(p) = Im_eQP(p) - eta*num/(eps**2 + eta**2)

        Re_dSig = Re_dSig - (eps**2 - eta**2)*num/(eps**2 + eta**2)**2
        Im_dSig = Im_dSig - 2d0*eta*eps*num/(eps**2 + eta**2)**2

      end do
    end do

    Re_ZQP(p) = exp(Re_dSig)*cos(Im_dSig)
    Im_ZQP(p) = exp(Re_dSig)*sin(Im_dSig)

  end do

! Dump results

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)' GW+C calculation '
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X)') &
            '|','#','|','e_GW (eV)','|','e_GW+C (eV)','|','Z_GW','|','Z_GW+C','|'
  write(*,*)'-------------------------------------------------------------------------------'

  do p=1,nBas
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
    '|',p,'|',eGW(p)*HaToeV,'|',Re_eQP(p)*HaToeV,'|',Z(p),'|',Re_ZQP(p),'|'
  end do

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)

! Initializatio

! GW+C satellites on hole branch

  if(do_hole_branch) then

!   do p=nC+1,nBas-nR
    do p=nC+1,nO
      do i=nC+1,nO
        do m=1,nS

          eps = de(p,i,m)
          num = 2d0*rho(p,i,m)**2

          Re_eSat(p,i,m) = Re_eQP(p) + eps
          Im_eSat(p,i,m) = Im_eQP(p) - eta

          Re_zeta = (eps**2 - eta**2)*num/(eps**2 + eta**2)**2
          Im_zeta = 2d0*eta*eps*num/(eps**2 + eta**2)**2

          Re_ZSat(p,i,m) = Re_ZQP(p)*Re_zeta - Im_ZQP(p)*Im_zeta
          Im_ZSat(p,i,m) = Re_ZQP(p)*Im_zeta + Im_ZQP(p)*Re_zeta

        end do
      end do
    end do
 
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)' Satellite series from GW+C on hole branch'
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(1X,A5,1X,A5,1X,A5,1X,A15,1X,A15,1X)') '#','i','m','e_Sat (eV)','Z_Sat'
 
    write(*,*)'-------------------------------------------------------------------------------'
!   do p=nC+1,nBas-nR
    do p=nC+1,nO
      do i=nC+1,nO
        do m=1,nS
          if(Re_ZSat(p,i,m) > cutoff) &
            write(*,'(1X,I5,1X,I5,1X,I5,F15.6,1X,F15.6,1X)') p,i,m,Re_eSat(p,i,m)*HaToeV,Re_ZSat(p,i,m)
        end do
      end do
    end do
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)

  end if

! GW+C satellites on electron branch

  if(do_electron_branch) then

    do p=nC+1,nBas-nR
      do a=nO+1,nBas-nR
        do m=1,nS

          eps = de(p,a,m)
          num = 2d0*rho(p,a,m)**2

          Re_eSat(p,a,m) = Re_eQP(p) + eps
          Im_eSat(p,a,m) = Im_eQP(p) - eta

          Re_zeta = (eps**2 - eta**2)*num/(eps**2 + eta**2)**2
          Im_zeta = 2d0*eta*eps*num/(eps**2 + eta**2)**2

          Re_ZSat(p,a,m) = Re_ZQP(p)*Re_zeta - Im_ZQP(p)*Im_zeta
          Im_ZSat(p,a,m) = Re_ZQP(p)*Im_zeta + Im_ZQP(p)*Re_zeta

        end do
      end do
    end do
 
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)' Satellite series from GW+C on electron branch'
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(1X,A5,1X,A5,1X,A5,1X,A15,1X,A15,1X)') '#','a','m','e_Sat (eV)','Z_Sat'
 
    write(*,*)'-------------------------------------------------------------------------------'
    do p=nC+1,nBas-nR
      do a=nO+1,nBas-nR
        do m=1,nS
          if(Re_ZSat(p,a,m) > cutoff) &
            write(*,'(1X,I5,I5,1X,1X,I5,F15.6,1X,F15.6,1X)') p,a,m,Re_eSat(p,a,m)*HaToeV,Re_ZSat(p,a,m)
        end do
      end do
    end do
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)

  end if

! Construct grid

  nGrid = 5000
  allocate(w(nGrid),AGWC(nBas,nGrid))

! Minimum and maximum frequency values

  wmin = -10d0
  wmax = 0d0
  dw = (wmax - wmin)/dble(ngrid)

  do g=1,nGrid
    w(g) = wmin + dble(g)*dw
  end do

! Compute QP part of the spectral function

  do g=1,nGrid
    do p=nC+1,nBas-nR

      AGWC(p,g) = (Re_ZQP(p)*Im_eQP(p) + Im_ZQP(p)*(w(g) - Re_eQP(p)))/((w(g) - Re_eQP(p))**2 + Im_eQP(p)**2)

    end do
  end do

  AGWC(:,:) = - AGWC(:,:)/pi

! Dump quantities in files as a function of w

  open(unit=11 ,file='GWC_AQP.dat')

  do g=1,nGrid
    write(11,*) w(g)*HaToeV,(AGWC(p,g),p=nC+1,nBas-nR)
  end do

! Closing files

  close(unit=11)

! Compute cumulant part of the spectral function 

  AGWC(:,:) = 0d0

  if(do_hole_branch) then

    do g=1,nGrid
      do p=nC+1,nO
!     do p=nC+1,nBas-nR
        do i=nC+1,nO
          do m=1,nS
            num = Re_ZSat(p,i,m)*Im_eSat(p,i,m) + Im_ZSat(p,i,m)*(w(g) - Re_eSat(p,i,m))
            eps = (w(g) - Re_eSat(p,i,m))**2 + Im_eSat(p,i,m)**2
            AGWC(p,g) = AGWC(p,g) + num/eps
          end do
        end do
      end do
    end do

  end if

  if(do_electron_branch) then

    do g=1,nGrid
      do p=nC+1,nBas-nR
        do a=nO+1,nBas-nR
          do m=1,nS
            num = Re_ZSat(p,a,m)*Im_eSat(p,a,m) + Im_ZSat(p,a,m)*(w(g) - Re_eSat(p,a,m))
            eps = (w(g) - Re_eSat(p,a,m))**2 + Im_eSat(p,a,m)**2
            AGWC(p,g) = AGWC(p,g) + num/eps
          end do
        end do
      end do
    end do

  end if

  AGWC(:,:) = - AGWC(:,:)/pi

! Dump quantities in files as a function of w

  open(unit=12 ,file='GWC_AC.dat')

  do g=1,nGrid
    write(12,*) w(g)*HaToeV,(AGWC(p,g),p=nC+1,nBas-nR)
  end do

! Closing files

  close(unit=12)

! Testing zone

! if(dotest) then

!   call dump_test_value('R','G0W0 correlation energy',EcRPA)

! end if

end subroutine 
