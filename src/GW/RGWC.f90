subroutine RGWC(dotest,nBas,nC,nO,nV,nR,nS,Om,rho,eHF,eGW,Z)

! Perform GW+C calculation

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dotest

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS

  double precision,intent(in)   :: Om(nS)
  double precision,intent(in)   :: rho(nBas,nBas,nS)
  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: eGW(nBas)
  double precision,intent(in)   :: Z(nBas)

! Local variables

  integer                       :: p,q,i,a,m
  integer                       :: iSat
  double precision,parameter    :: cutoff = 0d0

  logical,parameter             :: do_hole_branch = .true.
  logical,parameter             :: do_electron_branch = .false.

  double precision,allocatable  :: de(:,:,:)
  double precision,allocatable  :: ZC(:)
  double precision,allocatable  :: ZSat(:,:)
  double precision,allocatable  :: eSat(:,:)

  double precision              :: eta
  integer                       :: g
  integer                       :: nGrid
  double precision              :: wmin,wmax,dw
  double precision,external     :: GW_ReSigC,GW_ImSigC,GW_RedSigC,GW_ImdSigC
  double precision,allocatable  :: w(:)
  double precision,allocatable  :: ReSigC(:,:),ImSigC(:,:)
  double precision,allocatable  :: RedSigC(:,:),ImdSigC(:,:)
  double precision,allocatable  :: AGWC(:,:)

! Output variables

! Hello world

  write(*,*)
  write(*,*)'*******************************'
  write(*,*)'* Restricted GW+C Calculation *'
  write(*,*)'*******************************'
  write(*,*)

! Memory allocation

  allocate(ZC(nBas),eSat(nBas,nS),ZSat(nBas,nS))

! Useful quantities

  allocate(de(nBas,nBas,nS))

  do p=nC+1,nBas-nR
    do i=nC+1,nO  
      do m=1,nS
        de(p,i,m) = eGW(p) - eGW(i) + Om(m)
      end do
    end do
  end do

  do p=nC+1,nBas-nR
    do a=nO+1,nBas-nR
      do m=1,nS
        de(p,a,m) = eGW(p) - eGW(a) - Om(m)
      end do
    end do
  end do

! GW+C weights

  ZC(:) = 0d0

  do p=nC+1,nBas-nR
    do q=nC+1,nBas-nR
      do m=1,nS
        ZC(p) = ZC(p) + 2d0*rho(p,q,m)**2/de(p,q,m)**2
      end do
    end do
  end do

  ZC(:) = exp(-ZC(:))

! Dump results

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)' GW+C calculation '
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X)') &
            '|','#','|','e_GW (eV)','|','e_GW+C (eV)','|','Z_GW','|','Z_GW+C','|'
  write(*,*)'-------------------------------------------------------------------------------'

  do p=1,nBas
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
    '|',p,'|',eGW(p)*HaToeV,'|',eGW(p)*HaToeV,'|',Z(p),'|',ZC(p),'|'
  end do

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)

! GW+C satellites on hole branch

  if(do_hole_branch) then

    ZSat(:,:) = 0d0
    do i=nC+1,nO
      do m=1,nS
        eSat(i,m) = eGW(i) - Om(m)
        do q=nC+1,nBas-nR
          ZSat(i,m) = ZSat(i,m) + ZC(i)*2d0*rho(i,q,m)**2/de(i,q,m)**2
        end do
      end do
    end do
 
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)' Satellite series from GW+C on hole branch'
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(1X,A5,1X,A5,1X,A5,1X,A15,1X,A15,1X)') '#','i','m','e_Sat (eV)','Z_Sat'
 
    write(*,*)'-------------------------------------------------------------------------------'
    iSat = 0
    do i=nC+1,nO
      do m=1,nS
        iSat = iSat + 1
        if(ZSat(i,m) > cutoff) &
          write(*,'(1X,I5,1X,I5,1X,I5,F15.6,1X,F15.6,1X)') iSat,i,m,eSat(i,m)*HaToeV,ZSat(i,m)
      end do
    end do
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)

  end if

! GW+C satellites on electron branch

  if(do_electron_branch) then

    ZSat(:,:) = 0d0
    do a=nO+1,nBas-nR
      do m=1,nS
        eSat(a,m) = eGW(a) + Om(m)
        do q=nC+1,nBas-nR
          ZSat(a,m) = ZSat(a,m) + ZC(a)*2d0*rho(a,q,m)**2/de(a,q,m)**2
        end do
      end do
    end do
 
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)' Satellite series from GW+C on electron branch'
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(1X,A5,1X,A5,1X,A5,1X,A15,1X,A15,1X)') '#','a','m','e_Sat (eV)','Z_Sat'
 
    write(*,*)'-------------------------------------------------------------------------------'
    iSat = 0
    do a=nO+1,nBas-nR
      do m=1,nS
        iSat = iSat + 1
        if(ZSat(a,m) > cutoff) &
          write(*,'(1X,I5,1X,I5,1X,I5,F15.6,1X,F15.6,1X)') iSat,a,m,eSat(a,m)*HaToeV,ZSat(a,m)
      end do
    end do
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)

  end if

! Broadening parameter

  eta = 0.01d0

! Construct grid

  nGrid = 1000
  allocate(w(nGrid),AGWC(nBas,nGrid))

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
      AGWC(p,g) = ZC(p)*abs(ImSigC(p,g))/((w(g) - eGW(p))**2 + ImSigC(p,g)**2)
    end do
  end do

  AGWC(:,:) = AGWC(:,:)/pi

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
      do i=nC+1,nO
        do m=1,nS
          AGWC(i,g) = AGWC(i,g) + ZSat(i,m)*abs(ImSigC(i,g))/((w(g) - eSat(i,m))**2 + ImSigC(i,g)**2)
        end do
      end do
    end do

  end if

  if(do_electron_branch) then

    do g=1,nGrid
      do a=nO+1,nBas-nR
        do m=1,nS
          AGWC(a,g) = AGWC(a,g) + ZSat(a,m)*abs(ImSigC(a,g))/((w(g) - eSat(a,m))**2 + ImSigC(a,g)**2)
        end do
      end do
    end do

  end if

  AGWC(:,:) = AGWC(:,:)/pi

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
!   call dump_test_value('R','G0W0 HOMO energy',eGW(nO))
!   call dump_test_value('R','G0W0 LUMO energy',eGW(nO+1))

! end if

end subroutine 
