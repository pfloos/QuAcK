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
  integer,parameter             :: maxS = 50

  logical,parameter             :: do_hole_branch = .true.
  logical,parameter             :: do_electron_branch = .false.

  double precision,allocatable  :: de(:,:,:)
  double precision,allocatable  :: eQP(:)
  double precision,allocatable  :: ZQP(:)
  double precision,allocatable  :: ZSat(:,:,:)
  double precision,allocatable  :: eSat(:,:,:)

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

  allocate(eQP(nBas),ZQP(nBas),eSat(nBas,nBas,nS),ZSat(nBas,nBas,nS))

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

  eQP(:) = eHF(:)
  ZQP(:) = 0d0

  do p=nC+1,nBas-nR
    do q=nC+1,nBas-nR
      do m=1,nS
        num = 2d0*rho(p,q,m)**2
        eps = de(p,q,m)
        eQP(p) = eQP(p) - num*eps/(eps**2 + eta**2)
        ZQP(p) = ZQP(p) - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2
      end do
    end do
  end do

  ZQP(:) = exp(ZQP(:))

! Dump results

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)' GW+C calculation '
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X)') &
            '|','#','|','e_GW (eV)','|','e_GW+C (eV)','|','Z_GW','|','Z_GW+C','|'
  write(*,*)'-------------------------------------------------------------------------------'

  do p=1,nBas
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
    '|',p,'|',eGW(p)*HaToeV,'|',eQP(p)*HaToeV,'|',Z(p),'|',ZQP(p),'|'
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
          num = ZQP(p)*2d0*rho(p,i,m)**2
          eSat(p,i,m) = eQP(p) + eps
          ZSat(p,i,m) = num*(eps**2 - eta**2)/(eps**2 + eta**2)**2
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
        do m=1,maxS
          if(ZSat(p,i,m) > cutoff) &
            write(*,'(1X,I5,1X,I5,1X,I5,F15.6,1X,F15.6,1X)') p,i,m,eSat(p,i,m)*HaToeV,ZSat(p,i,m)
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
          num = ZQP(a)*2d0*rho(p,a,m)**2
          eSat(p,a,m) = eQP(p) + eps
          ZSat(p,a,m) = num*(eps**2 - eta**2)/(eps**2 + eta**2)**2
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
        do m=1,maxS
          if(ZSat(p,a,m) > cutoff) &
            write(*,'(1X,I5,I5,1X,1X,I5,F15.6,1X,F15.6,1X)') p,a,m,eSat(p,a,m)*HaToeV,ZSat(p,a,m)
        end do
      end do
    end do
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)

  end if

! Construct grid

  nGrid = 1000
  allocate(w(nGrid),AGWC(nBas,nGrid))

! Minimum and maximum frequency values

  wmin = -10d0
  wmax = 0d0
  dw = (wmax - wmin)/dble(ngrid)

  do g=1,nGrid
    w(g) = wmin + dble(g)*dw
  end do

! Compute QP part of the spectral function

  allocate(ReSigC(nBas,nGrid),ImSigC(nBas,nGrid))

  do g=1,nGrid
    do p=nC+1,nBas-nR
      ReSigC(p,g) = GW_ReSigC(p,eQP(p),eta,nBas,nC,nO,nV,nR,nS,eHF,Om,rho)
      ImSigC(p,g) = GW_ImSigC(p,eQP(p),eta,nBas,nC,nO,nV,nR,nS,eHF,Om,rho)
    end do
  end do

  do g=1,nGrid
    do p=nC+1,nBas-nR
      AGWC(p,g) = ZQP(p)*abs(ImSigC(p,g))/((w(g) - eQP(p))**2 + ImSigC(p,g)**2)
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
      do p=nC+1,nO
!     do p=nC+1,nBas-nR
        do i=nC+1,nO
          do m=1,maxS
            ReSigC(p,g) = GW_ReSigC(p,eSat(p,i,m),eta,nBas,nC,nO,nV,nR,nS,eHF,Om,rho)
            ImSigC(p,g) = GW_ImSigC(p,eSat(p,i,m),eta,nBas,nC,nO,nV,nR,nS,eHF,Om,rho)
            AGWC(p,g) = AGWC(p,g) + ZSat(p,i,m)*abs(ImSigC(p,g))/((w(g) - eSat(p,i,m))**2 + ImSigC(p,g)**2)
          end do
        end do
      end do
    end do

  end if

  if(do_electron_branch) then

    do g=1,nGrid
      do p=nC+1,nBas-nR
        do a=nO+1,nBas-nR
          do m=1,maxS
            ReSigC(p,g) = GW_ReSigC(p,eSat(p,a,m),eta,nBas,nC,nO,nV,nR,nS,eHF,Om,rho)
            ImSigC(p,g) = GW_ImSigC(p,eSat(p,a,m),eta,nBas,nC,nO,nV,nR,nS,eHF,Om,rho)
            AGWC(p,g) = AGWC(p,g) + ZSat(p,a,m)*abs(ImSigC(p,g))/((w(g) - eSat(p,a,m))**2 + ImSigC(p,g)**2)
          end do
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

! end if

end subroutine 
