subroutine RGWC(dotest,nBas,nC,nO,nR,nS,Om,rho,e,eGW,Z)

! Perform GW+C calculation

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dotest

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS

  double precision,intent(in)   :: Om(nS)
  double precision,intent(in)   :: rho(nBas,nBas,nS)
  double precision,intent(in)   :: e(nBas)
  double precision,intent(in)   :: eGW(nBas)
  double precision,intent(in)   :: Z(nBas)

! Local variables

  integer                       :: p,q,i,a,m
  integer                       :: iSat
  double precision,parameter    :: cutoff = 1d-3

  double precision,allocatable  :: de(:,:,:)
  double precision,allocatable  :: ZC(:)
  double precision,allocatable  :: ZSat(:,:)
  double precision,allocatable  :: eSat(:,:)

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

! GW+C satellites on electron branch

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

! Testing zone

! if(dotest) then

!   call dump_test_value('R','G0W0 correlation energy',EcRPA)
!   call dump_test_value('R','G0W0 HOMO energy',eGW(nO))
!   call dump_test_value('R','G0W0 LUMO energy',eGW(nO+1))

! end if

end subroutine 
