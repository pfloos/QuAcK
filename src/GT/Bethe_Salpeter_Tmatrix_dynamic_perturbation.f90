subroutine Bethe_Salpeter_Tmatrix_dynamic_perturbation(ispin,dTDA,eta,nBas,nC,nO,nV,nR,nS,nOOs,nVVs,nOOt,nVVt,         &
                                                       Omega1s,Omega2s,Omega1t,Omega2t,rho1s,rho2s,rho1t,rho2t,eT,eGT, &
                                                       dipole_int,OmBSE,XpY,XmY,TAs,TBs,TAt,TBt)

! Compute dynamical effects via perturbation theory for BSE@GT

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: ispin
  logical,intent(in)            :: dTDA 
  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS

  integer,intent(in)            :: nOOs
  integer,intent(in)            :: nVVs
  integer,intent(in)            :: nOOt
  integer,intent(in)            :: nVVt

  double precision,intent(in)   :: eT(nBas)
  double precision,intent(in)   :: eGT(nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)
  double precision,intent(in)   :: OmBSE(nS)
  double precision,intent(in)   :: XpY(nS,nS)
  double precision,intent(in)   :: XmY(nS,nS)

  double precision,intent(in)   :: Omega1s(nVVs)
  double precision,intent(in)   :: Omega2s(nOOs)
  double precision,intent(in)   :: rho1s(nBas,nBas,nVVs)
  double precision,intent(in)   :: rho2s(nBas,nBas,nOOs)
  double precision,intent(in)   :: Omega1t(nVVt)
  double precision,intent(in)   :: Omega2t(nOOt)
  double precision,intent(in)   :: rho1t(nBas,nBas,nVVt)
  double precision,intent(in)   :: rho2t(nBas,nBas,nOOt)

  double precision,intent(in)   :: TAs(nS,nS)
  double precision,intent(in)   :: TBs(nS,nS)
  double precision,intent(in)   :: TAt(nS,nS)
  double precision,intent(in)   :: TBt(nS,nS)

! Local variables

  integer                       :: ia

  integer                       :: maxS = 10
  double precision              :: gapGT

  double precision,allocatable  :: OmDyn(:)
  double precision,allocatable  :: ZDyn(:)
  double precision,allocatable  :: X(:)
  double precision,allocatable  :: Y(:)

  double precision,allocatable  :: dTAs(:,:)
  double precision,allocatable  :: ZAs(:,:)

  double precision,allocatable  :: dTAt(:,:)
  double precision,allocatable  :: ZAt(:,:)

! Memory allocation

  maxS = min(nS,maxS)
  allocate(OmDyn(maxS),ZDyn(maxS),X(nS),Y(nS),dTAs(nS,nS),ZAs(nS,nS),dTAt(nS,nS),ZAt(nS,nS))

  if(dTDA) then 
    write(*,*)
    write(*,*) '*** dynamical TDA activated ***'
    write(*,*)
  else
    print*, ' Beyond-TDA dynamical correction for BSE@GT NYI'
    return
  end if

  OmDyn(:) = 0d0
  ZDyn(:)  = 0d0

  do ia=1,maxS

    ! Compute dynamical T-matrix for alpha-beta block  

    call dynamic_Tmatrix_A(eta,nBas,nC,nO,nV,nR,nS,nOOs,nVVs,1d0,eGT,Omega1s,Omega2s,rho1s,rho2s,OmBSE(ia),dTAs,ZAs)
 
    ! Compute dynamical T-matrix for alpha-beta block  

    call dynamic_Tmatrix_A(eta,nBas,nC,nO,nV,nR,nS,nOOt,nVVt,1d0,eGT,Omega1t,Omega2t,rho1t,rho2t,OmBSE(ia),dTAt,ZAt)
 
    X(:) = 0.5d0*(XpY(ia,:) + XmY(ia,:))
    Y(:) = 0.5d0*(XpY(ia,:) - XmY(ia,:))

    ! First-order correction 
    
    if(ispin == 1) then 
      ZDyn(ia)  = - dot_product(X,matmul(ZAt+ZAs,X)) 
      OmDyn(ia) = - dot_product(X,matmul(dTAt+dTAs,X)) + dot_product(X,matmul(TAt+TAs,X)) 
    end if

    if(ispin == 2) then 
      ZDyn(ia)  = - dot_product(X,matmul(ZAt-ZAs,X)) 
      OmDyn(ia) = - dot_product(X,matmul(dTAt-dTAs,X)) + dot_product(X,matmul(TAt-TAs,X)) 
    end if
    
    ZDyn(ia)  = 1d0/(1d0 - ZDyn(ia))
    OmDyn(ia) = ZDyn(ia)*OmDyn(ia)

  end do

!--------------!
! Dump results !
!--------------!

  gapGT = eGT(nO+1) - eGT(nO) 

  write(*,*) '---------------------------------------------------------------------------------------------------'
  write(*,*) ' First-order dynamical correction to static Bethe-Salpeter excitation energies                     '
  write(*,*) '---------------------------------------------------------------------------------------------------'
  write(*,'(A57,F10.6,A3)') ' BSE neutral excitation must be lower than the GT gap = ',gapGT*HaToeV,' eV'
  write(*,*) '---------------------------------------------------------------------------------------------------'
  write(*,'(2X,A5,1X,A20,1X,A20,1X,A20,1X,A20)') '#','Static (eV)','Dynamic (eV)','Correction (eV)','Renorm. (eV)'
  write(*,*) '---------------------------------------------------------------------------------------------------'

  do ia=1,maxS
    write(*,'(2X,I5,5X,F15.6,5X,F15.6,5X,F15.6,5X,F15.6)') & 
      ia,OmBSE(ia)*HaToeV,(OmBSE(ia)+OmDyn(ia))*HaToeV,OmDyn(ia)*HaToeV,ZDyn(ia)
  end do

  write(*,*) '---------------------------------------------------------------------------------------------------'
  write(*,*) 

end subroutine Bethe_Salpeter_Tmatrix_dynamic_perturbation
