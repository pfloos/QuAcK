subroutine Bethe_Salpeter_Tmatrix_dynamic_perturbation(dTDA,eta,nBas,nC,nO,nV,nR,nS,nOO,nVV,Omega1,Omega2,rho1,rho2, &
                                                       eT,eGT,dipole_int,OmBSE,XpY,XmY)
! Compute dynamical effects via perturbation theory for BSE@GT

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dTDA 
  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS

  integer,intent(in)            :: nOO
  integer,intent(in)            :: nVV

  double precision,intent(in)   :: eT(nBas)
  double precision,intent(in)   :: eGT(nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)
  double precision,intent(in)   :: OmBSE(nS)
  double precision,intent(in)   :: XpY(nS,nS)
  double precision,intent(in)   :: XmY(nS,nS)

  double precision,intent(in)   :: Omega1(nVV)
  double precision,intent(in)   :: Omega2(nOO)
  double precision,intent(in)   :: rho1(nBas,nBas,nVV)
  double precision,intent(in)   :: rho2(nBas,nBas,nOO)

! Local variables

  integer                       :: ia

  integer,parameter             :: maxS = 10
  double precision              :: gapGT

  double precision,allocatable  :: OmDyn(:)
  double precision,allocatable  :: ZDyn(:)
  double precision,allocatable  :: X(:)
  double precision,allocatable  :: Y(:)

  double precision,allocatable  ::  Ap_dyn(:,:)
  double precision,allocatable  :: ZAp_dyn(:,:)

  double precision,allocatable  ::  Bp_dyn(:,:)
  double precision,allocatable  :: ZBp_dyn(:,:)

  double precision,allocatable  ::  Am_dyn(:,:)
  double precision,allocatable  :: ZAm_dyn(:,:)

  double precision,allocatable  ::  Bm_dyn(:,:)
  double precision,allocatable  :: ZBm_dyn(:,:)

! Memory allocation

  allocate(OmDyn(nS),ZDyn(nS),X(nS),Y(nS),Ap_dyn(nS,nS),ZAp_dyn(nS,nS))

  if(.not.dTDA) allocate(Am_dyn(nS,nS),ZAm_dyn(nS,nS),Bp_dyn(nS,nS),ZBp_dyn(nS,nS),Bm_dyn(nS,nS),ZBm_dyn(nS,nS))

  if(dTDA) then 
    write(*,*)
    write(*,*) '*** dynamical TDA activated ***'
    write(*,*)
  end if

  gapGT = eGT(nO+1) - eGT(nO) 

  write(*,*) '---------------------------------------------------------------------------------------------------'
  write(*,*) ' First-order dynamical correction to static Bethe-Salpeter excitation energies                     '
  write(*,*) '---------------------------------------------------------------------------------------------------'
  write(*,'(A57,F10.6,A3)') ' BSE neutral excitation must be lower than the GT gap = ',gapGT*HaToeV,' eV'
  write(*,*) '---------------------------------------------------------------------------------------------------'
  write(*,'(2X,A5,1X,A20,1X,A20,1X,A20,1X,A20)') '#','Static (eV)','Dynamic (eV)','Correction (eV)','Renorm. (eV)'
  write(*,*) '---------------------------------------------------------------------------------------------------'

  do ia=1,min(nS,maxS)

    X(:) = 0.5d0*(XpY(ia,:) + XmY(ia,:))
    Y(:) = 0.5d0*(XpY(ia,:) - XmY(ia,:))

    ! First-order correction 

    if(dTDA) then 

      ! Resonant part of the BSE correction for dynamical TDA

      call dynamic_Tmatrix_A(eta,nBas,nC,nO,nV,nR,nS,nOO,nVV,1d0,eGT,Omega1,Omega2,rho1,rho2,OmBSE(ia),Ap_dyn,Zap_dyn)

      ZDyn(ia)  = dot_product(X,matmul(ZAp_dyn,X))
      OmDyn(ia) = dot_product(X,matmul( Ap_dyn,X))

    else

      print*, ' Beyond-TDA dynamical correction for BSE@GT NYI'
      ! Resonant and anti-resonant part of the BSE correction

!     call dynamic_Tmatrix_TAB(eta,nBas,nC,nO,nV,nR,nS,nOO,nVV,1d0,eGT,Omega1,Omega2,rho1,rho2,OmBSE(ia), &
!                              Ap_dyn,Am_dyn,Bp_dyn,Bm_dyn)

      ! Renormalization factor of the resonant and anti-resonant parts

!     call dynamic_Tmatrix_ZAB(eta,nBas,nC,nO,nV,nR,nS,nOO,nVV,1d0,eGT,Omega1,Omega2,rho1,rho2,OmBSE(ia), &
!                              ZAp_dyn,ZAm_dyn,ZBp_dyn,ZBm_dyn)

      ZDyn(ia)  = dot_product(X,matmul(ZAp_dyn,X)) &
                - dot_product(Y,matmul(ZAm_dyn,Y)) &
                + dot_product(X,matmul(ZBp_dyn,Y)) & 
                - dot_product(Y,matmul(ZBm_dyn,X))  

      OmDyn(ia) = dot_product(X,matmul(Ap_dyn,X)) &
                - dot_product(Y,matmul(Am_dyn,Y)) &
                + dot_product(X,matmul(Bp_dyn,Y)) & 
                - dot_product(Y,matmul(Bm_dyn,X))  

    end if

    ZDyn(ia)  = 1d0/(1d0 - ZDyn(ia))
    OmDyn(ia) = ZDyn(ia)*OmDyn(ia)

    write(*,'(2X,I5,5X,F15.6,5X,F15.6,5X,F15.6,5X,F15.6)') & 
      ia,OmBSE(ia)*HaToeV,(OmBSE(ia)+OmDyn(ia))*HaToeV,OmDyn(ia)*HaToeV,ZDyn(ia)

  end do
  write(*,*) '---------------------------------------------------------------------------------------------------'
  write(*,*) 

end subroutine Bethe_Salpeter_Tmatrix_dynamic_perturbation
