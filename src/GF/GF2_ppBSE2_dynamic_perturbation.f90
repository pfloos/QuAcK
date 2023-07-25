subroutine GF2_ppBSE2_dynamic_perturbation(ispin,dTDA,eta,nBas,nC,nO,nV,nR,nS,nOO,nVV,eGF,ERI,dipole_int, & 
                                           Om1,X1,Y1,Om2,X2,Y2,KB_sta,KC_sta,KD_sta)

! Compute dynamical effects via perturbation theory for BSE

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
  integer,intent(in)            :: nOO
  integer,intent(in)            :: nVV

  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: eGF(nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)
  double precision,intent(in)   :: Om1(nVV)
  double precision,intent(in)   :: X1(nVV,nVV)
  double precision,intent(in)   :: Y1(nOO,nVV)
  double precision,intent(in)   :: Om2(nOO)
  double precision,intent(in)   :: X2(nVV,nOO)
  double precision,intent(in)   :: Y2(nOO,nOO)
  double precision,intent(in)   :: KB_sta(nVV,nOO)
  double precision,intent(in)   :: KC_sta(nVV,nVV)
  double precision,intent(in)   :: KD_sta(nOO,nOO)

! Local variables

  integer                       :: ab,ij

  integer                       :: maxOO = 10
  integer                       :: maxVV = 10

  double precision,allocatable  :: Om1Dyn(:)
  double precision,allocatable  :: Om2Dyn(:)
  double precision,allocatable  :: Z1Dyn(:)
  double precision,allocatable  :: Z2Dyn(:)

  double precision,allocatable  :: KB_dyn(:,:)
  double precision,allocatable  :: KC_dyn(:,:)
  double precision,allocatable  :: KD_dyn(:,:)
  double precision,allocatable  :: ZC_dyn(:,:)
  double precision,allocatable  :: ZD_dyn(:,:)

! Memory allocation

  allocate(Om1Dyn(maxVV),Om2Dyn(maxOO),Z1Dyn(maxVV),Z2Dyn(maxOO), & 
           KB_dyn(nVV,nOO),KC_dyn(nVV,nVV),KD_dyn(nOO,nOO),   & 
           ZC_dyn(nVV,nVV),ZD_dyn(nOO,nOO))

  if(dTDA) then 
    write(*,*)
    write(*,*) '*** dynamical TDA activated ***'
    write(*,*)
  end if

  write(*,*) '---------------------------------------------------------------------------------------------------'
  write(*,*) ' First-order dynamical correction to static ppBSE2 double electron attachment energies             '
  write(*,*) '---------------------------------------------------------------------------------------------------'
  write(*,'(2X,A5,1X,A20,1X,A20,1X,A20,1X,A20)') '#','Static (eV)','Dynamic (eV)','Correction (eV)','Renorm. (eV)'
  write(*,*) '---------------------------------------------------------------------------------------------------'

  do ab=1,min(nVV,maxVV)

!   if(.not.dTDA) call GF2_ppBSE2_dynamic_kernel_B(eta,nBas,nC,nO,nV,nR,nS,nOO,nVV,1d0,eGF,OmBSE(ab),KB_dyn)
    call GF2_ppBSE2_dynamic_kernel_C(ispin,eta,nBas,nC,nO,nV,nR,nS,nVV,1d0,eGF,Om1(ab),KC_dyn,ZC_dyn)

    Z1Dyn(ab)  = dot_product(X1(:,ab),matmul(ZC_dyn,X1(:,ab)))
    Om1Dyn(ab) = dot_product(X1(:,ab),matmul(KC_dyn - KC_sta,X1(:,ab)))

    Z1Dyn(ab)  = 1d0/(1d0 - Z1Dyn(ab))
    Om1Dyn(ab) = Z1Dyn(ab)*Om1Dyn(ab)

    write(*,'(2X,I5,5X,F15.6,5X,F15.6,5X,F15.6,5X,F15.6)') & 
      ab,Om1(ab)*HaToeV,(Om1(ab)+Om1Dyn(ab))*HaToeV,Om1Dyn(ab)*HaToeV,Z1Dyn(ab)

  end do
  write(*,*) '---------------------------------------------------------------------------------------------------'
  write(*,*) 

  write(*,*) '---------------------------------------------------------------------------------------------------'
  write(*,*) ' First-order dynamical correction to static ppBSE2 double electron detachment energies             '
  write(*,*) '---------------------------------------------------------------------------------------------------'
  write(*,'(2X,A5,1X,A20,1X,A20,1X,A20,1X,A20)') '#','Static (eV)','Dynamic (eV)','Correction (eV)','Renorm. (eV)'
  write(*,*) '---------------------------------------------------------------------------------------------------'

  do ij=1,min(nOO,maxOO)

!   if(.not.dTDA) call GF2_ppBSE2_dynamic_kernel_B(eta,nBas,nC,nO,nV,nR,nS,nOO,nVV,1d0,eGF,OmBSE(ab),KB_dyn)
    call GF2_ppBSE2_dynamic_kernel_D(ispin,eta,nBas,nC,nO,nV,nR,nS,nOO,1d0,eGF,Om2(ij),KD_dyn,ZD_dyn)

    Z2Dyn(ij)  = dot_product(Y2(:,ij),matmul(ZD_dyn,Y2(:,ij)))
    Om2Dyn(ij) = dot_product(Y2(:,ij),matmul(KD_dyn - KD_sta,Y2(:,ij)))

    Z2Dyn(ij)  = 1d0/(1d0 - Z2Dyn(ij))
    Om2Dyn(ij) = Z2Dyn(ij)*Om2Dyn(ij)

    write(*,'(2X,I5,5X,F15.6,5X,F15.6,5X,F15.6,5X,F15.6)') & 
      ij,Om2(ij)*HaToeV,(Om2(ij)+Om2Dyn(ij))*HaToeV,Om2Dyn(ij)*HaToeV,Z2Dyn(ij)

  end do
  write(*,*) '---------------------------------------------------------------------------------------------------'
  write(*,*) 

end subroutine 
