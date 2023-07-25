subroutine GF2_phBSE2_dynamic_perturbation(dTDA,ispin,eta,nBas,nC,nO,nV,nR,nS,ERI,dipole_int,eGF,KA_sta,KB_sta,OmBSE,XpY,XmY)

! Compute dynamical effects via perturbation theory for BSE

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dTDA
  integer,intent(in)            :: ispin
  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS

  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)
  double precision,intent(in)   :: eGF(nBas)
  double precision,intent(in)   :: KA_sta(nS,nS)
  double precision,intent(in)   :: KB_sta(nS,nS)
  double precision,intent(in)   :: OmBSE(nS)
  double precision,intent(in)   :: XpY(nS,nS)
  double precision,intent(in)   :: XmY(nS,nS)

! Local variables

  integer                       :: ia
  integer,parameter             :: maxS = 10
  double precision              :: gapGF

  double precision,allocatable  :: OmDyn(:)
  double precision,allocatable  :: ZDyn(:)
  double precision,allocatable  :: X(:)
  double precision,allocatable  :: Y(:)

  double precision,allocatable  :: KAp_dyn(:,:)
  double precision,allocatable  :: KAm_dyn(:,:)
  double precision,allocatable  :: ZAp_dyn(:,:)
  double precision,allocatable  :: ZAm_dyn(:,:)

  double precision,allocatable  :: KB_dyn(:,:)

! Memory allocation

  allocate(OmDyn(nS),ZDyn(nS),X(nS),Y(nS),KAp_dyn(nS,nS),ZAp_dyn(nS,nS), &
           KAm_dyn(nS,nS),ZAm_dyn(nS,nS),KB_dyn(nS,nS))

  if(dTDA) then
    write(*,*)
    write(*,*) '*** dynamical TDA activated ***'
    write(*,*)
  end if

  ! Print main components of transition vectors

  call print_transition_vectors_ph(.false.,nBas,nC,nO,nV,nR,nS,dipole_int,OmBSE,XpY,XmY)

  gapGF = eGF(nO+1) - eGF(nO) 

  write(*,*) '---------------------------------------------------------------------------------------------------'
  write(*,*) ' First-order dynamical correction to static 2nd-order Bethe-Salpeter excitation energies '
  write(*,*) '---------------------------------------------------------------------------------------------------'
  write(*,'(A58,F10.6,A3)') ' BSE neutral excitation must be lower than the GF2 gap = ',gapGF*HaToeV,' eV'
  write(*,*) '---------------------------------------------------------------------------------------------------'
  write(*,'(2X,A5,1X,A20,1X,A20,1X,A20,1X,A20)') '#','Static (eV)','Dynamic (eV)','Correction (eV)','Renorm. (eV)'
  write(*,*) '---------------------------------------------------------------------------------------------------'

  do ia=1,min(nS,maxS)

    X(:) = 0.5d0*(XpY(ia,:) + XmY(ia,:))
    Y(:) = 0.5d0*(XpY(ia,:) - XmY(ia,:))

    ! Resonant part of the BSE correction for dynamical TDA

    call GF2_phBSE2_dynamic_kernel_A(ispin,eta,nBas,nC,nO,nV,nR,nS,1d0,ERI,eGF,+OmBSE(ia),KAp_dyn,ZAp_dyn)

    if(dTDA) then 

      ZDyn(ia)  = dot_product(X,matmul(ZAp_dyn,X))
      OmDyn(ia) = dot_product(X,matmul(KAp_dyn - KA_sta,X))

    else

      ! Second part of the resonant and anti-resonant part of the BSE correction (frequency independent)

      call GF2_phBSE2_dynamic_kernel_A(ispin,eta,nBas,nC,nO,nV,nR,nS,1d0,ERI,eGF,-OmBSE(ia),KAm_dyn,ZAm_dyn)

      call GF2_phBSE2_dynamic_kernel_B(ispin,eta,nBas,nC,nO,nV,nR,nS,1d0,ERI,eGF,KB_dyn)

      ZDyn(ia)  = dot_product(X,matmul(ZAp_dyn,X)) &
                + dot_product(Y,matmul(ZAm_dyn,Y))  

      OmDyn(ia) = dot_product(X,matmul(KAp_dyn - KA_sta,X)) &
                - dot_product(Y,matmul(KAm_dyn - KA_sta,Y)) &
                + dot_product(X,matmul(KB_dyn  - KB_sta,Y)) & 
                - dot_product(Y,matmul(KB_dyn  - KB_sta,X))  

    end if

    ZDyn(ia) = 1d0/(1d0 - ZDyn(ia))
    OmDyn(ia) = ZDyn(ia)*OmDyn(ia)

    write(*,'(2X,I5,5X,F15.6,5X,F15.6,5X,F15.6,5X,F15.6)') & 
      ia,OmBSE(ia)*HaToeV,(OmBSE(ia)+OmDyn(ia))*HaToeV,OmDyn(ia)*HaToeV,ZDyn(ia)

  end do
  write(*,*) '---------------------------------------------------------------------------------------------------'
  write(*,*) 

end subroutine 
