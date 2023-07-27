subroutine UGW_phBSE_dynamic_perturbation(ispin,dTDA,eta,nBas,nC,nO,nV,nR,nS,nSa,nSb,nSt,nS_sc,eW,eGW, &
                                          ERI_aaaa,ERI_aabb,ERI_bbbb,dipole_int_aa,dipole_int_bb,      & 
                                          OmRPA,rho_RPA,OmBSE,XpY_BSE,XmY_BSE)

! Compute dynamical effects via perturbation theory for BSE

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: ispin
  logical,intent(in)            :: dTDA 
  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC(nspin)
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nV(nspin)
  integer,intent(in)            :: nR(nspin)
  integer,intent(in)            :: nS(nspin)
  integer,intent(in)            :: nSa
  integer,intent(in)            :: nSb
  integer,intent(in)            :: nSt
  integer,intent(in)            :: nS_sc

  double precision,intent(in)   :: eW(nBas,nspin)
  double precision,intent(in)   :: eGW(nBas,nspin)
  double precision,intent(in)   :: ERI_aaaa(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_aabb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_bbbb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int_aa(nBas,nBas,ncart)
  double precision,intent(in)   :: dipole_int_bb(nBas,nBas,ncart)
  double precision,intent(in)   :: OmRPA(nS_sc)
  double precision,intent(in)   :: rho_RPA(nBas,nBas,nS_sc,nspin)
  double precision,intent(in)   :: OmBSE(nSt)
  double precision,intent(in)   :: XpY_BSE(nSt,nSt)
  double precision,intent(in)   :: XmY_BSE(nSt,nSt)

! Local variables

  integer                       :: ia

  integer,parameter             :: maxS = 10
  double precision              :: gapGW

  double precision,allocatable  :: OmDyn(:)
  double precision,allocatable  :: ZDyn(:)
  double precision,allocatable  :: X(:)
  double precision,allocatable  :: Y(:)

  double precision,allocatable  ::  A_dyn(:,:)
  double precision,allocatable  :: ZA_dyn(:,:)

! Memory allocation

  allocate(OmDyn(nSt),ZDyn(nSt),X(nSt),Y(nSt),A_dyn(nSt,nSt),ZA_dyn(nSt,nSt))

! Print main components of transition vectors

  if(dTDA) then 
    write(*,*)
    write(*,*) '*** dynamical TDA activated ***'
    write(*,*)
  end if

  gapGW = min(eGW(nO(1)+1,1),eGW(nO(2)+1,2)) - max(eGW(nO(1),1),eGW(nO(2),2))

  write(*,*) '---------------------------------------------------------------------------------------------------'
  write(*,*) ' First-order dynamical correction to static Bethe-Salpeter excitation energies                     '
  write(*,*) '---------------------------------------------------------------------------------------------------'
  write(*,'(A57,F10.6,A3)') ' BSE neutral excitation must be lower than the GW gap = ',gapGW*HaToeV,' eV'
  write(*,*) '---------------------------------------------------------------------------------------------------'
  write(*,'(2X,A5,1X,A20,1X,A20,1X,A20,1X,A20)') '#','Static (eV)','Dynamic (eV)','Correction (eV)','Renorm. (eV)'
  write(*,*) '---------------------------------------------------------------------------------------------------'

  do ia=1,min(nSt,maxS)

    X(:) = 0.5d0*(XpY_BSE(ia,:) + XmY_BSE(ia,:))
    Y(:) = 0.5d0*(XpY_BSE(ia,:) - XmY_BSE(ia,:))

    ! First-order correction 

    if(dTDA) then 

      ! Resonant part of the BSE correction for dynamical TDA and its renormalization factor

      call UGW_phBSE_dynamic_kernel_A(ispin,eta,nBas,nC,nO,nV,nR,nSa,nSb,nSt,nS_sc,1d0,eGW, & 
                                      ERI_aaaa,ERI_aabb,ERI_bbbb,OmRPA,rho_RPA,OmBSE(ia),A_dyn,ZA_dyn)

      ZDyn(ia)  = dot_product(X,matmul(ZA_dyn,X))
      OmDyn(ia) = dot_product(X,matmul( A_dyn,X))

    end if

    ZDyn(ia)  = 1d0/(1d0 - ZDyn(ia))
    OmDyn(ia) = ZDyn(ia)*OmDyn(ia)

    write(*,'(2X,I5,5X,F15.6,5X,F15.6,5X,F15.6,5X,F15.6)') & 
      ia,OmBSE(ia)*HaToeV,(OmBSE(ia)+OmDyn(ia))*HaToeV,OmDyn(ia)*HaToeV,ZDyn(ia)

  end do
  write(*,*) '---------------------------------------------------------------------------------------------------'
  write(*,*) 

end subroutine 
