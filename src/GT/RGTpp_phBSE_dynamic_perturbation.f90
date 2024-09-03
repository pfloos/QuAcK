subroutine RGTpp_phBSE_dynamic_perturbation(ispin,dTDA,eta,nBas,nC,nO,nV,nR,nS,nOOab,nVVab,nOOaa,nVVaa,Om1ab,Om2ab,Om1aa,Om2aa, & 
                                           rho1ab,rho2ab,rho1aa,rho2aa,eT,eGT,dipole_int,OmBSE,XpY,XmY,TAab,TAaa)

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

  integer,intent(in)            :: nOOab
  integer,intent(in)            :: nVVab
  integer,intent(in)            :: nOOaa
  integer,intent(in)            :: nVVaa

  double precision,intent(in)   :: eT(nBas)
  double precision,intent(in)   :: eGT(nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)
  double precision,intent(in)   :: OmBSE(nS)
  double precision,intent(in)   :: XpY(nS,nS)
  double precision,intent(in)   :: XmY(nS,nS)

  double precision,intent(in)   :: Om1ab(nVVab)
  double precision,intent(in)   :: Om2ab(nOOab)
  double precision,intent(in)   :: rho1ab(nBas,nBas,nVVab)
  double precision,intent(in)   :: rho2ab(nBas,nBas,nOOab)
  double precision,intent(in)   :: Om1aa(nVVaa)
  double precision,intent(in)   :: Om2aa(nOOaa)
  double precision,intent(in)   :: rho1aa(nBas,nBas,nVVaa)
  double precision,intent(in)   :: rho2aa(nBas,nBas,nOOaa)

  double precision,intent(in)   :: TAab(nS,nS)
  double precision,intent(in)   :: TAaa(nS,nS)

! Local variables

  integer                       :: ia

  integer                       :: maxS = 10
  double precision              :: gapGT

  double precision,allocatable  :: OmDyn(:)
  double precision,allocatable  :: ZDyn(:)
  double precision,allocatable  :: X(:)
  double precision,allocatable  :: Y(:)

  double precision,allocatable  :: dTAab(:,:)
  double precision,allocatable  :: ZAab(:,:)

  double precision,allocatable  :: dTAaa(:,:)
  double precision,allocatable  :: ZAaa(:,:)

! Memory allocation

  maxS = min(nS,maxS)
  allocate(OmDyn(maxS),ZDyn(maxS),X(nS),Y(nS),dTAab(nS,nS),ZAab(nS,nS),dTAaa(nS,nS),ZAaa(nS,nS))

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

    call RGTpp_phBSE_dynamic_kernel_A(eta,nBas,nC,nO,nV,nR,nS,nOOab,nVVab,1d0,eGT,Om1ab,Om2ab,rho1ab,rho2ab,OmBSE(ia),dTAab,ZAab)
 
    ! Compute dynamical T-matrix for alpha-beta block  

    call RGTpp_phBSE_dynamic_kernel_A(eta,nBas,nC,nO,nV,nR,nS,nOOaa,nVVaa,1d0,eGT,Om1aa,Om2aa,rho1aa,rho2aa,OmBSE(ia),dTAaa,ZAaa)
 
    X(:) = 0.5d0*(XpY(ia,:) + XmY(ia,:))
    Y(:) = 0.5d0*(XpY(ia,:) - XmY(ia,:))

    ! First-order correction 
    
    if(ispin == 1) then 
      ZDyn(ia)  = - dot_product(X,matmul(ZAaa+ZAab,X)) 
      OmDyn(ia) = - dot_product(X,matmul(dTAaa+dTAab,X)) + dot_product(X,matmul(TAaa+TAab,X)) 
    end if

    if(ispin == 2) then 
      ZDyn(ia)  = - dot_product(X,matmul(ZAaa-ZAab,X)) 
      OmDyn(ia) = - dot_product(X,matmul(dTAaa-dTAab,X)) + dot_product(X,matmul(TAaa-TAab,X)) 
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

end subroutine 
