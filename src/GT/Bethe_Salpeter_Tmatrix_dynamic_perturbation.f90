subroutine Bethe_Salpeter_Tmatrix_dynamic_perturbation(singlet,triplet,dTDA,eta,nBas,nC,nO,nV,nR,nS,nOOs,nVVs,nOOt,nVVt, &
                                                       Omega1s,Omega2s,Omega1t,Omega2t,rho1s,rho2s,rho1t,rho2t,eT,eGT,   &
                                                       dipole_int,OmBSE,XpY,XmY,TAs,TBs,TAt,TBt)

! Compute dynamical effects via perturbation theory for BSE@GT

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: singlet
  logical,intent(in)            :: triplet
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
  double precision,intent(in)   :: OmBSE(nS,nspin)
  double precision,intent(in)   :: XpY(nS,nS,nspin)
  double precision,intent(in)   :: XmY(nS,nS,nspin)

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
  integer                       :: ispin

  integer                       :: maxS = 10
  double precision              :: gapGT

  double precision,allocatable  :: OmDyn(:,:)
  double precision,allocatable  :: ZDyn(:,:)
  double precision,allocatable  :: X(:)
  double precision,allocatable  :: Y(:)

  double precision,allocatable  :: dTAs(:,:)
  double precision,allocatable  :: ZAs(:,:)

  double precision,allocatable  :: dTAt(:,:)
  double precision,allocatable  :: ZAt(:,:)

! Memory allocation

  maxS = min(nS,maxS)
  allocate(OmDyn(maxS,nspin),ZDyn(maxS,nspin),X(nS),Y(nS),dTAs(nS,nS),ZAs(nS,nS),dTAt(nS,nS),ZAt(nS,nS))

  if(dTDA) then 
    write(*,*)
    write(*,*) '*** dynamical TDA activated ***'
    write(*,*)
  else
    print*, ' Beyond-TDA dynamical correction for BSE@GT NYI'
    return
  end if

  OmDyn(:,:) = 0d0
  ZDyn(:,:)  = 0d0

  do ia=1,maxS

  ! Compute dynamical T-matrix for alpha-beta block !

    ispin = 1
    call dynamic_Tmatrix_A(eta,nBas,nC,nO,nV,nR,nS,nOOs,nVVs,1d0,eGT,Omega1s,Omega2s,rho1s,rho2s,OmBSE(ia,ispin),dTAs,ZAs)
 
  ! Compute dynamical T-matrix for alpha-beta block !

    ispin = 2
    call dynamic_Tmatrix_A(eta,nBas,nC,nO,nV,nR,nS,nOOt,nVVt,1d0,eGT,Omega1t,Omega2t,rho1t,rho2t,OmBSE(ia,ispin),dTAt,ZAt)

    do ispin=1,nspin
 
      X(:) = 0.5d0*(XpY(ia,:,ispin) + XmY(ia,:,ispin))
      Y(:) = 0.5d0*(XpY(ia,:,ispin) - XmY(ia,:,ispin))

      ! First-order correction 
     
      if(ispin == 1) then 
        ZDyn(ia,ispin)  = dot_product(X,matmul(ZAt+ZAs,X)) 
        OmDyn(ia,ispin) = dot_product(X,matmul(dTAt+dTAs,X)) - dot_product(X,matmul(TAt+TAs,X))
      end if

      if(ispin == 2) then 
        ZDyn(ia,ispin)  = dot_product(X,matmul(ZAt-ZAs,X)) 
        OmDyn(ia,ispin) = dot_product(X,matmul(dTAt-dTAs,X)) - dot_product(X,matmul(TAt-TAs,X))
      end if
     
      ZDyn(ia,ispin)  = 1d0/(1d0 - ZDyn(ia,ispin))
      OmDyn(ia,ispin) = ZDyn(ia,ispin)*OmDyn(ia,ispin)

    end do

  end do

!--------------!
! Dump results !
!--------------!

  gapGT = eGT(nO+1) - eGT(nO) 

  if(singlet) then 

    ispin = 1

    write(*,*) '---------------------------------------------------------------------------------------------------'
    write(*,*) ' First-order dynamical correction to static singlet Bethe-Salpeter excitation energies             '
    write(*,*) '---------------------------------------------------------------------------------------------------'
    write(*,'(A57,F10.6,A3)') ' BSE neutral excitation must be lower than the GT gap = ',gapGT*HaToeV,' eV'
    write(*,*) '---------------------------------------------------------------------------------------------------'
    write(*,'(2X,A5,1X,A20,1X,A20,1X,A20,1X,A20)') '#','Static (eV)','Dynamic (eV)','Correction (eV)','Renorm. (eV)'
    write(*,*) '---------------------------------------------------------------------------------------------------'

    do ia=1,maxS
      write(*,'(2X,I5,5X,F15.6,5X,F15.6,5X,F15.6,5X,F15.6)') & 
        ia,OmBSE(ia,ispin)*HaToeV,(OmBSE(ia,ispin)+OmDyn(ia,ispin))*HaToeV,OmDyn(ia,ispin)*HaToeV,ZDyn(ia,ispin)
    end do

    write(*,*) '---------------------------------------------------------------------------------------------------'
    write(*,*) 

  end if

  if(triplet) then 

    ispin = 2

    write(*,*) '---------------------------------------------------------------------------------------------------'
    write(*,*) ' First-order dynamical correction to static triplet Bethe-Salpeter excitation energies             '
    write(*,*) '---------------------------------------------------------------------------------------------------'
    write(*,'(A57,F10.6,A3)') ' BSE neutral excitation must be lower than the GT gap = ',gapGT*HaToeV,' eV'
    write(*,*) '---------------------------------------------------------------------------------------------------'
    write(*,'(2X,A5,1X,A20,1X,A20,1X,A20,1X,A20)') '#','Static (eV)','Dynamic (eV)','Correction (eV)','Renorm. (eV)'
    write(*,*) '---------------------------------------------------------------------------------------------------'

    do ia=1,maxS
      write(*,'(2X,I5,5X,F15.6,5X,F15.6,5X,F15.6,5X,F15.6)') & 
        ia,OmBSE(ia,ispin)*HaToeV,(OmBSE(ia,ispin)+OmDyn(ia,ispin))*HaToeV,OmDyn(ia,ispin)*HaToeV,ZDyn(ia,ispin)
    end do

    write(*,*) '---------------------------------------------------------------------------------------------------'
    write(*,*) 

  end if

end subroutine Bethe_Salpeter_Tmatrix_dynamic_perturbation
