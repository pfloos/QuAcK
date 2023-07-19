subroutine GW_phBSE_dynamic_perturbation_iterative(dTDA,eta,nBas,nC,nO,nV,nR,nS,eGW,dipole_int,OmRPA,rho_RPA,OmBSE,XpY,XmY)

! Compute self-consistently the dynamical effects via perturbation theory for BSE

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

  double precision,intent(in)   :: eGW(nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)
  double precision,intent(in)   :: OmRPA(nS)
  double precision,intent(in)   :: rho_RPA(nBas,nBas,nS)
  double precision,intent(in)   :: OmBSE(nS)
  double precision,intent(in)   :: XpY(nS,nS)
  double precision,intent(in)   :: XmY(nS,nS)

! Local variables

  integer                       :: ia

  integer                       :: maxS = 10
  double precision              :: gapGW

  integer                       :: nSCF
  integer                       :: maxSCF = 10
  double precision              :: Conv
  double precision              :: thresh = 1d-3


  double precision,allocatable  :: OmDyn(:)
  double precision,allocatable  :: OmOld(:)
  double precision,allocatable  :: X(:)
  double precision,allocatable  :: Y(:)

  double precision,allocatable  ::  Ap_dyn(:,:)
  double precision,allocatable  :: ZAp_dyn(:,:)

  double precision,allocatable  ::  Bp_dyn(:,:)

  double precision,allocatable  ::  Am_dyn(:,:)
  double precision,allocatable  :: ZAm_dyn(:,:)

  double precision,allocatable  ::  Bm_dyn(:,:)

! Memory allocation

  maxS = min(nS,maxS)
  allocate(OmDyn(maxS),OmOld(maxS),X(nS),Y(nS),Ap_dyn(nS,nS),ZAp_dyn(nS,nS))
  allocate(Am_dyn(nS,nS),ZAm_dyn(nS,nS),Bp_dyn(nS,nS),Bm_dyn(nS,nS))

  if(dTDA) then
    write(*,*)
    write(*,*) '*** dynamical TDA activated ***'
    write(*,*)
  end if

  gapGW = eGW(nO+1) - eGW(nO) 

  Conv = 1d0
  nSCF = 0
  OmOld(1:maxS) = OmBSE(1:maxS)

  write(*,*) '---------------------------------------------------------------------------------------------------'
  write(*,*) ' First-order dynamical correction to static Bethe-Salpeter excitation energies                     '
  write(*,*) '---------------------------------------------------------------------------------------------------'
  write(*,'(A57,F10.6,A3)') ' BSE neutral excitation must be lower than the GW gap = ',gapGW*HaToeV,' eV'
  write(*,*) '---------------------------------------------------------------------------------------------------'
  write(*,*)

  do while(Conv > thresh .and. nSCF < maxSCF)

    nSCF = nSCF + 1

    write(*,*) '---------------------------------------------------------------------------------------------------'
    write(*,'(2X,A15,I3)') 'Iteration n.',nSCF
    write(*,*) '---------------------------------------------------------------------------------------------------'
    write(*,'(2X,A5,1X,A20,1X,A20,1X,A20,A20)') '#','Static (eV)','Dynamic (eV)','Correction (eV)','Convergence (eV)'
    write(*,*) '---------------------------------------------------------------------------------------------------'

    do ia=1,maxS

      X(:) = 0.5d0*(XpY(ia,:) + XmY(ia,:))
      Y(:) = 0.5d0*(XpY(ia,:) - XmY(ia,:))
 
      ! First-order correction 
 
      if(dTDA) then 

       ! Resonant part of the BSE correction
 
        call GW_phBSE_dynamic_kernel_A(eta,nBas,nC,nO,nV,nR,nS,1d0,eGW,OmRPA,rho_RPA,OmOld(ia),Ap_dyn,ZAp_dyn)
 
        OmDyn(ia) = dot_product(X(:),matmul(Ap_dyn(:,:),X(:)))
 
      else
 
        ! Anti-resonant part of the BSE correction
 
        call GW_phBSE_dynamic_kernel(eta,nBas,nC,nO,nV,nR,nS,1d0,eGW,OmRPA,rho_RPA,OmOld(ia),Ap_dyn,Am_dyn,Bp_dyn,Bm_dyn)
 
        OmDyn(ia) = dot_product(X(:),matmul(Ap_dyn(:,:),X(:))) &
                  - dot_product(Y(:),matmul(Am_dyn(:,:),Y(:))) &
                  + dot_product(X(:),matmul(Bp_dyn(:,:),Y(:))) &
                  - dot_product(Y(:),matmul(Bm_dyn(:,:),X(:)))

      end if
 
      write(*,'(2X,I5,5X,F15.6,5X,F15.6,5X,F15.6,5X,F15.6)') & 
        ia,OmBSE(ia)*HaToeV,(OmBSE(ia)+OmDyn(ia))*HaToeV,OmDyn(ia)*HaToeV,(OmBSE(ia) + OmDyn(ia) - OmOld(ia))*HaToeV
 
    end do

    Conv = maxval(abs(OmBSE(1:maxS) + OmDyn(:) - OmOld(:)))*HaToeV
    OmOld(:) = OmBSE(1:maxS) + OmDyn(:)

    write(*,*) '---------------------------------------------------------------------------------------------------'
    write(*,'(2X,A20,1X,F10.6)') ' Convergence = ',Conv
    write(*,*) '---------------------------------------------------------------------------------------------------'
    write(*,*) 

  end do

! Did it actually converge?

  if(nSCF == maxSCF) then

    write(*,*)
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)'                 Convergence failed                 '
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)

  endif

end subroutine 
