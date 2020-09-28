subroutine BSE2_dynamic_perturbation_iterative(dTDA,ispin,eta,nBas,nC,nO,nV,nR,nS,ERI,dipole_int, & 
                                               eHF,eGF,OmBSE,XpY,XmY)

! Compute self-consistently the dynamical effects via perturbation theory for BSE2

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
  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: eGF(nBas)
  double precision,intent(in)   :: OmBSE(nS)
  double precision,intent(in)   :: XpY(nS,nS)
  double precision,intent(in)   :: XmY(nS,nS)

! Local variables

  integer                       :: ia

  integer,parameter             :: maxS = 10
  double precision              :: gapGF

  integer                       :: nSCF
  integer                       :: maxSCF = 10
  double precision              :: Conv
  double precision              :: thresh = 1d-3


  double precision,allocatable  :: OmDyn(:)
  double precision,allocatable  :: ZDyn(:)
  double precision,allocatable  :: OmOld(:)
  double precision,allocatable  :: X(:)
  double precision,allocatable  :: Y(:)

  double precision,allocatable  ::  Ap_dyn(:,:)
  double precision,allocatable  ::  Am_dyn(:,:)
  double precision,allocatable  :: ZAp_dyn(:,:)
  double precision,allocatable  :: ZAm_dyn(:,:)

  double precision,allocatable  ::  B_dyn(:,:)

! Memory allocation

  allocate(OmDyn(nS),OmOld(nS),ZDyn(nS),X(nS),Y(nS),Ap_dyn(nS,nS),ZAp_dyn(nS,nS))

  if(.not.dTDA) allocate(Am_dyn(nS,nS),ZAm_dyn(nS,nS),B_dyn(nS,nS))

! Print main components of transition vectors

  call print_transition_vectors(.false.,nBas,nC,nO,nV,nR,nS,OmBSE,XpY,XmY)

  if(dTDA) then
    write(*,*)
    write(*,*) '*** dynamical TDA activated ***'
    write(*,*)
  end if

  gapGF = eGF(nO+1) - eGF(nO) 

  Conv = 1d0
  nSCF = 0
  OmOld(:) = OmBSE(:)

  write(*,*) '---------------------------------------------------------------------------------------------------'
  write(*,*) ' First-order dynamical correction to static Bethe-Salpeter excitation energies                     '
  write(*,*) '---------------------------------------------------------------------------------------------------'
  write(*,'(A57,F10.6,A3)') ' BSE2 neutral excitation must be lower than the GF gap = ',gapGF*HaToeV,' eV'
  write(*,*) '---------------------------------------------------------------------------------------------------'
  write(*,*)

  do while(Conv > thresh .and. nSCF < maxSCF)

    nSCF = nSCF + 1

    write(*,*) '---------------------------------------------------------------------------------------------------'
    write(*,'(2X,A15,I3)') 'Iteration n.',nSCF
    write(*,*) '---------------------------------------------------------------------------------------------------'
    write(*,'(2X,A5,1X,A20,1X,A20,1X,A20,A20)') '#','Static (eV)','Dynamic (eV)','Correction (eV)','Convergence (eV)'
    write(*,*) '---------------------------------------------------------------------------------------------------'

    do ia=1,min(nS,maxS)

 
      X(:) = 0.5d0*(XpY(ia,:) + XmY(ia,:))
      Y(:) = 0.5d0*(XpY(ia,:) - XmY(ia,:))
 
     ! Resonant part of the BSE correction
      call BSE2_A_matrix_dynamic(ispin,eta,nBas,nC,nO,nV,nR,nS,1d0,ERI,eGF,+OmOld(ia),Ap_dyn,ZAp_dyn)

      if(dTDA) then 
 
        OmDyn(ia) = dot_product(X,matmul(Ap_dyn,X))
        ZDyn(ia)  = dot_product(X,matmul(ZAp_dyn,X))
 
      else
 
        ! Anti-resonant part of the BSE correction
 
      call BSE2_A_matrix_dynamic(ispin,eta,nBas,nC,nO,nV,nR,nS,1d0,ERI,eGF,-OmOld(ia),Am_dyn,ZAm_dyn)

      call BSE2_B_matrix_dynamic(ispin,eta,nBas,nC,nO,nV,nR,nS,1d0,ERI,eGF,B_dyn)

      ZDyn(ia)  = dot_product(X,matmul(ZAp_dyn,X)) &
                + dot_product(Y,matmul(ZAm_dyn,Y))  

      OmDyn(ia) = dot_product(X,matmul(Ap_dyn,X)) &
                - dot_product(Y,matmul(Am_dyn,Y)) &
                + dot_product(X,matmul(B_dyn,Y))  &
                - dot_product(Y,matmul(B_dyn,X))

      end if
 
      write(*,'(2X,I5,5X,F15.6,5X,F15.6,5X,F15.6,5X,F15.6)') & 
        ia,OmBSE(ia)*HaToeV,(OmBSE(ia)+OmDyn(ia))*HaToeV,OmDyn(ia)*HaToeV,(OmBSE(ia) + OmDyn(ia) - OmOld(ia))*HaToeV
 
    end do

    Conv = maxval(abs(OmBSE(:) + OmDyn(:) - OmOld(:)))*HaToeV
    OmOld(:) = OmBSE(:) + OmDyn(:)

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

end subroutine BSE2_dynamic_perturbation_iterative
