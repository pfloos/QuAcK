subroutine Bethe_Salpeter_dynamic_perturbation_iterative(TDA,eta,nBas,nC,nO,nV,nR,nS,eGW,OmRPA,OmBSE,XpY,XmY,rho)

! Compute self-consistently the dynamical effects via perturbation theory for BSE

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: TDA
  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS

  double precision,intent(in)   :: eGW(nBas)
  double precision,intent(in)   :: OmRPA(nS)
  double precision,intent(in)   :: OmBSE(nS)
  double precision,intent(in)   :: XpY(nS,nS)
  double precision,intent(in)   :: XmY(nS,nS)
  double precision,intent(in)   :: rho(nBas,nBas,nS)

! Local variables

  logical                       :: TDA_dyn = .true.
  integer                       :: ia
  integer,parameter             :: maxS = 10
  double precision              :: gapGW

  integer                       :: nSCF
  integer                       :: maxSCF = 10
  double precision              :: Conv
  double precision              :: thresh = 1d-5


  double precision,allocatable  :: OmDyn(:)
  double precision,allocatable  :: OmOld(:)
  double precision,allocatable  :: X(:)
  double precision,allocatable  :: Y(:)
  double precision,allocatable  :: A_dyn(:,:)
  double precision,allocatable  :: B_dyn(:,:)
  double precision,allocatable  :: ZA_dyn(:,:)
  double precision,allocatable  :: ZB_dyn(:,:)

! Memory allocation

  allocate(OmDyn(nS),OmOld(nS),X(nS),Y(nS),A_dyn(nS,nS))

  if(.not.TDA_dyn) allocate(B_dyn(nS,nS))

  gapGW = eGW(nO+1) - eGW(nO) 

  Conv = 1d0
  nSCF = 0
  OmOld(:) = OmBSE(:)

  write(*,*) '---------------------------------------------------------------------------------------------------'
  write(*,*) ' First-order dynamical correction to static Bethe-Salpeter excitation energies                                 '
  write(*,*) '---------------------------------------------------------------------------------------------------'
  write(*,*)

  do while(Conv > thresh .and. nSCF < maxSCF)

    nSCF = nSCF + 1

    write(*,*) '---------------------------------------------------------------------------------------------------'
    write(*,'(2X,A15,I3)') 'Iteration n.',nSCF
    write(*,*) '---------------------------------------------------------------------------------------------------'
    write(*,'(2X,A5,1X,A20,1X,A20,1X,A20)') '#','Static (eV)','Dynamic (eV)','Correction (eV)'
    write(*,*) '---------------------------------------------------------------------------------------------------'

    do ia=1,min(nS,maxS)

 
      X(:) = 0.5d0*(XpY(ia,:) + XmY(ia,:))
      Y(:) = 0.5d0*(XpY(ia,:) - XmY(ia,:))
 
      ! Resonant part of the BSE correction
 
      call Bethe_Salpeter_A_matrix_dynamic(eta,nBas,nC,nO,nV,nR,nS,1d0,eGW(:),OmRPA(:),OmOld(ia),rho(:,:,:),A_dyn(:,:))
 
      ! First-order correction 
 
      if(TDA_dyn) then 
 
        OmDyn(ia) = dot_product(X(:),matmul(A_dyn(:,:),X(:)))
 
      else
 
        ! Anti-resonant part of the BSE correction
 
        call Bethe_Salpeter_B_matrix_dynamic(eta,nBas,nC,nO,nV,nR,nS,1d0,eGW(:),OmRPA(:),OmOld(ia),rho(:,:,:),B_dyn(:,:))
 
        OmDyn(ia) = dot_product(X(:),matmul(A_dyn(:,:),X(:))) &
                  - dot_product(Y(:),matmul(A_dyn(:,:),Y(:))) &
                  + dot_product(X(:),matmul(B_dyn(:,:),Y(:))) & 
                  - dot_product(Y(:),matmul(B_dyn(:,:),X(:)))  
 
      end if
 
      write(*,'(2X,I5,5X,F15.6,5X,F15.6,5X,F15.6)') & 
        ia,OmBSE(ia)*HaToeV,(OmBSE(ia)+OmDyn(ia))*HaToeV,OmDyn(ia)*HaToeV
 
      if(OmBSE(ia) > gapGW) write(*,*) ' !!! BSE neutral excitation larger than the GW gap !!! '
 
    end do

    Conv = maxval(abs(OmBSE(:) + OmDyn(:) - OmOld(:)))
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

end subroutine Bethe_Salpeter_dynamic_perturbation_iterative
