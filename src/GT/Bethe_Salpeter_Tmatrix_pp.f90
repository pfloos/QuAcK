subroutine Bethe_Salpeter_Tmatrix_pp(TDA_T,TDA,dBSE,dTDA,evDyn,singlet,triplet,eta,nBas,nC,nO,nV,nR,nOOab,nVVab,nOOaa,nVVaa, &
                                     Om1ab,X1ab,Y1ab,Om2ab,X2ab,Y2ab,rho1ab,rho2ab,Om1aa,X1aa,Y1aa,Om2aa,X2aa,Y2aa,rho1aa,rho2aa, & 
                                     ERI,dipole_int,eT,eGT,EcBSE)

! Compute the Bethe-Salpeter excitation energies with the T-matrix kernel

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: TDA_T
  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA
  logical,intent(in)            :: evDyn
  logical,intent(in)            :: singlet
  logical,intent(in)            :: triplet

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR

  integer,intent(in)            :: nOOab
  integer,intent(in)            :: nOOaa
  integer,intent(in)            :: nVVab
  integer,intent(in)            :: nVVaa

  double precision,intent(in)   :: eT(nBas)
  double precision,intent(in)   :: eGT(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)

  double precision,intent(in)   :: Om1ab(nVVab)
  double precision,intent(in)   :: X1ab(nVVab,nVVab)
  double precision,intent(in)   :: Y1ab(nOOab,nVVab)
  double precision,intent(in)   :: Om2ab(nOOab)
  double precision,intent(in)   :: X2ab(nVVab,nOOab)
  double precision,intent(in)   :: Y2ab(nOOab,nOOab)
  double precision,intent(in)   :: rho1ab(nBas,nBas,nVVab)
  double precision,intent(in)   :: rho2ab(nBas,nBas,nOOab)
  double precision,intent(in)   :: Om1aa(nVVaa)
  double precision,intent(in)   :: X1aa(nVVaa,nVVaa)
  double precision,intent(in)   :: Y1aa(nOOaa,nVVaa)
  double precision,intent(in)   :: Om2aa(nOOaa)
  double precision,intent(in)   :: X2aa(nVVaa,nOOaa)
  double precision,intent(in)   :: Y2aa(nOOaa,nOOaa) 
  double precision,intent(in)   :: rho1aa(nBas,nBas,nVVaa)
  double precision,intent(in)   :: rho2aa(nBas,nBas,nOOaa)

! Local variables

  integer                       :: ispin
  integer                       :: iblock

  integer                       :: nOOs
  integer                       :: nOOt
  integer                       :: nVVs
  integer                       :: nVVt

  double precision              :: EcRPA(nspin)
  double precision,allocatable  :: TBab(:,:),TCab(:,:),TDab(:,:)
  double precision,allocatable  :: TBaa(:,:),TCaa(:,:),TDaa(:,:)

  double precision,allocatable  :: Om1s(:),Om1t(:)
  double precision,allocatable  :: X1s(:,:),X1t(:,:)
  double precision,allocatable  :: Y1s(:,:),Y1t(:,:)
  double precision,allocatable  :: Om2s(:),Om2t(:)
  double precision,allocatable  :: X2s(:,:),X2t(:,:)
  double precision,allocatable  :: Y2s(:,:),Y2t(:,:)

! Output variables

  double precision,intent(out)  :: EcBSE(nspin)


!------------------!
! Singlet manifold !
!------------------!

  if(singlet) then

    ispin  = 1

    nOOs = nO*(nO+1)/2
    nVVs = nV*(nV+1)/2

  !---------------------------------------!
  ! Compute T-matrix for alpha-beta block !
  !---------------------------------------!

    iblock = 3
 
    EcRPA(ispin) = 0d0
 
    call linear_response_pp(iblock,TDA_T,nBas,nC,nO,nV,nR,nOOab,nVVab,1d0,eT,ERI,Om1ab,X1ab,Y1ab,Om2ab,X2ab,Y2ab,EcRPA(ispin))
 
    allocate(TBab(nVVs,nOOs),TCab(nVVs,nVVs),TDab(nOOs,nOOs))
 
    if(.not.TDA) call static_Tmatrix_B_pp(ispin,eta,nBas,nC,nO,nV,nR,nOOab,nVVab,nOOs,nVVs,1d0,Om1ab,rho1ab,Om2ab,rho2ab,TBab)
                 call static_Tmatrix_C_pp(ispin,eta,nBas,nC,nO,nV,nR,nOOab,nVVab,nOOs,nVVs,1d0,Om1ab,rho1ab,Om2ab,rho2ab,TCab)
                 call static_Tmatrix_D_pp(ispin,eta,nBas,nC,nO,nV,nR,nOOab,nVVab,nOOs,nVVs,1d0,Om1ab,rho1ab,Om2ab,rho2ab,TDab)

  !----------------------------------------!
  ! Compute T-matrix for alpha-alpha block !
  !----------------------------------------!

    iblock = 4
    EcRPA(ispin) = 0d0
 
    call linear_response_pp(iblock,TDA_T,nBas,nC,nO,nV,nR,nOOaa,nVVaa,1d0,eT,ERI,Om1aa,X1aa,Y1aa,Om2aa,X2aa,Y2aa,EcRPA(ispin))

    allocate(TBaa(nVVs,nOOs),TCaa(nVVs,nVVs),TDaa(nOOs,nOOs))
 
    if(.not.TDA) call static_Tmatrix_B_pp(ispin,eta,nBas,nC,nO,nV,nR,nOOaa,nVVaa,nOOs,nVVs,1d0,Om1aa,rho1aa,Om2aa,rho2aa,TBaa)
                 call static_Tmatrix_C_pp(ispin,eta,nBas,nC,nO,nV,nR,nOOaa,nVVaa,nOOs,nVVs,1d0,Om1aa,rho1aa,Om2aa,rho2aa,TCaa)
                 call static_Tmatrix_D_pp(ispin,eta,nBas,nC,nO,nV,nR,nOOaa,nVVaa,nOOs,nVVs,1d0,Om1aa,rho1aa,Om2aa,rho2aa,TDaa)

  !----------------------------------!
  ! pp/hh sectors for singlet states !
  !----------------------------------!

    EcBSE(ispin) = 0d0

    allocate(Om1s(nVVs),X1s(nVVs,nVVs),Y1s(nOOs,nVVs),Om2s(nOOs),X2s(nVVs,nOOs),Y2s(nOOs,nOOs))
  
    call linear_response_pp_BSE(ispin,TDA,.true.,nBas,nC,nO,nV,nR,nOOs,nVVs,1d0,eGT,ERI, &
                                TBaa+TBab,TCaa+TCab,TDaa+TDab,Om1s,X1s,Y1s,Om2s,X2s,Y2s,EcBSE(ispin))

    call print_transition_vectors_pp(.true.,nBas,nC,nO,nV,nR,nOOs,nVVs,dipole_int,Om1s,X1s,Y1s,Om2s,X2s,Y2s)

    deallocate(Om1s,X1s,Y1s,Om2s,X2s,Y2s,TBab,TCab,TDab,TBaa,TCaa,TDaa)

  end if

!------------------!
! Triplet manifold !
!------------------!

  if(triplet) then

    ispin  = 2

    nOOt = nO*(nO-1)/2
    nVVt = nV*(nV-1)/2

  !---------------------------------------!
  ! Compute T-matrix for alpha-beta block !
  !---------------------------------------!

    iblock = 3
    EcRPA(ispin) = 0d0
 
    call linear_response_pp(iblock,TDA_T,nBas,nC,nO,nV,nR,nOOab,nVVab,1d0,eT,ERI,Om1ab,X1ab,Y1ab,Om2ab,X2ab,Y2ab,EcRPA(ispin))
 
    allocate(TBab(nVVt,nOOt),TCab(nVVt,nVVt),TDab(nOOt,nOOt))
 
    if(.not.TDA) call static_Tmatrix_B_pp(ispin,eta,nBas,nC,nO,nV,nR,nOOab,nVVab,nOOt,nVVt,1d0,Om1ab,rho1ab,Om2ab,rho2ab,TBab)
                 call static_Tmatrix_C_pp(ispin,eta,nBas,nC,nO,nV,nR,nOOab,nVVab,nOOt,nVVt,1d0,Om1ab,rho1ab,Om2ab,rho2ab,TCab)
                 call static_Tmatrix_D_pp(ispin,eta,nBas,nC,nO,nV,nR,nOOab,nVVab,nOOt,nVVt,1d0,Om1ab,rho1ab,Om2ab,rho2ab,TDab)
 
  !----------------------------------------!
  ! Compute T-matrix for alpha-alpha block !
  !----------------------------------------!

    ispin  = 2
    iblock = 4
 
    EcRPA(ispin) = 0d0
 
    call linear_response_pp(iblock,TDA_T,nBas,nC,nO,nV,nR,nOOaa,nVVaa,1d0,eT,ERI,Om1aa,X1aa,Y1aa,Om2aa,X2aa,Y2aa,EcRPA(ispin))

    allocate(TBaa(nVVt,nOOt),TCaa(nVVt,nVVt),TDaa(nOOt,nOOt))
 
    if(.not.TDA) call static_Tmatrix_B_pp(ispin,eta,nBas,nC,nO,nV,nR,nOOaa,nVVaa,nOOt,nVVt,1d0,Om1aa,rho1aa,Om2aa,rho2aa,TBaa)
                 call static_Tmatrix_C_pp(ispin,eta,nBas,nC,nO,nV,nR,nOOaa,nVVaa,nOOt,nVVt,1d0,Om1aa,rho1aa,Om2aa,rho2aa,TCaa)
                 call static_Tmatrix_D_pp(ispin,eta,nBas,nC,nO,nV,nR,nOOaa,nVVaa,nOOt,nVVt,1d0,Om1aa,rho1aa,Om2aa,rho2aa,TDaa)

  !----------------------------------!
  ! pp/hh sectors for triplet states !
  !----------------------------------!

    EcBSE(ispin) = 0d0

    allocate(Om1t(nVVt),X1t(nVVt,nVVt),Y1t(nOOt,nVVt),Om2t(nOOt),X2t(nVVt,nOOt),Y2t(nOOt,nOOt))
  
    call linear_response_pp_BSE(ispin,TDA,.true.,nBas,nC,nO,nV,nR,nOOt,nVVt,1d0,eGT,ERI, &
                                TBaa-TBab,TCaa-TCab,TDaa-TDab,Om1t,X1t,Y1t,Om2t,X2t,Y2t,EcBSE(ispin))

    call print_transition_vectors_pp(.false.,nBas,nC,nO,nV,nR,nOOt,nVVt,dipole_int,Om1t,X1t,Y1t,Om2t,X2t,Y2t)

    deallocate(Om1t,X1t,Y1t,Om2t,X2t,Y2t,TBab,TCab,TDab,TBaa,TCaa,TDaa)

  end if

end subroutine Bethe_Salpeter_Tmatrix_pp
