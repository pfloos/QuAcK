subroutine Bethe_Salpeter_Tmatrix_so(eta,nBas,nC,nO,nV,nR,nS,nOO,nVV,Omega1,X1,Y1,Omega2,X2,Y2,rho1,rho2, &
                                     ERI,eT,eGT,EcBSE)

! Compute the Bethe-Salpeter excitation energies with the T-matrix kernel

  implicit none
  include 'parameters.h'

! Input variables

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
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

  double precision,intent(in)   :: Omega1(nVV)
  double precision,intent(in)   :: X1(nVV,nVV)
  double precision,intent(in)   :: Y1(nOO,nVV)
  double precision,intent(in)   :: Omega2(nOO)
  double precision,intent(in)   :: X2(nVV,nOO)
  double precision,intent(in)   :: Y2(nOO,nOO)
  double precision,intent(in)   :: rho1(nBas,nBas,nVV)
  double precision,intent(in)   :: rho2(nBas,nBas,nOO)

! Local variables

  integer                       :: ispin

  double precision              :: EcRPA
  double precision,allocatable  :: TA(:,:),TB(:,:)
  double precision,allocatable  :: OmBSE(:)
  double precision,allocatable  :: XpY_BSE(:,:)
  double precision,allocatable  :: XmY_BSE(:,:)

! Output variables

  double precision,intent(out)  :: EcBSE

! Memory allocation

  allocate(TA(nS,nS),TB(nS,nS),OmBSE(nS),XpY_BSE(nS,nS),XmY_BSE(nS,nS))

!------------------!
! Compute T-matrix !
!------------------!

  ispin = 4

  call linear_response_pp(ispin,.false.,nBas,nC,nO,nV,nR,nOO,nVV,1d0,eT,ERI, &
                          Omega1,X1,Y1,Omega2,X2,Y2,EcRPA)

  call static_Tmatrix_A(eta,nBas,nC,nO,nV,nR,nS,nOO,nVV,1d0,Omega1,rho1,Omega2,rho2,TA)
  call static_Tmatrix_B(eta,nBas,nC,nO,nV,nR,nS,nOO,nVV,1d0,Omega1,rho1,Omega2,rho2,TB)

!------------------!
! Singlet manifold !
!------------------!

    ispin = 3
    EcBSE = 0d0

    ! Compute BSE singlet excitation energies

    call linear_response_BSE(ispin,.false.,.false.,eta,nBas,nC,nO,nV,nR,nS,1d0,eGT,ERI,TA,TB, &
                             EcBSE,OmBSE,XpY_BSE,XmY_BSE)

    call print_excitation('BSE@GT      ',ispin,nS,OmBSE)

end subroutine Bethe_Salpeter_Tmatrix_so
