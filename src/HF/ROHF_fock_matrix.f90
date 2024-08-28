
! ---

subroutine ROHF_fock_matrix(nBas_AOs, nBas_MOs, nOa, nOb, S, c, FaAO, FbAO, FAO)

! Construct the ROHF Fock matrix in the AO basis
! For open shells, the ROHF Fock matrix in the MO basis reads
!
!       |   F-K    |  F + K/2  |    F     |
!       |---------------------------------|
!       | F + K/2  |     F     |  F - K/2 |
!       |---------------------------------|
!       |    F     |  F - K/2  |  F + K   |
!
! with F = 1/2 (Fa + Fb) and K = Fb - Fa
!

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas_AOs, nBas_MOs
  integer,intent(in)            :: nOa
  integer,intent(in)            :: nOb

  double precision,intent(in)   :: S(nBas_AOs,nBas_AOs)
  double precision,intent(in)   :: c(nBas_AOs,nBas_MOs)
  double precision,intent(inout):: FaAO(nBas_AOs,nBas_AOs)
  double precision,intent(inout):: FbAO(nBas_AOs,nBas_AOs)

! Local variables

  double precision              :: aC,bC
  double precision              :: aO,bO
  double precision              :: aV,bV

  integer                       :: nC
  integer                       :: nO
  integer                       :: nV

  double precision,allocatable  :: F(:,:)
  double precision,allocatable  :: Fa(:,:)
  double precision,allocatable  :: Fb(:,:)

! Output variables

  double precision,intent(out)  :: FAO(nBas_AOs,nBas_AOs)

! Memory allocation

  allocate(F(nBas_MOs,nBas_MOs), Fa(nBas_MOs,nBas_MOs), Fb(nBas_MOs,nBas_MOs))

! Roothan canonicalization parameters

  aC = +1.5d0
  bC = -0.5d0

  aO = +0.5d0
  bO = +0.5d0

  aV = -0.5d0
  bV = +1.5d0

! Number of closed, open, and virtual orbitals

  nC = min(nOa, nOb)
  nO = abs(nOa - nOb)
  nV = nBas_AOs - nC - nO

! Block-by-block Fock matrix 

  call AOtoMO(nBas_AOs, nBas_MOs, c, FaAO, Fa)
  call AOtoMO(nBas_AOs, nBas_MOs, c, FbAO, Fb)

  F(1:nC,      1:nC      ) =    aC*Fa(1:nC,      1:nC      ) +    bC*Fb(1:nC,      1:nC      )
  F(1:nC,   nC+1:nC+nO   ) =                                         Fb(1:nC,   nC+1:nC+nO   )
  F(1:nC,nO+nC+1:nC+nO+nV) = 0.5d0*Fa(1:nC,nO+nC+1:nC+nO+nV) + 0.5d0*Fb(1:nC,nO+nC+1:nC+nO+nV)

  F(nC+1:nO+nC,      1:nC      ) =                                         Fb(nC+1:nO+nC,   1:nC   )
  F(nC+1:nO+nC,   nC+1:nC+nO   ) = aO*Fa(nC+1:nO+nC,   nC+1:nC+nO   ) + bO*Fb(nC+1:nO+nC,nC+1:nC+nO)
  F(nC+1:nO+nC,nO+nC+1:nC+nO+nV) =    Fa(nC+1:nO+nC,nO+nC+1:nC+nO+nV)

  F(nO+nC+1:nC+nO+nV,      1:nC      ) = 0.5d0*Fa(nO+nC+1:nC+nO+nV,      1:nC      ) + 0.5d0*Fb(nO+nC+1:nC+nO+nV,      1:nC      )
  F(nO+nC+1:nC+nO+nV,   nC+1:nC+nO   ) =       Fa(nO+nC+1:nC+nO+nV,   nC+1:nC+nO   )
  F(nO+nC+1:nC+nO+nV,nO+nC+1:nC+nO+nV) =    aV*Fa(nO+nC+1:nC+nO+nV,nO+nC+1:nC+nO+nV) +    bV*Fb(nO+nC+1:nC+nO+nV,nO+nC+1:nC+nO+nV)

  call MOtoAO(nBas_AOs, nBas_MOs, S, c, F, FAO) 
  call MOtoAO(nBas_AOs, nBas_MOs, S, c, Fa, FaAO)
  call MOtoAO(nBas_AOs, nBas_MOs, S, c, Fb, FbAO)

  deallocate(F, Fa, Fb)

end subroutine
