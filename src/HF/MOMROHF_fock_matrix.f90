
! ---

subroutine MOMROHF_fock_matrix(nBas, nOrb, nOa, nOb, S, c, FaAO, FbAO, FAO,occupations)

! Construct the ROHF Fock matrix in the AO basis for given occupations
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

  integer,intent(in)            :: nBas, nOrb
  integer,intent(in)            :: nOa
  integer,intent(in)            :: nOb
  integer,intent(in)            :: occupations(max(nOa,nOb),nspin)

  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: c(nBas,nOrb)
  double precision,intent(inout):: FaAO(nBas,nBas)
  double precision,intent(inout):: FbAO(nBas,nBas)

! Local variables

  double precision              :: aC,bC
  double precision              :: aO,bO
  double precision              :: aV,bV

  integer                       :: nD,nSA,nSB,nV
  integer,allocatable           :: doublyOccupied(:)
  integer,allocatable           :: virtual(:)
  integer,allocatable           :: singlyOccupiedA(:),singlyOccupiedB(:)

  double precision,allocatable  :: F(:,:)
  double precision,allocatable  :: Fa(:,:)
  double precision,allocatable  :: Fb(:,:)
  
  logical                       :: debug = .false.

! Output variables

  double precision,intent(out)  :: FAO(nBas,nBas)

! Memory allocation

 allocate(F(nOrb,nOrb), Fa(nOrb,nOrb), Fb(nOrb,nOrb),                           &
         doublyOccupied(max(nOa,nOb)),singlyOccupiedA(nOa),singlyOccupiedB(nOb),&
         virtual(nOrb - max(nOa,nOb)))

! Roothan canonicalization parameters

  aC = +1.5d0
  bC = -0.5d0

  aO = +0.5d0
  bO = +0.5d0

  aV = -0.5d0
  bV = +1.5d0

  call MOMROHF_determine_occupations(nOa, nOb,nOrb,occupations,                       &
                                         nSA,nSB,nD,nV,                               &
                                         singlyOccupiedA,singlyOccupiedB,             &
                                         doublyOccupied,virtual)

! Block-by-block Fock matrix 

  call AOtoMO(nBas, nOrb, c, FaAO, Fa)
  call AOtoMO(nBas, nOrb, c, FbAO, Fb)

  
  ! Initialize F to zero
  F = 0.0d0
  
  ! Doubly occupied - Doubly occupied
  F(doublyOccupied(1:nD), doublyOccupied(1:nD))     = aC * Fa(doublyOccupied(1:nD), doublyOccupied(1:nD)) &
                                                    + bC * Fb(doublyOccupied(1:nD), doublyOccupied(1:nD))
  
  ! Doubly occupied - Singly occupied alpha
  F(doublyOccupied(1:nD), singlyOccupiedA(1:nSA))   = Fb(doublyOccupied(1:nD), singlyOccupiedA(1:nSA))
  
  ! Doubly occupied - Singly occupied beta
  F(doublyOccupied(1:nD), singlyOccupiedB(1:nSB))   = Fb(doublyOccupied(1:nD), singlyOccupiedB(1:nSB))
  
  ! Doubly occupied - Virtual
  F(doublyOccupied(1:nD), virtual(1:nV))            = 0.5d0 * Fa(doublyOccupied(1:nD), virtual(1:nV)) &
                                                    + 0.5d0 * Fb(doublyOccupied(1:nD), virtual(1:nV))
  
  ! Singly occupied alpha - Doubly occupied
  F(singlyOccupiedA(1:nSA), doublyOccupied(1:nD))   = Fb(singlyOccupiedA(1:nSA), doublyOccupied(1:nD))
  
  ! Singly occupied alpha - Singly occupied alpha
  F(singlyOccupiedA(1:nSA), singlyOccupiedA(1:nSA)) = aO * Fa(singlyOccupiedA(1:nSA), singlyOccupiedA(1:nSA)) &
                                                    + bO * Fb(singlyOccupiedA(1:nSA), singlyOccupiedA(1:nSA))
  
  ! Singly occupied alpha - Singly occupied beta
  F(singlyOccupiedA(1:nSA), singlyOccupiedB(1:nSB)) = Fa(singlyOccupiedA(1:nSA), singlyOccupiedB(1:nSB))
  
  ! Singly occupied alpha - Virtual
  F(singlyOccupiedA(1:nSA), virtual(1:nV))          = Fa(singlyOccupiedA(1:nSA), virtual(1:nV))
  
  ! Singly occupied beta - Doubly occupied
  F(singlyOccupiedB(1:nSB), doublyOccupied(1:nD))   = Fb(singlyOccupiedB(1:nSB), doublyOccupied(1:nD))
  
  ! Singly occupied beta - Singly occupied alpha
  F(singlyOccupiedB(1:nSB), singlyOccupiedA(1:nSA)) = Fb(singlyOccupiedB(1:nSB), singlyOccupiedA(1:nSA))
  
  ! Singly occupied beta - Singly occupied beta
  F(singlyOccupiedB(1:nSB), singlyOccupiedB(1:nSB)) = aO * Fa(singlyOccupiedB(1:nSB), singlyOccupiedB(1:nSB)) &
                                                    + bO * Fb(singlyOccupiedB(1:nSB), singlyOccupiedB(1:nSB))
  
  ! Singly occupied beta - Virtual
  F(singlyOccupiedB(1:nSB), virtual(1:nV))          = Fa(singlyOccupiedB(1:nSB), virtual(1:nV))
  
  ! Virtual - Doubly occupied
  F(virtual(1:nV), doublyOccupied(1:nD))            = 0.5d0 * Fa(virtual(1:nV), doublyOccupied(1:nD)) &
                                                    + 0.5d0 * Fb(virtual(1:nV), doublyOccupied(1:nD))
  ! Virtual - Singly occupied alpha
  F(virtual(1:nV), singlyOccupiedA(1:nSA))          = Fa(virtual(1:nV), singlyOccupiedA(1:nSA))
  
  ! Virtual - Singly occupied beta
  F(virtual(1:nV), singlyOccupiedB(1:nSB))          = Fa(virtual(1:nV), singlyOccupiedB(1:nSB))
  
  ! Virtual - Virtual
  F(virtual(1:nV), virtual(1:nV))                   = aV * Fa(virtual(1:nV), virtual(1:nV)) &
                                                    + bV * Fb(virtual(1:nV), virtual(1:nV))
  call MOtoAO(nBas, nOrb, S, c, F, FAO) 
  call MOtoAO(nBas, nOrb, S, c, Fa, FaAO)
  call MOtoAO(nBas, nOrb, S, c, Fb, FbAO)

  deallocate(F, Fa, Fb)

  if(debug) then
            
    print *, "-----------------------------------"
    print *, "Double occupied"
    print *, doublyOccupied(1:nD)
    print *, "Only alpha occupied"
    print *, singlyOccupiedA(1:nSA)
    print *, "Only beta occupied"
    print *, singlyOccupiedB(1:nSB)
    print *, "Virtual"
    print *, virtual(1:nV)
    print *, "-----------------------------------"

  end if

end subroutine
