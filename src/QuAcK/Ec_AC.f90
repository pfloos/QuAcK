subroutine Ec_AC(ispin,dRPA,nBas,nC,nO,nV,nR,nS,ERI,XpY,EcAC)

! Compute the correlation energy via the adiabatic connection formula

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: ispin
  logical,intent(in)            :: dRPA
  integer,intent(in)            :: nBas,nC,nO,nV,nR,nS
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: XpY(nS,nS)

! Local variables

  integer                       :: i,j,a,b
  integer                       :: ia,jb,kc
  double precision              :: delta_spin
  double precision              :: delta_dRPA
  double precision,allocatable  :: P(:,:)
  double precision,allocatable  :: V(:,:)
  double precision,external     :: trace_matrix

! Output variables

  double precision,intent(out)  :: EcAC

! Singlet or triplet manifold?

  delta_spin = 0d0
  if(ispin == 1) delta_spin = +1d0
  if(ispin == 2) delta_spin = -1d0

! Direct RPA

  delta_dRPA = 0d0
  if(dRPA) delta_dRPA = 1d0

! Memory allocation

  allocate(P(nS,nS),V(nS,nS))

! Compute P = (X+Y)(X+Y) - 1

  P(:,:) = matmul(transpose(XpY),XpY)

  do ia=1,nS
    P(ia,ia) = P(ia,ia) - 1d0
  enddo

! Compute Viajb = (ia|bj)

  ia = 0
  do i=nC+1,nO
    do a=nO+1,nBas-nR
      ia = ia + 1
      jb = 0
      do j=nC+1,nO
        do b=nO+1,nBas-nR
          jb = jb + 1

            V(ia,jb) = (1d0 + delta_spin)*ERI(i,b,a,j) 

        enddo
      enddo
    enddo
  enddo

! Compute Tr(VP)

  EcAC = trace_matrix(nS,matmul(V,P))

end subroutine Ec_AC

