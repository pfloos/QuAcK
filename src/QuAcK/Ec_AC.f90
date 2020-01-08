subroutine Ec_AC(ispin,nBas,nC,nO,nV,nR,nS,ERI,XpY,EcAC)

! Compute the correlation energy via the adiabatic connection formula

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: ispin
  integer,intent(in)            :: nBas,nC,nO,nV,nR,nS
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: XpY(nS,nS)

! Local variables

  integer                       :: i,j,a,b
  integer                       :: ia,jb,kc
  double precision              :: delta_spin
  double precision,allocatable  :: P(:,:)

! Output variables

  double precision,intent(out)  :: EcAC

! Singlet or triplet manifold?

  delta_spin = 0d0
  if(ispin == 1) delta_spin = +1d0
  if(ispin == 2) delta_spin = -1d0

! Memory allocation

  allocate(P(nS,nS))

! Compute P = (X+Y)(X+Y) - 1

  P(:,:) = 0d0

  do ia=1,nS
    do jb=1,nS
      do kc=1,nS

        P(ia,jb) = P(ia,jb) + XpY(ia,kc)*XpY(kc,jb)

      enddo
    enddo

    P(ia,ia) = P(ia,ia) - 1d0

  enddo

! Compute Tr[VP]

  EcAC = 0d0

  ia = 0
  do i=nC+1,nO
    do a=nO+1,nBas-nR
      ia = ia + 1
      jb = 0
      do j=nC+1,nO
        do b=nO+1,nBas-nR
          jb = jb + 1

            EcAC = EcAC + (1d0 + delta_spin)*ERI(i,b,a,j)*P(jb,ia)

        enddo
      enddo
    enddo
  enddo

end subroutine Ec_AC

