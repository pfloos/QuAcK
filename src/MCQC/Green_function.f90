subroutine Green_function(nBas,nO,nV,nWalk,nWP,cO,cV,eO_Quad,eV_Quad,AO, & 
                          o1MO,o2MO,v1MO,v2MO,o11,o12,o21,o22,v11,v12,v21,v22)

! Calculate the Green functions

  implicit none

  include 'parameters.h'
  include 'quadrature.h'

! Input variables

  integer,intent(in)            :: nBas,nO,nV,nWalk,nWP
  double precision,intent(in)   :: AO(nWalk,2,nBas),cO(nBas,nO),cV(nBas,nV), &
                                   eO_Quad(nQuad,nO),eV_Quad(nQuad,nV)

! Local variables

  integer                       :: kW,lW,klW,i,a,q
  double precision              :: o1MO(nWalk,nO),o2MO(nWalk,nO),v1MO(nWalk,nV),v2MO(nWalk,nV)

! Output variables

  double precision,intent(out)  :: o11(nQuad,nWP),o12(nQuad,nWP),o21(nQuad,nWP),o22(nQuad,nWP)
  double precision,intent(out)  :: v11(nQuad,nWP),v12(nQuad,nWP),v21(nQuad,nWP),v22(nQuad,nWP)

! Calculate occupied and virtual MOs

  o1MO = matmul(AO(:,1,:),cO)
  o2MO = matmul(AO(:,2,:),cO)
  v1MO = matmul(AO(:,1,:),cV)
  v2MO = matmul(AO(:,2,:),cV)

! Compute occupied Green functions
  o11 = 0d0
  o12 = 0d0
  o21 = 0d0
  o22 = 0d0
  v11 = 0d0
  v12 = 0d0
  v21 = 0d0
  v22 = 0d0

  do q=1,nQuad
    klW = 0
    do kW=1,nWalk-1
      do lW=kW+1,nWalk
        klW = klW + 1
        do i=1,nO
          o11(q,klW) = o11(q,klW) + o1MO(kW,i)*o1MO(lW,i)*eO_Quad(q,i)
          o12(q,klW) = o12(q,klW) + o1MO(kW,i)*o2MO(lW,i)*eO_Quad(q,i)
          o21(q,klW) = o21(q,klW) + o2MO(kW,i)*o1MO(lW,i)*eO_Quad(q,i)
          o22(q,klW) = o22(q,klW) + o2MO(kW,i)*o2MO(lW,i)*eO_Quad(q,i)
        enddo
        do a=1,nV
          v11(q,klW) = v11(q,klW) + v1MO(kW,a)*v1MO(lW,a)*eV_Quad(q,a)
          v12(q,klW) = v12(q,klW) + v1MO(kW,a)*v2MO(lW,a)*eV_Quad(q,a)
          v21(q,klW) = v21(q,klW) + v2MO(kW,a)*v1MO(lW,a)*eV_Quad(q,a)
          v22(q,klW) = v22(q,klW) + v2MO(kW,a)*v2MO(lW,a)*eV_Quad(q,a)
        enddo
      enddo
    enddo
  enddo

end subroutine Green_function
