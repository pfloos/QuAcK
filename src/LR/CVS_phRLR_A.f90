subroutine CVS_phRLR_A(ispin,dRPA,nBas,nC,nO,nV,nR,nS,nCVS,nFC,occupations,virtuals,lambda,e,ERI,Aph)

! Compute resonant block of the ph channel

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)           :: dRPA
  integer,intent(in)           :: ispin
  integer,intent(in)           :: nBas
  integer,intent(in)           :: nC
  integer,intent(in)           :: nO
  integer,intent(in)           :: nV
  integer,intent(in)           :: nR
  integer,intent(in)           :: nS
  integer,intent(in)           :: nCVS
  integer,intent(in)           :: nFC
  integer,intent(in)           :: occupations(nO-nFC)
  integer,intent(in)           :: virtuals(nBas - nO )
  double precision,intent(in)  :: lambda
  double precision,intent(in)  :: e(nBas)
  double precision,intent(in)  :: ERI(nBas,nBas,nBas,nBas) 
  
! Local variables

  double precision             :: delta_dRPA
  double precision,external    :: Kronecker_delta

  integer                      :: i,j,a,b,ia,jb
  integer                      :: nn,jb0
  logical                      :: i_eq_j
  double precision             :: ct1,ct2

! Output variables

  double precision,intent(out) :: Aph(nS,nS)

! Direct RPA

  delta_dRPA = 0d0
  if(dRPA) delta_dRPA = 1d0

! Build A matrix for single manifold

  if(ispin == 1) then 

    ia = 0
    do i=1,nO-nFC
      do a=1+nCVS,nBas-nO
        ia = ia + 1
        jb = 0
        do j=1,nO-nFC
          do b=nCVS+1,nBas-nO
            jb = jb + 1
            Aph(ia,jb) = (e(virtuals(a)) - e(occupations(i)))*Kronecker_delta(i,j)*Kronecker_delta(a,b) &
                       + 2d0*lambda*ERI(occupations(i), &
                       virtuals(b),                     &
                       virtuals(a),                     &
                       occupations(j))                  &
                       - (1d0 - delta_dRPA)*lambda*ERI(occupations(i),virtuals(b),occupations(j),virtuals(a))
          end do
        end do
      end do
    end do

  end if

! Build A matrix for triplet manifold

  if(ispin == 2) then 

    ia = 0
    do i=1,nO - nFC
      do a=nCVS+1,nBas-nO
        ia = ia + 1
        jb = 0
        do j=1,nO -nFC
          do b=nCVS+1,nBas-nO
            jb = jb + 1
            Aph(ia,jb) = (e(virtuals(a)) - e(occupations(i)))*Kronecker_delta(i,j)*Kronecker_delta(a,b) &
                       - (1d0 - delta_dRPA)*lambda*ERI(occupations(i),virtuals(b),occupations(j),virtuals(a))
          end do
        end do
      end do
    end do

  end if

end subroutine 
