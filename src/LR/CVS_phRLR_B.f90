subroutine CVS_phRLR_B(ispin,dRPA,nBas,nC,nO,nV,nR,nS,nCVS,nFC,occupations,virtuals,lambda,ERI,Bph)

! Compute the coupling block of the ph channel

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)           :: dRPA
  integer,intent(in)           :: ispin,nBas,nC,nO,nV,nR,nS,nCVS,nFC
  integer,intent(in)           :: occupations(nO-nFC),virtuals(nV)
  double precision,intent(in)  :: lambda
  double precision,intent(in)  :: ERI(nBas,nBas,nBas,nBas)
  
! Local variables

  double precision             :: delta_dRPA

  integer                      :: i,j,a,b,ia,jb
  integer                      :: nn,jb0
  double precision             :: ct1,ct2

! Output variables

  double precision,intent(out) :: Bph(nS,nS)

! Direct RPA

  delta_dRPA = 0d0
  if(dRPA) delta_dRPA = 1d0

! Build B matrix for singlet manifold

  if(ispin == 1) then

    ia = 0
    do i=1,nO-nFC
      do a=nCVS+1,nBas-nO
        ia = ia + 1
        jb = 0
        do j=1,nO-nFC
          do b=nCVS+1,nBas-nO
            jb = jb + 1
            Bph(ia,jb) = 2d0*lambda*ERI(occupations(i),occupations(j),virtuals(a),virtuals(b)) &
                    - (1d0 - delta_dRPA)*lambda*ERI(occupations(i),occupations(j),virtuals(b),virtuals(a))
          end do
        end do
      end do
    end do

  end if

! Build B matrix for triplet manifold

  if(ispin == 2) then

    ia = 0
    do i=1,nO-nFC
      do a=nCVS+1,nBas-nO
        ia = ia + 1
        jb = 0
        do j=1,nO-nFC
          do b=nCVS+1,nBas-nO
            jb = jb + 1
            Bph(ia,jb) = - (1d0 - delta_dRPA)*lambda*ERI(occupations(i),occupations(j),virtuals(b),virtuals(a))
          end do
        end do
      end do
    end do

  end if

end subroutine 
