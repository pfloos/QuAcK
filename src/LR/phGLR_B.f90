subroutine phGLR_B(dRPA,nOrb,nC,nO,nV,nR,nS,lambda,ERI,Bph)

! Compute the coupling block of the ph channel

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dRPA
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)
  
! Local variables

  double precision              :: delta_dRPA

  integer                       :: i,j,a,b,ia,jb
  integer                       :: nn,jb0
  double precision              :: ct1,ct2
  
! Output variables

  double precision,intent(out)  :: Bph(nS,nS)

! Direct RPA

  delta_dRPA = 0d0
  if(dRPA) delta_dRPA = 1d0

! Build B matrix for spin orbitals 
  ! nn = nOrb - nR - nO
  ! ct1 = lambda
  ! ct2 = - (1d0 - delta_dRPA) * lambda
  
  ! !$OMP PARALLEL DEFAULT(NONE)                    &
  ! !$OMP PRIVATE (i, a, j, b, ia, jb0, jb) &
  ! !$OMP SHARED (nC, nO, nR, nOrb, nn, ct1, ct2, ERI, Bph)
  ! !$OMP DO COLLAPSE(2)
  ! do i = nC+1, nO
  !    do a = nO+1, nOrb-nR
  !       ia = a - nO + (i - nC - 1) * nn

  !       do j = nC+1, nO
  !          jb0 = (j - nC - 1) * nn - nO
  !          do b = nO+1, nOrb-nR
  !             jb = b + jb0

  !             Bph(ia,jb) = ct1 * ERI(i,j,a,b) + ct2 * ERI(i,j,b,a)
              
  !          enddo 
  !       enddo
       
  !    enddo
  ! enddo
  ! !$OMP END DO
  ! !$OMP END PARALLEL
  ia = 0
  do i=nC+1,nO
    do a=nO+1,nOrb-nR
      ia = ia + 1
      jb = 0
      do j=nC+1,nO
        do b=nO+1,nOrb-nR
          jb = jb + 1

          Bph(ia,jb) = lambda*ERI(i,j,a,b) - (1d0 - delta_dRPA)*lambda*ERI(i,j,b,a)

        end do
      end do
    end do
  end do

end subroutine 
