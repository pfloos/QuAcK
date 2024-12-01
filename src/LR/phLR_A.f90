subroutine phLR_A(ispin,dRPA,nBas,nC,nO,nV,nR,nS,lambda,e,ERI,Aph)

! Compute resonant block of the ph channel

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dRPA
  integer,intent(in)            :: ispin
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: e(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas) 
  
! Local variables

  double precision              :: delta_dRPA
  double precision,external     :: Kronecker_delta

  integer                       :: i,j,a,b,ia,jb
  integer                       :: nn,jb0
  logical                       :: i_eq_j
  double precision              :: ct1,ct2

! Output variables

  double precision,intent(out)  :: Aph(nS,nS)

! Direct RPA

  delta_dRPA = 0d0
  if(dRPA) delta_dRPA = 1d0

! Build A matrix for single manifold

  if(ispin == 1) then 

    nn = nBas - nR - nO
    ct1 = 2d0 * lambda
    ct2 = - (1d0 - delta_dRPA) * lambda
    !$OMP PARALLEL DEFAULT(NONE)                    &
    !$OMP PRIVATE (i, a, j, b, i_eq_j, ia, jb0, jb) &
    !$OMP SHARED (nC, nO, nR, nBas, nn, ct1, ct2, e, ERI, Aph)
    !$OMP DO COLLAPSE(2)
    do i = nC+1, nO
      do a = nO+1, nBas-nR
        ia = a - nO + (i - nC - 1) * nn

        do j = nC+1, nO
          i_eq_j = i == j
          jb0 = (j - nC - 1) * nn - nO

          do b = nO+1, nBas-nR
            jb = b + jb0

            Aph(ia,jb) = ct1 * ERI(i,b,a,j) + ct2 * ERI(i,b,j,a)
            if(i_eq_j) then
              if(a == b) Aph(ia,jb) = Aph(ia,jb) + e(a) - e(i)
            endif
          enddo
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

    !ia = 0
    !do i=nC+1,nO
    !  do a=nO+1,nBas-nR
    !    ia = ia + 1
    !    jb = 0
    !    do j=nC+1,nO
    !      do b=nO+1,nBas-nR
    !        jb = jb + 1
    !        Aph(ia,jb) = (e(a) - e(i))*Kronecker_delta(i,j)*Kronecker_delta(a,b) &
    !                   + 2d0*lambda*ERI(i,b,a,j) - (1d0 - delta_dRPA)*lambda*ERI(i,b,j,a)
    !      end do
    !    end do
    !  end do
    !end do

  end if

! Build A matrix for triplet manifold

  if(ispin == 2) then 

    nn = nBas - nR - nO
    ct2 = - (1d0 - delta_dRPA) * lambda
    !$OMP PARALLEL DEFAULT(NONE)                    &
    !$OMP PRIVATE (i, a, j, b, i_eq_j, ia, jb0, jb) &
    !$OMP SHARED (nC, nO, nR, nBas, nn, ct2, e, ERI, Aph)
    !$OMP DO COLLAPSE(2)
    do i = nC+1, nO
      do a = nO+1, nBas-nR
        ia = a - nO + (i - nC - 1) * nn

        do j = nC+1, nO
          i_eq_j = i == j
          jb0 = (j - nC - 1) * nn - nO

          do b = nO+1, nBas-nR
            jb = b + jb0

            Aph(ia,jb) = ct2 * ERI(i,b,j,a)
            if(i_eq_j) then
              if(a == b) Aph(ia,jb) = Aph(ia,jb) + e(a) - e(i)
            endif
          enddo
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

!    ia = 0
!    do i=nC+1,nO
!      do a=nO+1,nBas-nR
!        ia = ia + 1
!        jb = 0
!        do j=nC+1,nO
!          do b=nO+1,nBas-nR
!            jb = jb + 1
!            Aph(ia,jb) = (e(a) - e(i))*Kronecker_delta(i,j)*Kronecker_delta(a,b) &
!                       - (1d0 - delta_dRPA)*lambda*ERI(i,b,j,a)
!          end do
!        end do
!      end do
!    end do

  end if

end subroutine 
