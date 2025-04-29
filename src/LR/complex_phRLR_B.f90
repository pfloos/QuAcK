subroutine complex_phRLR_B(ispin,dRPA,nBas,nC,nO,nV,nR,nS,lambda,ERI,Bph)

! Compute the coupling block of the ph channel

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)           :: dRPA
  integer,intent(in)           :: ispin,nBas,nC,nO,nV,nR,nS
  double precision,intent(in)  :: lambda
  complex*16,intent(in)        :: ERI(nBas,nBas,nBas,nBas)
  
! Local variables

  double precision             :: delta_dRPA

  integer                      :: i,j,a,b,ia,jb
  integer                      :: nn,jb0
  complex*16                   :: ct1,ct2

! Output variables

  complex*16,intent(out)       :: Bph(nS,nS)

! Direct RPA

  delta_dRPA = 0d0
  if(dRPA) delta_dRPA = 1d0

! Build B matrix for singlet manifold

  if(ispin == 1) then

    nn = nBas - nR - nO
    ct1 = 2d0 * lambda
    ct2 = - (1d0 - delta_dRPA) * lambda
    !$OMP PARALLEL DEFAULT(NONE)            &
    !$OMP PRIVATE (i, a, j, b, ia, jb0, jb) &
    !$OMP SHARED (nC, nO, nR, nBas, nn, ct1, ct2, ERI, Bph)
    !$OMP DO COLLAPSE(2)
    do i = nC+1, nO
      do a = nO+1, nBas-nR
        ia = a - nO + (i - nC - 1) * nn

        do j = nC+1, nO
          jb0 = (j - nC - 1) * nn - nO

          do b = nO+1, nBas-nR
            jb = b + jb0

            Bph(ia,jb) = ct1 * ERI(b,i,j,a) + ct2 * ERI(b,j,i,a)
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
    !        Bph(ia,jb) = 2d0*lambda*ERI(i,j,a,b) - (1d0 - delta_dRPA)*lambda*ERI(i,j,b,a)
    !      end do
    !    end do
    !  end do
    !end do

  end if

! Build B matrix for triplet manifold

  if(ispin == 2) then

    nn = nBas - nR - nO
    ct2 = - (1d0 - delta_dRPA) * lambda
    !$OMP PARALLEL DEFAULT(NONE)            &
    !$OMP PRIVATE (i, a, j, b, ia, jb0, jb) &
    !$OMP SHARED (nC, nO, nR, nBas, nn, ct2, ERI, Bph)
    !$OMP DO COLLAPSE(2)
    do i = nC+1, nO
      do a = nO+1, nBas-nR
        ia = a - nO + (i - nC - 1) * nn

        do j = nC+1, nO
          jb0 = (j - nC - 1) * nn - nO

          do b = nO+1, nBas-nR
            jb = b + jb0

            Bph(ia,jb) = ct2 * ERI(b,j,i,a)
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
!            Bph(ia,jb) = - (1d0 - delta_dRPA)*lambda*ERI(i,j,b,a)
!          end do
!        end do
!      end do
!    end do

  end if

end subroutine 
