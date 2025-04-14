subroutine AOtoMO_GHF(nBas,nOrb,Ca,Cb,A,B)

! Perform AO to MO transformation of a matrix A for given coefficients c

  implicit none

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  double precision,intent(in)   :: Ca(nBas,nOrb)
  double precision,intent(in)   :: Cb(nBas,nOrb)
  double precision,intent(in)   :: A(nBas,nBas)

! Local variables

  double precision,allocatable  :: AC(:,:)
! double precision,allocatable  :: Ba(:,:)

! Output variables

  double precision,intent(out)  :: B(nOrb,nOrb)

  allocate(AC(nBas,nOrb))
! allocate(Ba(nOrb,nOrb))

  AC = matmul(A,Ca)
  B  = matmul(transpose(Ca),AC)

! call dgemm("N", "N", nBas, nOrb, nBas, 1.d0, &
!            A(1,1), nBas, Ca(1,1), nBas,      &
!            0.d0, AC(1,1), nBas)

! call dgemm("T", "N", nOrb, nOrb, nBas, 1.d0, &
!            Ca(1,1), nBas, AC(1,1), nBas,     &
!            0.d0, Ba(1,1), nOrb)

  AC = matmul(A,Cb)
  B  = B + matmul(transpose(Cb),AC)

! call dgemm("N", "N", nBas, nOrb, nBas, 1.d0, &
!            A(1,1), nBas, Cb(1,1), nBas,      &
!            0.d0, AC(1,1), nBas)

! call dgemm("T", "N", nOrb, nOrb, nBas, 1.d0, &
!            Cb(1,1), nBas, AC(1,1), nBas,     &
!            0.d0, B(1,1), nOrb)

! B(:,:) = Ba(:,:) + B(:,:)

end subroutine 
