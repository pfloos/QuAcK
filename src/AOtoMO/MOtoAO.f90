subroutine MOtoAO(nBas_AOs, nBas_MOs, S, C, M_MOs, M_AOs)

! Perform MO to AO transformation of a matrix M_AOs for a given metric S
! and coefficients c
! 
! M_AOs = S C M_MOs (S C).T
!

  implicit none

! Input variables

  integer,intent(in)            :: nBas_AOs, nBas_MOs
  double precision,intent(in)   :: S(nBas_AOs,nBas_AOs)
  double precision,intent(in)   :: C(nBas_AOs,nBas_MOs)
  double precision,intent(in)   :: M_MOs(nBas_MOs,nBas_MOs)

! Local variables

  double precision,allocatable  :: SC(:,:),BSC(:,:)

! Output variables

  double precision,intent(out)  :: M_AOs(nBas_AOs,nBas_AOs)

! Memory allocation

  allocate(SC(nBas_AOs,nBas_MOs), BSC(nBas_MOs,nBas_AOs))

  !SC  = matmul(S, C)
  !BSC = matmul(M_MOs, transpose(SC))
  !M_AOs = matmul(SC, BSC)

  call dgemm("N", "N", nBas_AOs, nBas_MOs, nBas_AOs, 0.d0, &
             S(1,1), nBas_AOs, C(1,1), nBas_AOs,           &
             1.d0, SC(1,1), nBas_AOs)

  call dgemm("N", "T", nBas_MOs, nBas_AOs, nBas_MOs, 0.d0, &
             M_MOs(1,1), nBas_MOs, SC(1,1), nBas_AOs,      &
             1.d0, BSC(1,1), nBas_MOs)

  call dgemm("N", "N", nBas_AOs, nBas_AOs, nBas_MOs, 0.d0, &
             SC(1,1), nBas_AOs, BSC(1,1), nBas_MOs,        &
             1.d0, M_AOs(1,1), nBas_AOs)

  deallocate(SC, BSC)

end subroutine 
