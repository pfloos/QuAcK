subroutine MOtoAO(nBas, nOrb, S, C, M_MOs, M_AOs)

  ! Perform MO to AO transformation of a matrix M_AOs for a given metric S
  ! and coefficients c
  ! 
  ! M_AOs = S C M_MOs (S C).T

  implicit none

  integer,          intent(in)  :: nBas, nOrb
  double precision, intent(in)  :: S(nBas,nBas)
  double precision, intent(in)  :: C(nBas,nOrb)
  double precision, intent(in)  :: M_MOs(nOrb,nOrb)
  double precision, intent(out) :: M_AOs(nBas,nBas)

  double precision, allocatable :: SC(:,:),BSC(:,:)


  allocate(SC(nBas,nOrb), BSC(nOrb,nBas))

  !SC  = matmul(S, C)
  !BSC = matmul(M_MOs, transpose(SC))
  !M_AOs = matmul(SC, BSC)

  call dgemm("N", "N", nBas, nOrb, nBas, 1.d0, &
             S(1,1), nBas, C(1,1), nBas,       &
             0.d0, SC(1,1), nBas)

  call dgemm("N", "T", nOrb, nBas, nOrb, 1.d0, &
             M_MOs(1,1), nOrb, SC(1,1), nBas,  &
             0.d0, BSC(1,1), nOrb)

  call dgemm("N", "N", nBas, nBas, nOrb, 1.d0, &
             SC(1,1), nBas, BSC(1,1), nOrb,    &
             0.d0, M_AOs(1,1), nBas)

  deallocate(SC, BSC)

end subroutine 
