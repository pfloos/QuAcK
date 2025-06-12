subroutine complex_MOtoAO(nBas, nOrb, S, C, M_MOs, M_AOs)

  ! Perform MO to AO transformation of a matrix M_AOs for a given metric S
  ! and coefficients c
  ! 
  ! M_AOs = S C M_MOs (S C).T

  implicit none

  integer,          intent(in)  :: nBas, nOrb
  double precision, intent(in)  :: S(nBas,nBas)
  complex*16, intent(in)        :: C(nBas,nOrb)
  complex*16, intent(in)        :: M_MOs(nOrb,nOrb)
  complex*16, intent(out)       :: M_AOs(nBas,nBas)

  complex*16, allocatable       :: SC(:,:),BSC(:,:),cS(:,:)


  allocate(SC(nBas,nOrb), BSC(nOrb,nBas),cS(nBas,nBas))
  cS(:,:) = (1d0,0d0)*S(:,:)
  !SC  = matmul(S, C)
  !BSC = matmul(M_MOs, transpose(SC))
  !M_AOs = matmul(SC, BSC)

  call zgemm("N", "N", nBas, nOrb, nBas, 1.d0, &
             cS(1,1), nBas, C(1,1), nBas,       &
             0.d0, SC(1,1), nBas)

  call zgemm("N", "T", nOrb, nBas, nOrb, 1.d0, &
             M_MOs(1,1), nOrb, SC(1,1), nBas,  &
             0.d0, BSC(1,1), nOrb)

  call zgemm("N", "N", nBas, nBas, nOrb, 1.d0, &
             SC(1,1), nBas, BSC(1,1), nOrb,    &
             0.d0, M_AOs(1,1), nBas)

  deallocate(SC, BSC)

end subroutine 
