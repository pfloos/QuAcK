subroutine AOtoMO(nBas,nOrb,C,M_AOs,M_MOs)

  ! Perform AO to MO transformation of a matrix M_AOs for given coefficients c
  ! M_MOs = C.T M_AOs C

  implicit none

  integer,intent(in)            :: nBas, nOrb
  double precision,intent(in)   :: C(nBas,nOrb)
  double precision,intent(in)   :: M_AOs(nBas,nBas)

  double precision,intent(out)  :: M_MOs(nOrb,nOrb)

  double precision,allocatable  :: AC(:,:)

  allocate(AC(nBas,nOrb))

  !AC = matmul(M_AOs, C)
  !M_MOs = matmul(transpose(C), AC)

  call dgemm("N", "N", nBas, nOrb, nBas, 1.d0, &
             M_AOs(1,1), nBas, C(1,1), nBas,   &
             0.d0, AC(1,1), nBas)

  call dgemm("T", "N", nOrb, nOrb, nBas, 1.d0, &
             C(1,1), nBas, AC(1,1), nBas,      &
             0.d0, M_MOs(1,1), nOrb)

  deallocate(AC)

end subroutine 
