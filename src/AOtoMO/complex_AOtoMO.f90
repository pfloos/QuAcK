subroutine complex_AOtoMO(nBas, nOrb, C, M_AOs, M_MOs)

  ! Perform AO to MO transformation of a matrix M_AOs for given coefficients c
  ! M_MOs = C.T M_AOs C

  implicit none

  integer,          intent(in)  :: nBas, nOrb
  complex*16, intent(in)        :: C(nBas,nOrb)
  double precision, intent(in)  :: M_AOs(nBas,nBas)

  complex*16, intent(out)       :: M_MOs(nOrb,nOrb)

  complex*16, allocatable       :: AC(:,:)
  complex*16, allocatable       :: complex_C(:,:)

  allocate(AC(nBas,nOrb))

  AC = matmul(M_AOs, C)
  M_MOs = matmul(transpose(C), AC)
  
  deallocate(AC)

end subroutine 
