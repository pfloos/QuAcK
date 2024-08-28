subroutine AOtoMO(nBas_AOs, nBas_MOs, C, M_AOs, M_MOs)

! Perform AO to MO transformation of a matrix M_AOs for given coefficients c
! M_MOs = C.T M_AOs C

  implicit none

! Input variables

  integer,intent(in)            :: nBas_AOs, nBas_MOs
  double precision,intent(in)   :: C(nBas_AOs,nBas_MOs)
  double precision,intent(in)   :: M_AOs(nBas_AOs,nBas_AOs)

! Local variables

  double precision,allocatable  :: AC(:,:)

! Output variables

  double precision,intent(out)  :: M_MOs(nBas_MOs,nBas_MOs)

  allocate(AC(nBas_AOs,nBas_MOs))

  !AC = matmul(M_AOs, C)
  !M_MOs = matmul(transpose(C), AC)

  call dgemm("N", "N", nBas_AOs, nBas_MOs, nBas_AOs, 1.d0, &
             M_AOs(1,1), nBas_AOs, C(1,1), nBas_AOs,       &
             0.d0, AC(1,1), nBas_AOs)

  call dgemm("T", "N", nBas_MOs, nBas_MOs, nBas_AOs, 1.d0, &
             C(1,1), nBas_AOs, AC(1,1), nBas_AOs,          &
             0.d0, M_MOs(1,1), nBas_MOs)

  deallocate(AC)

end subroutine 
