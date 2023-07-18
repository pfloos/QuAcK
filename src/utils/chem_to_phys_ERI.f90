subroutine chem_to_phys_ERI(nBas,ERI)

! Antisymmetrize ERIs

  implicit none

! Input variables

  integer,intent(in)            :: nBas
  double precision,intent(inout):: ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: p,q,r,s
  double precision,allocatable  :: pERI(:,:,:,:)

  allocate(pERI(nBas,nBas,nBas,nBas))

  do p=1,nBas
    do q=1,nBas
      do r=1,nBas
        do s=1,nBas

          pERI(p,q,r,s) = ERI(p,r,q,s)

        enddo
      enddo
    enddo
  enddo

  ERI(:,:,:,:) = pERI(:,:,:,:)

end subroutine 
