subroutine CCSD_correlation_energy(nC,nO,nV,nR,OOVV,tau,EcCCSD)

! Compute the CCSD energy 

  implicit none

! Input variables

  integer,intent(in)            :: nC,nO,nV,nR
  double precision,intent(in)   :: OOVV(nO,nO,nV,nV)
  double precision,intent(in)   :: tau(nO,nO,nV,nV)

! Local variables

  integer                       :: i,j,a,b

! Output variables

  double precision,intent(out)  :: EcCCSD

  EcCCSD = 0d0
  do i=nC+1,nO
    do j=nC+1,nO
      do a=1,nV-nR
        do b=1,nV-nR

          EcCCSD = EcCCSD + OOVV(i,j,a,b)*tau(i,j,a,b)

        enddo
      enddo
    enddo
  enddo

  EcCCSD = 0.5d0*EcCCSD

end subroutine CCSD_correlation_energy
