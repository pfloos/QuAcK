subroutine CCD_correlation_energy(nC,nO,nV,nR,OOVV,t2,EcCCD)

! Compute the CCD energy 

  implicit none

! Input variables

  integer,intent(in)            :: nC,nO,nV,nR
  double precision,intent(in)   :: OOVV(nO,nO,nV,nV)
  double precision,intent(in)   :: t2(nO,nO,nV,nV)

! Local variables

  integer                       :: i,j,a,b

! Output variables

  double precision,intent(out)  :: EcCCD

  EcCCD = 0d0
  do i=nC+1,nO
    do j=nC+1,nO
      do a=1,nV-nR
        do b=1,nV-nR

          EcCCD = EcCCD + OOVV(i,j,a,b)*t2(i,j,a,b)

        enddo
      enddo
    enddo
  enddo

  EcCCD = 0.25d0*EcCCD

end subroutine CCD_correlation_energy
