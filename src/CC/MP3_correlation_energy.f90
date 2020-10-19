subroutine MP3_correlation_energy(nC,nO,nV,nR,OOVV,t2,v,delta_OOVV,EcMP3)

! Compute the MP3 correlation energy 

  implicit none

! Input variables

  integer,intent(in)            :: nC,nO,nV,nR
  double precision,intent(in)   :: OOVV(nO-nC,nO-nC,nV-nR,nV-nR)
  double precision,intent(in)   :: t2(nO-nC,nO-nC,nV-nR,nV-nR)
  double precision,intent(in)   :: v(nO-nC,nO-nC,nV-nR,nV-nR)
  double precision,intent(in)   :: delta_OOVV(nO-nC,nO-nC,nV-nR,nV-nR)

! Local variables

  integer                       :: i,j,a,b

! Output variables

  double precision,intent(out)  :: EcMP3

  EcMP3 = 0d0
  do i=1,nO-nC
    do j=1,nO-nC
      do a=1,nV-nR
        do b=1,nV-nR

          EcMP3 = EcMP3 + OOVV(i,j,a,b)*(t2(i,j,a,b) + v(i,j,a,b)/delta_OOVV(i,j,a,b))

        enddo
      enddo
    enddo
  enddo

  EcMP3 = 0.25d0*EcMP3

end subroutine MP3_correlation_energy
