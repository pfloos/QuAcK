subroutine MP2(nBas,nC,nO,nV,nR,ERI,ENuc,EHF,e,EcMP2)

! Perform second-order Moller-Plesset calculation

  implicit none

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: EHF
  double precision,intent(in)   :: e(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: i,j,a,b
  double precision              :: eps,E2a,E2b

! Output variables

  double precision,intent(out)  :: EcMP2

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'|    Moller-Plesset second-order calculation   |'
  write(*,*)'************************************************'
  write(*,*)

! Compute MP2 energy

  E2a = 0d0
  E2b = 0d0
  do i=nC+1,nO
  do j=nC+1,nO
    do a=nO+1,nBas-nR
    do b=nO+1,nBas-nR

      eps = e(i) + e(j) - e(a) - e(b) 

!     Second-order ring diagram

      E2a = E2a + ERI(i,j,a,b)*ERI(i,j,a,b)/eps

!     Second-order exchange diagram

      E2b = E2b + ERI(i,j,a,b)*ERI(i,j,b,a)/eps

    enddo
    enddo
  enddo
  enddo

  EcMP2 = 2d0*E2a - E2b

  write(*,*)
  write(*,'(A32)')           '--------------------------'
  write(*,'(A32)')           ' MP2 calculation          '
  write(*,'(A32)')           '--------------------------'
  write(*,'(A32,1X,F16.10)') ' MP2 correlation energy = ',EcMP2
  write(*,'(A32,1X,F16.10)') ' Direct part            = ',2d0*E2a
  write(*,'(A32,1X,F16.10)') ' Exchange part          = ',-E2b
  write(*,'(A32)')           '--------------------------'
  write(*,'(A32,1X,F16.10)') ' MP2 electronic  energy = ',EHF + EcMP2
  write(*,'(A32,1X,F16.10)') ' MP2 total       energy = ',ENuc + EHF + EcMP2
  write(*,'(A32)')           '--------------------------'
  write(*,*)

end subroutine MP2
