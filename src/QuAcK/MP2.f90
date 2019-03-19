subroutine MP2(nBas,nC,nO,nV,nR,ERI,ENuc,EHF,e,EcMP2)

! Perform third-order Moller-Plesset calculation

  implicit none

! Input variables

  integer,intent(in)            :: nBas,nC,nO,nV,nR
  double precision,intent(in)   :: ENuc,EHF
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas),e(nBas)

! Local variables

  integer                       :: i,j,a,b
  double precision              :: eps,E2a,E2b

! Output variables

  double precision,intent(out)  :: EcMP2(3)

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

  EcMP2(2) = 2d0*E2a
  EcMP2(3) = -E2b
  EcMP2(1) = EcMP2(2) + EcMP2(3)

  write(*,*)
  write(*,'(A32)')           '-----------------------'
  write(*,'(A32)')           ' MP2 calculation       '
  write(*,'(A32)')           '-----------------------'
  write(*,'(A32,1X,F16.10)') ' MP2 correlation energy',EcMP2(1)
  write(*,'(A32,1X,F16.10)') ' Direct part           ',EcMP2(2)
  write(*,'(A32,1X,F16.10)') ' Exchange part         ',EcMP2(3)
  write(*,'(A32)')           '-----------------------'
  write(*,'(A32,1X,F16.10)') ' MP2 electronic  energy',EHF + EcMP2(1)
  write(*,'(A32,1X,F16.10)') ' MP2 total       energy',ENuc + EHF + EcMP2(1)
  write(*,'(A32)')           '-----------------------'
  write(*,*)

end subroutine MP2
