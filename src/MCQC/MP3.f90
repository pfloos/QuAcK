subroutine MP3(nBas,nC,nO,nV,nR,V,EHF,EcMP2,e,EcMP3)

! Perform third-order Moller-Plesset calculation

  implicit none

! Input variables

  integer,intent(in)            :: nBas,nC,nO,nV,nR
  double precision,intent(in)   :: EHF,EcMP2
  double precision,intent(in)   :: V(nBas,nBas,nBas,nBas),e(nBas)

! Local variables

  double precision              :: eps1,eps2,E3a,E3b,E3c

  integer                       :: i,j,k,l,a,b,c,d

! Output variables

  double precision,intent(out)  :: EcMP3

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'|    Moller-Plesset third-order calculation    |'
  write(*,*)'************************************************'
  write(*,*)

! Compute MP3 energy

  E3a = 0d0
  do i=nC+1,nO
  do j=nC+1,nO
  do k=nC+1,nO
  do l=nC+1,nO
    do a=nO+1,nBas-nR
    do b=nO+1,nBas-nR

      eps1 = e(i) + e(j) - e(a) - e(b) 
      eps2 = e(k) + e(l) - e(a) - e(b) 

      E3a = E3a + (V(i,j,a,b) - V(i,j,b,a))* &
                  (V(k,l,i,j) - V(k,l,j,i))* &
                  (V(a,b,k,l) - V(a,b,l,k))/(eps1*eps2)

    enddo
    enddo
  enddo
  enddo
  enddo
  enddo

  E3b = 0d0
  do i=nC+1,nO
  do j=nC+1,nO
    do a=nO+1,nBas-nR
    do b=nO+1,nBas-nR
    do c=nO+1,nBas-nR
    do d=nO+1,nBas-nR

      eps1 = e(i) + e(j) - e(a) - e(b) 
      eps2 = e(i) + e(j) - e(c) - e(d) 

      E3b = E3b + (V(i,j,a,b) - V(i,j,b,a))* &
                  (V(a,b,c,d) - V(a,b,d,c))* &
                  (V(c,d,i,j) - V(c,d,j,i))/(eps1*eps2)

    enddo
    enddo
    enddo
    enddo
  enddo
  enddo

  E3c = 0d0
  do i=nC+1,nO
  do j=nC+1,nO
  do k=nC+1,nO
    do a=nO+1,nBas-nR
    do b=nO+1,nBas-nR
    do c=nO+1,nBas-nR

      eps1 = e(i) + e(j) - e(a) - e(b) 
      eps2 = e(i) + e(k) - e(a) - e(c) 

      E3c = E3c + (V(i,j,a,b) - V(i,j,b,a))* &
                  (V(k,b,c,j) - V(k,b,j,c))* &
                  (V(a,c,i,k) - V(a,c,k,i))/(eps1*eps2)

    enddo
    enddo
    enddo
  enddo
  enddo
  enddo

  EcMP3 = 0.25d0*E3a + 0.25d0*E3b + E3c

  write(*,*)
  write(*,'(A32)')           '-----------------------'
  write(*,'(A32)')           ' MP3 calculation       '
  write(*,'(A32)')           '-----------------------'
  write(*,'(A32,1X,F16.10)') ' MP2 contribution      ',EcMP2
  write(*,'(A32,1X,F16.10)') ' MP3 contribution      ',EcMP3
  write(*,'(A32)')           '-----------------------'
  write(*,'(A32,1X,F16.10)') ' MP3 correlation energy',      EcMP2 + EcMP3
  write(*,'(A32,1X,F16.10)') ' MP3 total       energy',EHF + EcMP2 + EcMP3
  write(*,'(A32)')           '-----------------------'
  write(*,*)

end subroutine MP3
