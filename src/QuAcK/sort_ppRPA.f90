subroutine sort_ppRPA(nOO,nVV,Omega,Z,Omega1,X1,Y1,Omega2,X2,Y2)

! Compute the metric matrix for pp-RPA

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nOO
  integer,intent(in)            :: nVV
  double precision,intent(in)   :: Omega(nOO+nVV)
  double precision,intent(in)   :: Z(nOO+nVV,nOO+nVV)
  
! Local variables

  integer                       :: pq,ab,ij

! Output variables

  double precision,intent(out)  :: Omega1(nVV)
  double precision,intent(out)  :: X1(nVV,nVV)
  double precision,intent(out)  :: Y1(nOO,nVV)
  double precision,intent(out)  :: Omega2(nOO)
  double precision,intent(out)  :: X2(nVV,nOO)
  double precision,intent(out)  :: Y2(nOO,nOO)

! Initializatiom

  Omega1(:) = 0d0
  X1(:,:)   = 0d0
  Y1(:,:)   = 0d0

  Omega2(:) = 0d0
  X2(:,:)   = 0d0
  Y2(:,:)   = 0d0

  ab = 0
  ij = 0

  do pq=1,nOO+nVV

    if(Omega(pq) > 0d0) then 

      ab = ab + 1
      Omega1(ab) = Omega(pq)
      X1(1:nVV,ab) = Z(    1:    nVV,pq)
      Y1(1:nOO,ab) = Z(nVV+1:nOO+nVV,pq)

    else

      ij = ij + 1
      Omega2(ij) = Omega(pq)
      X2(1:nVV,ij) = Z(    1:    nVV,pq)
      Y2(1:nOO,ij) = Z(nVV+1:nOO+nVV,pq)

    end if

  end do

  if(minval(Omega1(:)) < 0d0 .or. ab /= nVV) call print_warning('You may have instabilities in pp-RPA!!')
  if(maxval(Omega2(:)) > 0d0 .or. ij /= nOO) call print_warning('You may have instabilities in pp-RPA!!')

  write(*,*) 'pp-RPA positive excitation energies'
  call matout(nVV,1,Omega1(:))
  write(*,*)

  write(*,*) 'pp-RPA negative excitation energies'
  call matout(nOO,1,Omega2(:))
  write(*,*)

! write(*,*) 'X1.X1 - Y1.Y1'
! call matout(nVV,nVV,matmul(transpose(X1),X1) - matmul(transpose(Y1),Y1))
! write(*,*) 'X2.X2 - Y2.Y2'
! call matout(nOO,nOO,matmul(transpose(X2),X2) - matmul(transpose(Y2),Y2))
! write(*,*) 'X1.X2 - Y1.Y2'
! call matout(nVV,nOO,matmul(transpose(X1),X2) - matmul(transpose(Y1),Y2))


end subroutine sort_ppRPA
