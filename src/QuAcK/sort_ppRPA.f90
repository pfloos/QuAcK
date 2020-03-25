subroutine sort_ppRPA(ortho_eigvec,nOO,nVV,Omega,Z,Omega1,X1,Y1,Omega2,X2,Y2)

! Compute the metric matrix for pp-RPA

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: ortho_eigvec
  integer,intent(in)            :: nOO
  integer,intent(in)            :: nVV
  double precision,intent(in)   :: Omega(nOO+nVV)
  double precision,intent(in)   :: Z(nOO+nVV,nOO+nVV)
  
! Local variables

  integer                       :: pq,ab,ij
  integer                       :: deg1,ab_start,ab_end
  integer                       :: deg2,ij_start,ij_end
  double precision,allocatable  :: M(:,:)
  double precision,allocatable  :: Z1(:,:)
  double precision,allocatable  :: Z2(:,:)
  double precision,allocatable  :: S1(:,:)
  double precision,allocatable  :: S2(:,:)
  double precision,allocatable  :: O1(:,:)
  double precision,allocatable  :: O2(:,:)

  integer,allocatable           :: order1(:)
  integer,allocatable           :: order2(:)

! Output variables

  double precision,intent(out)  :: Omega1(nVV)
  double precision,intent(out)  :: X1(nVV,nVV)
  double precision,intent(out)  :: Y1(nOO,nVV)
  double precision,intent(out)  :: Omega2(nOO)
  double precision,intent(out)  :: X2(nVV,nOO)
  double precision,intent(out)  :: Y2(nOO,nOO)

! Memory allocation

 allocate(M(nOO+nVV,nOO+nVV),              &
          Z1(nOO+nVV,nVV),Z2(nOO+nVV,nOO), &
          order1(nVV),order2(nOO))

! Initializatiom

  Omega1(:) = 0d0
  X1(:,:)   = 0d0
  Y1(:,:)   = 0d0

  Omega2(:) = 0d0
  X2(:,:)   = 0d0
  Y2(:,:)   = 0d0

! Compute metric 

  M(:,:) = 0d0
  
  do ab=1,nVV
    M(ab,ab) = 1d0
  end do

  do ij=1,nOO
    M(nVV+ij,nVV+ij) = -1d0
  end do

! Start sorting eigenvectors

  ab = 0
  ij = 0

  do pq=1,nOO+nVV

    if(Omega(pq) > 0d0) then 

      ab = ab + 1
      Omega1(ab) = Omega(pq)
      Z1(1:nOO+nVV,ab) = Z(1:nOO+nVV,pq)

    else

      ij = ij + 1
      Omega2(ij) = Omega(pq)
      Z2(1:nOO+nVV,ij) = Z(1:nOO+nVV,pq)

    end if

  end do

  if(minval(Omega1(:)) < 0d0 .or. ab /= nVV) call print_warning('You may have instabilities in pp-RPA!!')
  if(maxval(Omega2(:)) > 0d0 .or. ij /= nOO) call print_warning('You may have instabilities in pp-RPA!!')

  do ab=1,nVV
    order1(ab) = ab
  end do

  do ij=1,nOO
    order2(ij) = ij
  end do

  call quick_sort(Omega1(:),order1(:),nVV)
  call set_order(Z1(:,:),order1(:),nOO+nVV,nVV)

  call quick_sort(Omega2(:),order2(:),nOO)
  call set_order(Z2(:,:),order2(:),nOO+nVV,nOO)

! write(*,*) 'pp-RPA positive excitation energies'
! call matout(nVV,1,Omega1(:))
! write(*,*)

! write(*,*) 'pp-RPA negative excitation energies'
! call matout(nOO,1,Omega2(:))
! write(*,*)

! Orthogonalize eigenvectors

  if(ortho_eigvec) then

    ! Find degenerate eigenvalues

      deg1 = 1
      ab_start = 1

      do ab=1,nVV

        if(ab < nVV .and. abs(Omega1(ab) - Omega1(ab+1)) < 1d-10) then 

          if(deg1 == 1) ab_start = ab
          deg1 = deg1 + 1

        else

          ab_end = ab 
!         print*,'deg = ',deg1,ab_start,ab_end

          allocate(S1(deg1,deg1),O1(deg1,deg1))

          S1 = matmul(transpose(Z1(:,ab_start:ab_end)),matmul(M,Z1(:,ab_start:ab_end)))
          call orthogonalization_matrix(1,deg1,S1,O1)
          Z1(:,ab_start:ab_end) = matmul(Z1(:,ab_start:ab_end),O1)

          deallocate(S1,O1)

          deg1 = 1
          ab_start = ab + 1

        end if
      end do

      deg2 = 1
      ij_start = 1

      do ij=1,nOO

        if(ij < nOO .and. abs(Omega2(ij) - Omega2(ij+1)) < 1d-10) then 

          if(deg2 == 1) ij_start = ij 
          deg2 = deg2 + 1

        else

          ij_end = ij 
!         print*,'deg = ',deg2,ij_start,ij_end

          allocate(S2(deg2,deg2),O2(deg2,deg2))

          S2 = - matmul(transpose(Z2(:,ij_start:ij_end)),matmul(M,Z2(:,ij_start:ij_end)))
          call orthogonalization_matrix(1,deg2,S2,O2)
          Z2(:,ij_start:ij_end) = matmul(Z2(:,ij_start:ij_end),O2)

          deallocate(S2,O2)

          deg2 = 1
          ij_start = ij + 1

        end if
      end do

!   allocate(S1(nVV,nVV),S2(nOO,nOO),O1(nVV,nVV),O2(nOO,nOO))
!   S1 = + matmul(transpose(Z1),matmul(M,Z1))
!   S2 = - matmul(transpose(Z2),matmul(M,Z2))

!   if(nVV > 0) call orthogonalization_matrix(1,nVV,S1,O1)
!   if(nOO > 0) call orthogonalization_matrix(1,nOO,S2,O2)

!   Z1 = matmul(Z1,O1)
!   Z2 = matmul(Z2,O2)

  end if

! Define submatrices X1, Y1, X2, & Y2

  X1(1:nVV,1:nVV) = + Z1(    1:    nVV,1:nVV)
  Y1(1:nOO,1:nVV) = - Z1(nVV+1:nOO+nVV,1:nVV)

  X2(1:nVV,1:nOO) = + Z2(    1:    nVV,1:nOO)
  Y2(1:nOO,1:nOO) = - Z2(nVV+1:nOO+nVV,1:nOO)

! write(*,*) 'Z1t.M.Z1'
! call matout(nVV,nVV,matmul(matmul(transpose(Z1),M),Z1))
! write(*,*) 'Z2t.M.Z2'
! call matout(nOO,nOO,matmul(matmul(transpose(Z2),M),Z2))

! write(*,*) 'X1t.X1 - Y1t.Y1'
! call matout(nVV,nVV,matmul(transpose(X1),X1) - matmul(transpose(Y1),Y1))
! write(*,*) 'X2t.X2 - Y2t.Y2'
! call matout(nOO,nOO,matmul(transpose(X2),X2) - matmul(transpose(Y2),Y2))
! write(*,*) 'X1t.X2 - Y1t.Y2'
! call matout(nVV,nOO,matmul(transpose(X1),X2) - matmul(transpose(Y1),Y2))


end subroutine sort_ppRPA
