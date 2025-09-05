subroutine sort_ppRPA(nOO,nVV,nPP,Om,Z,Om1,X1,Y1,Om2,X2,Y2)

! Compute the metric matrix for pp-RPA

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nOO
  integer,intent(in)            :: nVV
  integer,intent(in)            :: nPP
  double precision,intent(in)   :: Om(nPP)
  double precision,intent(in)   :: Z(nPP,nPP)
  
! Local variables

  integer                       :: pq,ab,ij
  double precision,allocatable  :: M(:,:)
  double precision,allocatable  :: Z1(:,:)
  double precision,allocatable  :: Z2(:,:)
  double precision,allocatable  :: S1(:,:)
  double precision,allocatable  :: S2(:,:)
  double precision,allocatable  :: O1(:,:)
  double precision,allocatable  :: O2(:,:)  
  double precision,allocatable  :: tmp1(:,:)
  double precision,allocatable  :: tmp2(:,:)

  integer,allocatable           :: order1(:)
  integer,allocatable           :: order2(:)

! Output variables

  double precision,intent(out)  :: Om1(nVV)
  double precision,intent(out)  :: X1(nVV,nVV)
  double precision,intent(out)  :: Y1(nOO,nVV)
  double precision,intent(out)  :: Om2(nOO)
  double precision,intent(out)  :: X2(nVV,nOO)
  double precision,intent(out)  :: Y2(nOO,nOO)


! Memory allocation

  allocate(M(nPP,nPP),Z1(nPP,nVV),Z2(nPP,nOO),order1(nVV),order2(nOO))

! Initializatiom

  Om1(:)  = 0d0
  X1(:,:) = 0d0
  Y1(:,:) = 0d0

  Om2(:)  = 0d0
  X2(:,:) = 0d0
  Y2(:,:) = 0d0

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

  do pq=1,nPP

    if(Om(pq) > 0d0) then 

      ab = ab + 1
      if(ab <= nVV) then
        Om1(ab) = Om(pq)
        Z1(1:nPP,ab) = Z(1:nPP,pq)
      endif

    else

      ij = ij + 1
      if(ij <= nOO) then
        Om2(ij) = Om(pq)
        Z2(1:nPP,ij) = Z(1:nPP,pq)
      endif

    end if

  end do

  if(ab /= nVV) call print_warning('You may have instabilities in pp-RPA [in virt-virt] !!')
  if(ij /= nOO) call print_warning('You may have instabilities in pp-RPA [in occ-occ] !!')
  if(ab /= nVV .or. ij /= nOO) then
    call print_warning('pp-RPA excitation energies (incl. instabilities)')
    do pq=1,nPP
      write(*,'(i5,f10.5)') pq,Om(pq)
    enddo
  endif

  if(nVV > 0) then 

    do ab=1,nVV
      order1(ab) = ab
    end do

    call quick_sort(Om1,order1,nVV)
    call set_order(Z1,order1,nPP,nVV)

  end if

  if(nOO > 0) then

    do ij=1,nOO
      order2(ij) = ij
    end do

    call quick_sort(Om2,order2,nOO)
    call set_order(Z2,order2,nPP,nOO)

  end if
 
  allocate(S1(nVV,nVV),S2(nOO,nOO),O1(nVV,nVV),O2(nOO,nOO))
  allocate(tmp1(nPP,nVV),tmp2(nPP,nOO))

  if(nVV > 0) call dgemm ('N','N',nPP,nVV,nPP,1d0, M,nPP,Z1,nPP,0d0,tmp1,nPP)
  if(nVV > 0) call dgemm ('T','N',nVV,nVV,nPP,1d0,Z1,nPP,tmp1,nPP,0d0,S1,nVV)
  if(nOO > 0) call dgemm ('N','N',nPP,nOO,nPP,1d0, M,nPP,-1d0*Z2,nPP,0d0,tmp2, nPP)
  if(nOO > 0) call dgemm ('T','N',nOO,nOO,nPP,1d0,Z2,nPP,tmp2,nPP,0d0,S2,nOO)

! S1 = + matmul(transpose(Z1),matmul(M,Z1))
! S2 = - matmul(transpose(Z2),matmul(M,Z2))

  if(nVV > 0) call orthogonalize_matrix(1,nVV,S1,O1)
  if(nOO > 0) call orthogonalize_matrix(1,nOO,S2,O2)

  if(nVV > 0) call dgemm ('N','N',nPP,nVV,nVV,1d0,Z1,nPP,O1,nVV,0d0,tmp1,nPP)
  Z1 = tmp1
  if(nOO > 0) call dgemm ('N','N',nPP,nOO,nOO,1d0,Z2,nPP,O2,nOO,0d0,tmp2,nPP)
  Z2 = tmp2

! Z1 = matmul(Z1,O1)
! Z2 = matmul(Z2,O2)

! Define submatrices X1, Y1, X2, & Y2

  X1(1:nVV,1:nVV) = Z1(    1:    nVV,1:nVV)
  Y1(1:nOO,1:nVV) = Z1(nVV+1:nPP,1:nVV)

  X2(1:nVV,1:nOO) = Z2(    1:    nVV,1:nOO)
  Y2(1:nOO,1:nOO) = Z2(nVV+1:nPP,1:nOO)

! call matout(nVV,nVV,X1)
! call matout(nOO,nVV,Y1)

! call matout(nVV,nOO,X2)
! call matout(nOO,nOO,Y2)

! Check orthonormality

! call matout(nVV,nVV,matmul(transpose(X1),X1) - matmul(transpose(Y1),Y1))
! call matout(nOO,nOO,matmul(transpose(X2),X2) - matmul(transpose(Y2),Y2))
! call matout(nVV,nOO,matmul(transpose(X1),X2) - matmul(transpose(Y1),Y2))
! call matout(nOO,nVV,matmul(transpose(X2),X1) - matmul(transpose(Y2),Y1))
 
 deallocate(M,Z1,Z2,order1,order2)

end subroutine 
