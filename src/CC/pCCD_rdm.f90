subroutine pCCD_rdm(O,V,N,ENuc,h,ERI_MO,t2,z2,rdm1,rdm2,ECC)
      
! Compute the 1RDM and 2RDM at the pCCD level

  implicit none
      
! Input variables

  integer,intent(in)            :: O
  integer,intent(in)            :: V
  integer,intent(in)            :: N
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: h(N,N)
  double precision,intent(in)   :: ERI_MO(N,N,N,N)
  double precision,intent(in)   :: t2(O,V)
  double precision,intent(in)   :: z2(O,V)
      
! Local variables

  integer                       :: p,q,r,s
  integer                       :: i,j,a,b

  logical,parameter             :: debug = .false.
  double precision,parameter    :: thresh = 1d-6

  double precision              :: tr_1rdm
  double precision              :: tr_2rdm
  double precision              :: E1
  double precision              :: E2

  double precision,allocatable  :: xOO(:,:)
  double precision,allocatable  :: xVV(:,:)
  double precision,allocatable  :: xOV(:,:)

  double precision,external     :: trace_matrix

! Output variables

  double precision,intent(out)  :: rdm1(N,N)
  double precision,intent(out)  :: rdm2(N,N,N,N)
  double precision,intent(out)  :: ECC

! Allocate memory

  allocate(xOO(O,O),xVV(V,V),xOV(O,V))

! Build intermediates

  xOO(:,:) = matmul(t2,transpose(z2))
  xVV(:,:) = matmul(transpose(z2),t2)
  xOV(:,:) = matmul(t2,matmul(transpose(z2),t2))

! Form 1RDM

  rdm1(:,:) = 0d0

  do i=1,O
    rdm1(i,i) = 2d0*(1d0 - xOO(i,i))
  end do

  do a=1,V
    rdm1(O+a,O+a) = 2d0*xVV(a,a)
  end do

! Check 1RDM

  tr_1rdm = trace_matrix(N,rdm1)
  write(*,'(A25,F16.10)') ' --> Trace of the 1RDM = ',tr_1rdm

  if( abs(dble(2*O) - tr_1rdm) > thresh ) & 
  write(*,*) ' !!! Your 1RDM seems broken !!! '
  write(*,*)

  if(debug) then

    write(*,*) '1RDM is diagonal at the pCCD level:'
    call matout(N,N,rdm1)

  end if

! Form 2RM

  rdm2(:,:,:,:) = 0d0

  ! iijj

  do i=1,O
    do j=1,O
      rdm2(i,i,j,j) = 2d0*xOO(i,j)
    end do
  end do

  ! iiaa

  do i=1,O
    do a=1,V
      rdm2(i,i,O+a,O+a) = 2d0*(t2(i,a) + xOV(i,a) - 2d0*t2(i,a)*(xVV(a,a) + xOO(i,i) - t2(i,a)*z2(i,a)))
    end do
  end do

  ! aaii

  do i=1,O
    do a=1,V
      rdm2(O+a,O+a,i,i) = 2d0*z2(i,a)
    end do
  end do

  ! aabb

  do a=1,V
    do b=1,V
      rdm2(O+a,O+a,O+b,O+b) = 2d0*xVV(a,b)
    end do
  end do

  ! ijij

  do i=1,O
    do j=1,O
      rdm2(i,j,i,j) = 4d0*(1d0 - xOO(i,i) - xOO(j,j))
    end do
  end do

  ! ijji

  do i=1,O
    do j=1,O
      rdm2(i,j,j,i) = - 2d0*(1d0 - xOO(i,i) - xOO(j,j))
    end do
  end do

  ! iiii

  do i=1,O
    rdm2(i,i,i,i) = 2d0*(1d0 - xOO(i,i))
  end do

  ! iaia

  do i=1,O
    do a=1,V
      rdm2(i,O+a,i,O+a) = 4d0*(xVV(a,a) - t2(i,a)*z2(i,a))
    end do
  end do

  ! iaai

  do i=1,O
    do a=1,V
      rdm2(i,O+a,O+a,i) = - 2d0*(xVV(a,a) - t2(i,a)*z2(i,a))
    end do
  end do

  ! aiai

  do i=1,O
    do a=1,V
      rdm2(O+a,i,O+a,i) = 4d0*(xVV(a,a) - t2(i,a)*z2(i,a))
    end do
  end do

  ! aiia

  do i=1,O
    do a=1,V
      rdm2(O+a,i,i,O+a) = - 2d0*(xVV(a,a) - t2(i,a)*z2(i,a))
    end do
  end do

  ! abab

  do a=1,V
    rdm2(O+a,O+a,O+a,O+a) = 2d0*xVV(a,a)
  end do

! Check 2RDM

  tr_2rdm = trace_matrix(N**2,rdm2)
  write(*,'(A25,F16.10)') ' --> Trace of the 2RDM = ',tr_2rdm

  if( abs(dble(2*O*(2*O-1)) - tr_2rdm) > thresh ) & 
  write(*,*) ' !!! Your 2RDM seems broken !!! '
  write(*,*)

  if(debug) then

    write(*,*) '2RDM is not diagonal at the pCCD level:'
    call matout(N**2,N**2,rdm2)

  endif

  deallocate(xOO,xVV,xOV)

! Compute electronic energy

  E1 = 0d0
  E2 = 0d0

  do p=1,N
    do q=1,N
      E1 = E1 + rdm1(p,q)*h(p,q)
      do r=1,N
        do s=1,N
          E2 = E2 + rdm2(p,q,r,s)*ERI_MO(p,q,r,s)
        end do
      end do
    end do
  end do

  E2 = 0.5d0*E2

  ECC = E1 + E2

  write(*,'(A25,F16.10)') ' One-electron energy = ',E1
  write(*,'(A25,F16.10)') ' Two-electron energy = ',E2
  write(*,'(A25,F16.10)') ' Electronic   energy = ',ECC
  write(*,'(A25,F16.10)') ' Total pCCD   energy = ',ECC + ENuc
  write(*,*)

end
