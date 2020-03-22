subroutine pCCD(maxSCF,thresh,max_diis,nBas,nO,nV,ERI,ENuc,ERHF,eHF)

! pair CCD module

  implicit none

! Input variables

  integer,intent(in)            :: maxSCF
  integer,intent(in)            :: max_diis
  double precision,intent(in)   :: thresh

  integer,intent(in)            :: nBas,nO,nV
  double precision,intent(in)   :: ENuc,ERHF
  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: i,j,a,b

  integer                       :: nSCF
  double precision              :: Conv
  double precision              :: ECCD,EcCCD

  double precision,allocatable  :: delta_OOVV(:,:)

  double precision,allocatable  :: OOOO(:,:)
  double precision,allocatable  :: OOVV(:,:)
  double precision,allocatable  :: OVOV(:,:)
  double precision,allocatable  :: OVVO(:,:)
  double precision,allocatable  :: VVVV(:,:)

  double precision,allocatable  :: X(:,:)
  double precision,allocatable  :: Y(:,:)

  double precision,allocatable  :: r(:,:)
  double precision,allocatable  :: t(:,:)

  double precision,external     :: trace_matrix
           
! Hello world

  write(*,*)
  write(*,*)'**************************************'
  write(*,*)'|     pair CCD calculation           |'
  write(*,*)'**************************************'
  write(*,*)

! Form energy denominator

  allocate(delta_OOVV(nO,nV))

  do i=1,nO
    do a=1,nV
      delta_OOVV(i,a) = 2d0*(eHF(nO+a) - eHF(i))
    enddo
  enddo

! Create integral batches

  allocate(OOOO(nO,nO),OOVV(nO,nV),OVOV(nO,nV),OVVO(nO,nV),VVVV(nV,nV))

  do i=1,nO
    do j=1,nO
      OOOO(i,j) = ERI(i,i,j,j)
    end do
  end do

  do i=1,nO
    do a=1,nV
      OOVV(i,a) = ERI(i,i,nO+a,nO+a)
      OVOV(i,a) = ERI(i,nO+a,i,nO+a)
      OVVO(i,a) = ERI(i,nO+a,nO+a,i)
    end do
  end do

  do a=1,nV
    do b=1,nV
      VVVV(a,b) = ERI(nO+a,nO+a,nO+b,nO+b)
    end do
  end do

! MP2 guess amplitudes

  allocate(t(nO,nV))

  t(:,:) = - OOVV(:,:)/delta_OOVV(:,:)

! Initialization

  allocate(r(nO,nV),X(nV,nV),Y(nO,nO))

  Conv = 1d0
  nSCF = 0

!------------------------------------------------------------------------
! Main SCF loop
!------------------------------------------------------------------------
  write(*,*)
  write(*,*)'----------------------------------------------------'
  write(*,*)'| pair CCD calculation                             |'
  write(*,*)'----------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A16,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X)') &
            '|','#','|','E(pCCD)','|','Ec(pCCD)','|','Conv','|'
  write(*,*)'----------------------------------------------------'

  do while(Conv > thresh .and. nSCF < maxSCF)

  ! Increment 

    nSCF = nSCF + 1

  ! Form intermediate array
    
    X(:,:) = matmul(transpose(OOVV(:,:)),t(:,:))
    Y(:,:) = matmul(t(:,:),transpose(OOVV(:,:)))
    
  ! Compute residual

    do i=1,nO
      do a=1,nV
        r(i,a) = - 2d0*(X(a,a) + Y(i,i))*t(i,a)
      end do
    end do

    r(:,:) = r(:,:) + OOVV(:,:) + delta_OOVV(:,:)*t(:,:)               &
           - 2d0*(2d0*OVOV(:,:) - OVVO(:,:) - OOVV(:,:)*t(:,:))*t(:,:) &
           + matmul(t(:,:),transpose(VVVV(:,:)))                       &
           + matmul(transpose(OOOO(:,:)),t(:,:))                       &
           + matmul(Y(:,:),t)

   ! Check convergence 

    Conv = maxval(abs(r(:,:)))
  
   ! Update amplitudes

    t(:,:) = t(:,:) - r(:,:)/delta_OOVV(:,:)

   ! Compute correlation energy

    EcCCD = trace_matrix(nO,matmul(t(:,:),transpose(OOVV(:,:))))

!   Dump results

    ECCD = ERHF + EcCCD

    write(*,'(1X,A1,1X,I3,1X,A1,1X,F16.10,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X)') &
      '|',nSCF,'|',ECCD+ENuc,'|',EcCCD,'|',Conv,'|'

  enddo
  write(*,*)'----------------------------------------------------'
!------------------------------------------------------------------------
! End of SCF loop
!------------------------------------------------------------------------

! Did it actually converge?

  if(nSCF == maxSCF) then

    write(*,*)
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)'                 Convergence failed                 '
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)

    stop

  endif

end subroutine pCCD
