subroutine pCCD(maxSCF,thresh,max_diis,nBas,nC,nO,nV,nR,ERI,ENuc,ERHF,eHF)

! pair CCD module

  implicit none

! Input variables

  integer,intent(in)            :: maxSCF
  integer,intent(in)            :: max_diis
  double precision,intent(in)   :: thresh

  integer,intent(in)            :: nBas,nC,nO,nV,nR
  double precision,intent(in)   :: ENuc,ERHF
  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: i,j,a,b

  integer                       :: nSCF
  double precision              :: Conv
  double precision              :: ECCD,EcCCD

  double precision,allocatable  :: eO(:)
  double precision,allocatable  :: eV(:)
  double precision,allocatable  :: delta_OV(:,:)

  double precision,allocatable  :: OOOO(:,:)
  double precision,allocatable  :: OOVV(:,:)
  double precision,allocatable  :: OVOV(:,:)
  double precision,allocatable  :: OVVO(:,:)
  double precision,allocatable  :: VVVV(:,:)

  double precision,allocatable  :: y(:,:)

  double precision,allocatable  :: r(:,:)
  double precision,allocatable  :: t(:,:)

  integer                       :: n_diis
  double precision              :: rcond
  double precision,allocatable  :: error_diis(:,:)
  double precision,allocatable  :: t_diis(:,:)
  double precision,external     :: trace_matrix
          
! Hello world

  write(*,*)
  write(*,*)'**************************************'
  write(*,*)'|     pair CCD calculation           |'
  write(*,*)'**************************************'
  write(*,*)

! Form energy denominator

  allocate(eO(nO-nC),eV(nV-nR),delta_OV(nO-nC,nV-nR))

  eO(:) = eHF(nC+1:nO)
  eV(:) = eHF(nO+1:nBas-nR)

  call form_delta_OV(nC,nO,nV,nR,eO,eV,delta_OV)

! Create integral batches

  allocate(OOOO(nO-nC,nO-nC),OOVV(nO-nC,nV-nR),OVOV(nO-nC,nV-nR),OVVO(nO-nC,nV-nR),VVVV(nV-nR,nV-nR))

  do i=1,nO-nC
    do j=1,nO-nC
      OOOO(i,j) = ERI(nC+i,nC+i,nC+j,nC+j)
    end do
  end do

  do i=1,nO-nC
    do a=1,nV-nR
      OOVV(i,a) = ERI(nC+i,nC+i,nO+a,nO+a)
      OVOV(i,a) = ERI(nC+i,nO+a,nC+i,nO+a)
      OVVO(i,a) = ERI(nC+i,nO+a,nO+a,nC+i)
    end do
  end do

  do a=1,nV-nR
    do b=1,nV-nR
      VVVV(a,b) = ERI(nO+a,nO+a,nO+b,nO+b)
    end do
  end do

! MP2 guess amplitudes

  allocate(t(nO-nC,nV-nR))

  t(:,:) = -0.5d0*OOVV(:,:)/delta_OV(:,:)

  EcCCD = 0d0
  do i=1,nO-nC
    do a=1,nV-nR
      EcCCD = EcCCD + OOVV(i,a)*t(i,a)
    end do
  end do
  print*,'Ec = ',EcCCD

! Memory allocation for DIIS

  allocate(error_diis((nO-nC)*(nV-nR),max_diis),t_diis((nO-nC)*(nV-nR),max_diis))

! Initialization

  allocate(r(nO-nC,nV-nR),y(nO-nC,nO-nC))

  Conv = 1d0
  nSCF = 0

  n_diis          = 0
  t_diis(:,:)     = 0d0
  error_diis(:,:) = 0d0

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
    
   y(:,:) = 0d0
   do i=1,nO-nC
     do j=1,nO-nC
       do b=1,nV-nR
         y(i,j) = y(i,j) + OOVV(j,b)*t(i,b)
       end do
     end do
   end do
    
   ! Compute residual

    do i=1,nO-nC
      do a=1,nV-nR

        r(i,a) = OOVV(i,a) + 2d0*delta_OV(i,a)*t(i,a) & 
               - 2d0*(2d0*OVOV(i,a) - OVVO(i,a) - OOVV(i,a)*t(i,a))*t(i,a)

        do j=1,nO-nC
          r(i,a) = r(i,a) - 2d0*OOVV(j,a)*t(j,a)*t(i,a) + OOOO(j,i)*t(j,a) + y(i,j)*t(j,a) 
        end do 

        do b=1,nV-nR
          r(i,a) = r(i,a) - 2d0*OOVV(i,b)*t(i,b)*t(i,a) + VVVV(a,b)*t(i,b)
        end do 

      end do
    end do

   ! Check convergence 

    Conv = maxval(abs(r(:,:)))
  
   ! Update amplitudes

   t(:,:) = t(:,:) - 0.5d0*r(:,:)/delta_OV(:,:)

   ! Compute correlation energy

    EcCCD = 0d0
    do i=1,nO-nC    
      do a=1,nV-nR
        EcCCD = EcCCD + OOVV(i,a)*t(i,a)
     end do
   end do

   ! Dump results

    ECCD = ERHF + EcCCD

    ! DIIS extrapolation

!   n_diis = min(n_diis+1,max_diis)
!   call DIIS_extrapolation(rcond,nO*nV,nO*nV,n_diis,error_diis,t_diis,-0.5d0*r/delta_OV,t)

    !  Reset DIIS if required

!   if(abs(rcond) < 1d-15) n_diis = 0

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

    stop

  endif

end subroutine pCCD
