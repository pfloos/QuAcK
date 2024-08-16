subroutine pCCD(dotest,maxSCF,thresh,max_diis,nBas,nC,nO,nV,nR,ERI,ENuc,ERHF,eHF)

! pair CCD module

  implicit none

! Input variables

  logical,intent(in)            :: dotest

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
  double precision              :: ECC,EcCC

  double precision,allocatable  :: eO(:)
  double precision,allocatable  :: eV(:)
  double precision,allocatable  :: delta_OV(:,:)

  double precision,allocatable  :: OOOO(:,:)
  double precision,allocatable  :: OOVV(:,:)
  double precision,allocatable  :: OVOV(:,:)
  double precision,allocatable  :: OVVO(:,:)
  double precision,allocatable  :: VVVV(:,:)

  double precision,allocatable  :: yO(:,:),yV(:,:)

  double precision,allocatable  :: r(:,:)
  double precision,allocatable  :: t(:,:)
  double precision,allocatable  :: z(:,:)

  double precision,allocatable  :: rdm1(:,:)
  double precision,allocatable  :: rdm2(:,:,:,:)
  double precision,allocatable  :: xOO(:,:)
  double precision,allocatable  :: xVV(:,:)
  double precision,allocatable  :: xOV(:,:)
  double precision              :: tr_1rdm
  double precision              :: tr_2rdm

  integer                       :: O,V
  integer                       :: n_diis
  double precision              :: rcond
  double precision,allocatable  :: err_diis(:,:)
  double precision,allocatable  :: t_diis(:,:)
  double precision,allocatable  :: z_diis(:,:)
  double precision,external     :: trace_matrix
          
! Hello world

  write(*,*)
  write(*,*)'**************************************'
  write(*,*)'|     pair CCD calculation           |'
  write(*,*)'**************************************'
  write(*,*)

! Useful quantities

  O = nO - nC
  V = nV - NR

! Form energy denominator

  allocate(eO(O),eV(V),delta_OV(O,V))

  eO(:) = eHF(nC+1:nO)
  eV(:) = eHF(nO+1:nBas-nR)

  call form_delta_OV(nC,nO,nV,nR,eO,eV,delta_OV)

! Create integral batches

  allocate(OOOO(O,O),OOVV(O,V),OVOV(O,V),OVVO(O,V),VVVV(V,V))

  do i=1,O
    do j=1,O
      OOOO(i,j) = ERI(nC+i,nC+i,nC+j,nC+j)
    end do
  end do

  do i=1,O
    do a=1,V
      OOVV(i,a) = ERI(nC+i,nC+i,nO+a,nO+a)
      OVOV(i,a) = ERI(nC+i,nO+a,nC+i,nO+a)
      OVVO(i,a) = ERI(nC+i,nO+a,nO+a,nC+i)
    end do
  end do

  do a=1,V
    do b=1,V
      VVVV(a,b) = ERI(nO+a,nO+a,nO+b,nO+b)
    end do
  end do

! Initialization

  allocate(t(O,V),r(O,V),yO(O,O),yV(V,V))

! Memory allocation for DIIS

  allocate(err_diis(O*V,max_diis),t_diis(O*V,max_diis))

!------------------------------------------------------------------------
! Compute t ampltiudes
!------------------------------------------------------------------------

  Conv = 1d0
  nSCF = 0
  ECC  = ERHF
  EcCC = 0d0

  n_diis        = 0
  t(:,:)        = 0d0
  t_diis(:,:)   = 0d0
  err_diis(:,:) = 0d0

!------------------------------------------------------------------------
! Main SCF loop
!------------------------------------------------------------------------
  write(*,*)
  write(*,*)'----------------------------------------------------'
  write(*,*)'| pCCD calculation: t amplitudes                   |'
  write(*,*)'----------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A16,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X)') &
            '|','#','|','E(pCCD)','|','Ec(pCCD)','|','Conv','|'
  write(*,*)'----------------------------------------------------'

  do while(Conv > thresh .and. nSCF < maxSCF)

  ! Increment 

    nSCF = nSCF + 1

  ! Form intermediate array
    
   yO(:,:) = matmul(t,transpose(OOVV))
    
   ! Compute residual


    r(:,:) = OOVV(:,:) + 2d0*delta_OV(:,:)*t(:,:) & 
           - 2d0*(2d0*OVOV(:,:) - OVVO(:,:) - OOVV(:,:)*t(:,:))*t(:,:)

    do i=1,O
      do a=1,V

        do j=1,O
          r(i,a) = r(i,a) - 2d0*OOVV(j,a)*t(j,a)*t(i,a) + OOOO(j,i)*t(j,a) + yO(i,j)*t(j,a) 
        end do 

        do b=1,V
          r(i,a) = r(i,a) - 2d0*OOVV(i,b)*t(i,b)*t(i,a) + VVVV(a,b)*t(i,b)
        end do 

      end do
    end do

   ! Check convergence 

    Conv = maxval(abs(r(:,:)))
  
   ! Update amplitudes

   t(:,:) = t(:,:) - 0.5d0*r(:,:)/delta_OV(:,:)

   ! Compute correlation energy

    EcCC = trace_matrix(V,matmul(transpose(OOVV),t))

   ! Dump results

    ECC = ERHF + EcCC

   ! DIIS extrapolation

   if(max_diis > 1) then

     n_diis = min(n_diis+1,max_diis)
     call DIIS_extrapolation(rcond,nO*nV,nO*nV,n_diis,err_diis,t_diis,-0.5d0*r/delta_OV,t)

    end if

    write(*,'(1X,A1,1X,I3,1X,A1,1X,F16.10,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X)') &
      '|',nSCF,'|',ECC+ENuc,'|',EcCC,'|',Conv,'|'

  end do
  write(*,*)'----------------------------------------------------'
!------------------------------------------------------------------------
! End of SCF loop
!------------------------------------------------------------------------

! Did it actually converge?

  if(nSCF == maxSCF) then

    write(*,*)
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)'        Convergence failed for t ampitudes          '
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

    stop
 
  end if

! Deallocate memory

  deallocate(err_diis,t_diis)

! Memory allocation 

  allocate(z(O,V))

! Memory allocation for DIIS

  allocate(err_diis(O*V,max_diis),z_diis(O*V,max_diis))

!------------------------------------------------------------------------
! Compute z ampltiudes
!------------------------------------------------------------------------

  Conv = 1d0
  nSCF = 0

  n_diis          = 0
  z_diis(:,:)     = 0d0
  err_diis(:,:) = 0d0

!------------------------------------------------------------------------
! Main SCF loop
!------------------------------------------------------------------------
  write(*,*)
  write(*,*)'----------------------------------------------------'
  write(*,*)'| pCCD calculation: z amplitudes                   |'
  write(*,*)'----------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A16,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X)') &
            '|','#','|','E(pCCD)','|','Ec(pCCD)','|','Conv','|'
  write(*,*)'----------------------------------------------------'

  do while(Conv > thresh .and. nSCF < maxSCF)

  ! Increment 

    nSCF = nSCF + 1

  ! Form intermediate array
  
   yO(:,:) = matmul(OOVV,transpose(t))
   yV(:,:) = matmul(transpose(OOVV),t)

   ! Compute residual

    r(:,:) = OOVV(:,:) + 2d0*delta_OV(:,:)*z(:,:) &
           - 2d0*(2d0*OVOV(:,:) - OVVO(:,:) - 2d0*OOVV(:,:)*t(:,:))*z(:,:)

    do i=1,O
      do a=1,V

        do j=1,O
          r(i,a) = r(i,a) - 2d0*OOVV(j,a)*t(j,a)*z(i,a) &
                          - 2d0*OOVV(i,a)*z(j,a)*t(j,a) & 
                          + OOOO(i,j)*z(j,a)            &
                          + yO(i,j)*z(j,a)
        end do

        do b=1,V
          r(i,a) = r(i,a) - 2d0*OOVV(i,b)*t(i,b)*z(i,a) & 
                          - 2d0*OOVV(i,a)*z(i,b)*t(i,b) & 
                          + VVVV(b,a)*z(i,b)            & 
                          + yV(a,b)*z(i,b)
        end do

      end do
    end do

   ! Check convergence 

    Conv = maxval(abs(r(:,:)))
  
   ! Update amplitudes

   z(:,:) = z(:,:) - 0.5d0*r(:,:)/delta_OV(:,:)

   ! DIIS extrapolation

   if(max_diis > 1) then

     n_diis = min(n_diis+1,max_diis)
     call DIIS_extrapolation(rcond,O*V,O*V,n_diis,err_diis,z_diis,-0.5d0*r/delta_OV,z)

    end if

    write(*,'(1X,A1,1X,I3,1X,A1,1X,F16.10,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X)') &
      '|',nSCF,'|',ECC+ENuc,'|',EcCC,'|',Conv,'|'

  end do
  write(*,*)'----------------------------------------------------'
  write(*,*)
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
 
  end if

! Deallocate memory

  deallocate(err_diis,z_diis,r)

!--------------------------!
! Compute density matrices !
!--------------------------!

  allocate(xOO(O,O),xVV(V,V),xOV(O,V))

  xOO(:,:) = matmul(t,transpose(z))
  xVV(:,:) = matmul(transpose(z),t)
  xOV(:,:) = matmul(t,matmul(transpose(z),t))

! Form 1RDM

  allocate(rdm1(O+V,O+V))

  rdm1(:,:) = 0d0

  do i=1,O
    rdm1(i,i) = 2d0*(1d0 - xOO(i,i))
  end do

  do a=1,V
    rdm1(O+a,O+a) = 2d0*xVV(a,a)
  end do

! Check 1RDM

  tr_1rdm = trace_matrix(O+V,rdm1)
  write(*,*) ' --> Trace of the 1RDM = ',tr_1rdm

  if( abs(dble(2*O) - tr_1rdm) > thresh ) & 
  write(*,*) ' !!! Your 1RDM seems broken !!! '
  write(*,*)

! Form 2RM

  allocate(rdm2(O+V,O+V,O+V,O+V))

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
      rdm2(i,i,O+a,O+a) = 2d0*(t(i,a) + xOV(i,a) - 2d0*t(i,a)*(xVV(a,a) + xOO(i,i) - t(i,a)*z(i,a)))
    end do
  end do

  ! aaii

  do i=1,O
    do a=1,V
      rdm2(O+a,O+a,i,i) = 2d0*z(a,i)
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
      rdm2(i,O+a,i,O+a) = 4d0*(xVV(a,a) - t(i,a)*z(i,a))
    end do
  end do

  ! iaai

  do i=1,O
    do a=1,V
      rdm2(i,O+a,O+a,i) = - 2d0*(xVV(a,a) - t(i,a)*z(i,a))
    end do
  end do

  ! aiai

  do i=1,O
    do a=1,V
      rdm2(O+a,i,O+a,i) = 4d0*(xVV(a,a) - t(i,a)*z(i,a))
    end do
  end do

  ! aiia

  do i=1,O
    do a=1,V
      rdm2(O+a,i,i,O+a) = - 2d0*(xVV(a,a) - t(i,a)*z(i,a))
    end do
  end do

  ! abab

  do a=1,V
    rdm2(O+a,O+a,O+a,O+a) = 2d0*xVV(a,a)
  end do

! Check 2RDM

  tr_2rdm = trace_matrix((O+V)**2,rdm2)
  write(*,*) ' --> Trace of the 2RDM = ',tr_2rdm

  if( abs(dble(2*O*(2*O-1)) - tr_2rdm) > thresh ) & 
  write(*,*) ' !!! Your 2RDM seems broken !!! '
  write(*,*)

! Testing zone

  if(dotest) then

    call dump_test_value('R','pCCD correlation energy',EcCC)

  end if

end subroutine 
