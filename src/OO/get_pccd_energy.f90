subroutine get_pccd_energy(O,V,N,nS,Hc,c,ERI_AO,ECC)
      
! Compute pccd energy for given orbitals c

  implicit none
  include 'parameters.h'
      
! Input variables

  integer,intent(in)            :: O
  integer,intent(in)            :: V
  integer,intent(in)            :: N
  integer,intent(in)            :: nS
  double precision,intent(in)   :: Hc(N,N),c(N,N)
  double precision,intent(in)   :: ERI_AO(N,N,N,N)

! Local variables
  double precision,allocatable  :: rdm1(:,:)
  double precision,allocatable  :: rdm2(:,:,:,:)
  double precision,allocatable  :: PHF(:,:),JHF(:,:),FHF(:,:),K(:,:),h(:,:)
  double precision,allocatable  :: ERI_MO(:,:,:,:)

  integer                       :: mu,nu
  integer                       :: p,q,r,s
  integer                       :: pq,rs
  integer                       :: i,j,a,b
  integer                       :: max_diis = 5

  integer                       :: nItAmp
  integer                       :: nItOrb
  double precision              :: CvgAmp
  double precision              :: CvgOrb
  double precision              :: ENuc = 0d0
  double precision              :: EcCC
  double precision              :: EOld

  double precision,allocatable  :: eO(:)
  double precision,allocatable  :: eV(:)
  double precision,allocatable  :: delta_OV(:,:)

  double precision,allocatable  :: OOOO(:,:)
  double precision,allocatable  :: OOVV(:,:)
  double precision,allocatable  :: OVOV(:,:)
  double precision,allocatable  :: OVVO(:,:)
  double precision,allocatable  :: VVVV(:,:)

  double precision,allocatable  :: r2(:,:)
  double precision,allocatable  :: t2(:,:)
  double precision,allocatable  :: z2(:,:)

  integer                       :: n_diis
  double precision              :: rcond
  double precision              :: thresh = 1d-3
  integer                       :: maxIt = 256 
  double precision,allocatable  :: err_diis(:,:)
  double precision,allocatable  :: t2_diis(:,:)
  double precision,allocatable  :: z2_diis(:,:)

! Output variables

  double precision,intent(out)  :: ECC

allocate(rdm1(N,N), rdm2(N,N,N,N), PHF(N,N), JHF(N,N), K(N,N), FHF(N,N), &
         ERI_MO(N,N,N,N), h(N,N))
allocate(eO(O),eV(V),delta_OV(O,V))
allocate(OOOO(O,O),OOVV(O,V),OVOV(O,V),OVVO(O,V),VVVV(V,V))

  ! Compute Fock operator and AO to MO trafos
  PHF(:,:) = 2d0 * matmul(c(:,1:O), transpose(c(:,1:O))) 
  JHF(:,:) = 0d0
  K(:,:) = 0d0
  call Hartree_matrix_AO_basis(N,PHF,ERI_AO,JHF)
  call exchange_matrix_AO_basis(N,PHF,ERI_AO,K)
  FHF(:,:) = Hc(:,:) + JHF(:,:) + 0.5d0*K(:,:)
  call AOtoMO(N,N,c,Hc,h)
  call AOtoMO_ERI_RHF(N,N,c,ERI_AO,ERI_MO)

  ! Form energy denominator

  eO(:) = 0d0
  eV(:) = 0d0

  do mu=1,N
    do nu=1,N

      do i=1,O
        eO(i) = eO(i) + c(mu,i)*FHF(mu,nu)*c(nu,i)
      end do

      do a=1,V
        eV(a) = eV(a) + c(mu,O+a)*FHF(mu,nu)*c(nu,O+a)
      end do

    end do
  end do

  do i=1,O
    do a=1,V
      delta_OV(i,a) = eV(a) - eO(i)
    end do
  end do

  ! Create integral batches

  do i=1,O
    do j=1,O
      OOOO(i,j) = ERI_MO(i,i,j,j)
    end do
  end do

  do i=1,O
    do a=1,V
      OOVV(i,a) = ERI_MO(i,i,O+a,O+a)
      OVOV(i,a) = ERI_MO(i,O+a,i,O+a)
      OVVO(i,a) = ERI_MO(i,O+a,O+a,i)
    end do
  end do

  do a=1,V
    do b=1,V
      VVVV(a,b) = ERI_MO(O+a,O+a,O+b,O+b)
    end do
  end do

  !----------------------------!
  ! Star Loop for t amplitudes !
  !----------------------------!

  allocate(t2(O,V),r2(O,V))
  allocate(err_diis(O*V,max_diis),t2_diis(O*V,max_diis))

  CvgAmp = 1d0
  nItAmp = 0

  n_diis        = 0
  t2(:,:)       = 0d0
  t2_diis(:,:)  = 0d0
  err_diis(:,:) = 0d0

  write(*,*)
  write(*,*)'---------------------------------------'
  write(*,*)'| pCCD calculation: t amplitudes      |'
  write(*,*)'---------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A16,1X,A1,1X,A10,1X,A1,1X)') &
            '|','#','|','Ec(pCCD)','|','Conv','|'
  write(*,*)'---------------------------------------'

  do while(CvgAmp > thresh .and. nItAmp < maxIt)

    ! Increment 

    nItAmp = nItAmp + 1

    ! Compute residual for t amplitudes

    call pCCD_t_residual(O,V,N,OOVV,OVOV,OVVO,OOOO,VVVV,delta_OV,t2,r2)

    ! Check convergence 

    CvgAmp = maxval(abs(r2(:,:)))
  
    ! Update amplitudes

    t2(:,:) = t2(:,:) - 0.5d0*r2(:,:)/delta_OV(:,:)

    ! Compute correlation energy from t amplitudes

    EcCC = 0d0
    do i=1,O
      do a=1,V
        EcCC = EcCC + OOVV(i,a)*t2(i,a) 
      end do
    end do

    ! DIIS extrapolation

    if(max_diis > 1) then

      n_diis = min(n_diis+1,max_diis)
      call DIIS_extrapolation(rcond,O*V,O*V,n_diis,err_diis,t2_diis,-0.5d0*r2/delta_OV,t2)

    end if

    write(*,'(1X,A1,1X,I3,1X,A1,1X,F16.10,1X,A1,1X,F10.6,1X,A1,1X)') &
      '|',nItAmp,'|',EcCC,'|',CvgAmp,'|'

  end do
  write(*,*)'---------------------------------------'

  !---------------------------!
  ! End Loop for t amplitudes !
  !---------------------------!

  deallocate(r2)
  deallocate(err_diis,t2_diis)

  ! Did it actually converge?

  if(nItAmp == maxIt) then

    write(*,*)
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)'! Convergence failed for t ampitudes !'
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

    stop

  end if

  !-----------------------------!
  ! Start Loop for z amplitudes !
  !-----------------------------!

  allocate(z2(O,V),r2(O,V))
  allocate(err_diis(O*V,max_diis),z2_diis(O*V,max_diis))

  CvgAmp = 1d0
  nItAmp = 0

  n_diis        = 0
  z2_diis(:,:)  = 0d0
  err_diis(:,:) = 0d0

  write(*,*)
  write(*,*)'---------------------------------------'
  write(*,*)'| pCCD calculation: z amplitudes      |'
  write(*,*)'---------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A16,1X,A1,1X,A10,1X,A1,1X)') &
            '|','#','|','Ec(pCCD)','|','Conv','|'
  write(*,*)'---------------------------------------'

  do while(CvgAmp > thresh .and. nItAmp < maxIt)

   ! Increment 

    nItAmp = nItAmp + 1

   ! Compute residual for the z amplitudes

   call pCCD_z_residual(O,V,N,OOVV,OVOV,OVVO,OOOO,VVVV,delta_OV,t2,z2,r2)

   ! Check convergence 

    CvgAmp = maxval(abs(r2(:,:)))
  
   ! Update amplitudes

   z2(:,:) = z2(:,:) - 0.5d0*r2(:,:)/delta_OV(:,:)

    ! Compute correlation energy

    EcCC = 0d0
    do i=1,O
      do a=1,V
        EcCC = EcCC + OOVV(i,a)*z2(i,a)
      end do
    end do

   ! DIIS extrapolation

   if(max_diis > 1) then

     n_diis = min(n_diis+1,max_diis)
     call DIIS_extrapolation(rcond,O*V,O*V,n_diis,err_diis,z2_diis,-0.5d0*r2/delta_OV,z2)

    end if

    write(*,'(1X,A1,1X,I3,1X,A1,1X,F16.10,1X,A1,1X,F10.6,1X,A1,1X)') &
      '|',nItAmp,'|',EcCC,'|',CvgAmp,'|'

  end do
  write(*,*)'---------------------------------------'
  write(*,*)

  !---------------------------!
  ! End Loop for z ampltiudes !
  !---------------------------!

  deallocate(r2)
  deallocate(err_diis,z2_diis)

  ! Did it actually converge?

  if(nItAmp == maxIt) then

    write(*,*)
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)'! Convergence failed for z ampltiudes !'
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

    stop

  end if

  !--------------------------!
  ! Compute density matrices !
  !--------------------------!

  allocate(rdm1(N,N),rdm2(N,N,N,N))

  call pCCD_rdm(O,V,N,ENuc,h,ERI_MO,t2,z2,rdm1,rdm2,ECC)

  deallocate(t2,z2)


end subroutine
