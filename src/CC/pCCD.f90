subroutine pCCD(dotest,maxIt,thresh,max_diis,nBas,nOrb,nC,nO,nV,nR, & 
                Hc,ERI_AO,ENuc,ERHF,eHF,cHF,PHF,FHF)

! pair CCD module

  implicit none

! Input variables

  logical,intent(in)            :: dotest

  integer,intent(in)            :: maxIt
  integer,intent(in)            :: max_diis
  double precision,intent(in)   :: thresh

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: ERHF
  double precision,intent(in)   :: eHF(nOrb)
  double precision,intent(in)   :: cHF(nBas,nOrb)
  double precision,intent(in)   :: PHF(nBas,nBas)
  double precision,intent(in)   :: FHF(nBas,nBas)
  double precision,intent(in)   :: Hc(nBas,nBas)
  double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: mu,nu
  integer                       :: p,q,r,s
  integer                       :: pq,rs
  integer                       :: i,j,a,b

  integer                       :: nItAmp
  integer                       :: nItOrb
  double precision              :: CvgAmp
  double precision              :: CvgOrb
  double precision              :: ECC
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

  double precision,allocatable  :: rdm1(:,:)
  double precision,allocatable  :: rdm2(:,:,:,:)

  double precision              :: E1,E2
  double precision,allocatable  :: c(:,:)
  double precision,allocatable  :: h(:,:)
  double precision,allocatable  :: ERI_MO(:,:,:,:)
  double precision,allocatable  :: grad(:)
  double precision,allocatable  :: hess(:,:)
  double precision,allocatable  :: hessInv(:,:)
  double precision,allocatable  :: Kap(:,:)
  double precision,allocatable  :: ExpKap(:,:)

  integer                       :: O,V,N
  integer                       :: Np
  integer                       :: n_diis
  double precision              :: rcond
  double precision,allocatable  :: err_diis(:,:)
  double precision,allocatable  :: t2_diis(:,:)
  double precision,allocatable  :: z2_diis(:,:)
          
! Hello world

  write(*,*)
  write(*,*)'*******************************'
  write(*,*)'* Restricted pCCD Calculation *'
  write(*,*)'*******************************'
  write(*,*)

! Useful quantities

  O = nO - nC
  V = nV - nR
  N = O + V

  Np = N*N

  !------------------------------------!
  ! Star Loop for orbital optimization !
  !------------------------------------!

  allocate(ERI_MO(N,N,N,N))
  allocate(c(nBas,N),h(N,N))
  allocate(eO(O),eV(V),delta_OV(O,V))
  allocate(OOOO(O,O),OOVV(O,V),OVOV(O,V),OVVO(O,V),VVVV(V,V))

  do p=1,N
    c(:,p) = cHF(:,nC+p)
  enddo

  CvgOrb = 1d0
  nItOrb = 0
  EOld   = ERHF

  write(*,*)
  write(*,*)'---------------------------------------'
  write(*,*)'| Orbital Optimization for pCCD       |'
  write(*,*)'---------------------------------------'

  do while(CvgOrb > thresh .and. nItOrb < maxIt)

    nItOrb = nItOrb + 1

    ! Transform integrals

    h = matmul(transpose(c),matmul(Hc,c))

    call AOtoMO_ERI_RHF(nBas,N,c,ERI_AO,ERI_MO)

    ! Form energy denominator

    eO(:) = 0d0
    eV(:) = 0d0

    do mu=1,nBas
      do nu=1,nBas

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
        call DIIS_extrapolation(rcond,nO*nV,nO*nV,n_diis,err_diis,t2_diis,-0.5d0*r2/delta_OV,t2)
 
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
 
    !--------------------------!
    ! Compute orbital gradient !
    !--------------------------!
 
    allocate(grad(Np))
 
    call pCCD_orbital_gradient(O,V,N,Np,h,ERI_MO,rdm1,rdm2,grad)

   ! Check convergence of orbital optimization
 
    CvgOrb = maxval(abs(grad))
    write(*,*) '----------------------------------------------------------'
    write(*,'(A10,I4,A30)') ' Iteration',nItOrb,'for pCCD orbital optimization'
    write(*,*) '----------------------------------------------------------'
    write(*,'(A40,F16.10,A3)') ' Convergence of orbital gradient = ',CvgOrb,' au'
    write(*,'(A40,F16.10,A3)') ' Energy difference = ',ECC-EOld,' au'
    write(*,*) '----------------------------------------------------------'
    write(*,*)

    EOld = ECC

    !-------------------------!
    ! Compute orbital Hessian !
    !-------------------------!
 
    allocate(hess(Np,Np))

    call pCCD_orbital_hessian(O,V,N,Np,h,ERI_MO,rdm1,rdm2,hess)

    deallocate(rdm1,rdm2)
 
    allocate(hessInv(Np,Np))
 
    call inverse_matrix(Np,hess,hessInv)

    deallocate(hess)
 
    allocate(Kap(N,N))
 
    Kap(:,:) = 0d0
 
    pq = 0
    do p=1,N
      do q=1,N
 
        pq = pq + 1
 
        rs = 0
        do r=1,N
          do s=1,N
 
            rs = rs + 1
 
              Kap(p,q) = Kap(p,q) - hessInv(pq,rs)*grad(rs)
 
          end do
        end do
 
      end do
    end do
 
    deallocate(hessInv,grad)

    write(*,*) 'kappa'
    call matout(N,N,Kap)
    write(*,*)
 
    allocate(ExpKap(N,N))
    call matrix_exponential(N,Kap,ExpKap)
    deallocate(Kap)
 
    write(*,*) 'e^kappa'
    call matout(N,N,ExpKap)
    write(*,*)
 
    write(*,*) 'Old orbitals'
    call matout(nBas,N,c)
    write(*,*)
 
    c = matmul(c,ExpKap)
    deallocate(ExpKap)
 
    write(*,*) 'Rotated orbitals'
    call matout(nBas,N,c)
    write(*,*)

  end do

  !-----------------------------------!
  ! End Loop for orbital optimization !
  !-----------------------------------!

    ! Did it actually converge?

    if(nItOrb == maxIt) then

      write(*,*)
      write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      write(*,*)'! Convergence failed for orbital optimization !'
      write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

      stop

    end if

! Testing zone

  if(dotest) then

    call dump_test_value('R','pCCD correlation energy',EcCC)

  end if

end subroutine 
