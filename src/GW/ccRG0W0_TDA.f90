subroutine ccRG0W0_TDA(maxSCF,thresh,max_diis,nBas,nOrb,nC,nO,nV,nR,ERI,ENuc,ERHF,eHF)

! CC-based GW module (TDA screening)

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: maxSCF
  double precision,intent(in)   :: thresh
  integer,intent(in)            :: max_diis

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: ERHF
  double precision,intent(in)   :: eHF(nOrb)
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)

! Local variables

  integer                       :: p,q,r,s
  integer                       :: i,j,k,l
  integer                       :: a,b,c,d

  integer                       :: nSCF
  double precision              :: Conv

  double precision,allocatable  :: Sig(:)
  double precision,allocatable  :: Z(:)

  double precision,allocatable  :: del(:,:,:)
  double precision,allocatable  :: vec(:,:,:)
  double precision,allocatable  :: res(:,:,:)
  double precision,allocatable  :: amp(:,:,:)

  integer                       :: n_diis
  double precision              :: rcond
  double precision,allocatable  :: r_diis(:,:)
  double precision,allocatable  :: t_diis(:,:)


! Hello world

  write(*,*)
  write(*,*)'*****************************'
  write(*,*)'* CC-based G0W0 Calculation *'
  write(*,*)'*****************************'
  write(*,*)

! Memory allocation

  allocate(del(nO,nV,nOrb))
  allocate(vec(nO,nV,nOrb))
  allocate(res(nO,nV,nOrb))
  allocate(amp(nO,nV,nOrb))

  allocate(Sig(nOrb))
  allocate(Z(nOrb))

  allocate(r_diis(nO*nV*nOrb,max_diis))
  allocate(t_diis(nO*nV*nOrb,max_diis))

! Initialization

  Sig(:) = 0d0
  Z(:)   = 1d0

!-------------------------!
! Main loop over orbitals !
!-------------------------!

  do p=nO,nO+1

    ! Initialization
 
    Conv = 1d0
    nSCF =  0

    n_diis      = 0
    t_diis(:,:) = 0d0
    r_diis(:,:) = 0d0
    rcond       = 0d0

    amp(:,:,:) = 0d0
    res(:,:,:) = 0d0
 
  ! Compute energy differences
 
    do i=nC+1,nO
      do j=nC+1,nO
        do a=1,nV-nR
    
          del(i,a,j) = eHF(i) + eHF(j) - eHF(nO+a) - eHF(p)
    
        end do
      end do
    end do
    
    do i=nC+1,nO
      do a=1,nV-nR
        do b=1,nV-nR
    
          del(i,a,nO+b) = eHF(nO+a) + eHF(nO+b) - eHF(i) - eHF(p)
    
        end do
      end do
    end do

    do i=nC+1,nO
      do a=1,nV-nR
        do j=nC+1,nO
 
          vec(i,a,j) = sqrt(2d0)*ERI(p,nO+a,i,j)
 
        end do
      end do
    end do
 
    do i=nC+1,nO
      do a=1,nV-nR
        do b=1,nV-nR
 
          vec(i,a,nO+b) = sqrt(2d0)*ERI(p,i,nO+b,nO+a)
 
        end do
      end do
    end do

   !----------------------!
   ! Loop over amplitudes !
   !----------------------!

    write(*,*)
    write(*,*)'-------------------------------------------------------------'
    write(*,*)'| CC-based G0W0 calculation                                 |'
    write(*,*)'-------------------------------------------------------------'
    write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X)') &
              '|','#','|','Sig_c (eV)','|','e_GW (eV)','|','Conv','|'
    write(*,*)'-------------------------------------------------------------'
 
    do while(Conv > thresh .and. nSCF < maxSCF)
 
     ! Increment 
 
      nSCF = nSCF + 1
 
      ! Compute residual for 2h1p sector
 
      res(:,:,:) = vec(:,:,:) + (del(:,:,:) - Sig(p))*amp(:,:,:)
 
      do i=nC+1,nO
        do a=1,nV-nR
          do j=nC+1,nO
 
            do k=nC+1,nO
              do c=1,nV-nR
 
                res(i,a,j) = res(i,a,j) - 2d0*ERI(j,nO+c,nO+a,k)*amp(i,c,k)  
!                                       - 2d0*ERI(i,nO+c,nO+a,k)*amp(k,c,j)
 
              end do
            end do
 
          end do
        end do
      end do
 
      ! Compute residual for 2p1h sector
 
      do i=nC+1,nO
        do a=1,nV-nR
          do b=1,nV-nR
 
            do k=nC+1,nO
              do c=1,nV-nR
 
                res(i,a,nO+b) = res(i,a,nO+b) + 2d0*ERI(nO+a,k,i,nO+c)*amp(k,c,nO+b) 
!                                             + 2d0*ERI(nO+b,k,i,nO+c)*amp(k,a,nO+c) 
 
              end do
            end do
 
          end do
        end do
      end do
  
      ! Check convergence 
 
      Conv = maxval(abs(res))
    
      ! Update amplitudes

      amp(:,:,:) = amp(:,:,:) - res(:,:,:)/del(:,:,:)
 
      ! DIIS extrapolation

      if(max_diis > 1) then
     
        n_diis = min(n_diis+1,max_diis)
        call DIIS_extrapolation(rcond,nO*nV*nOrb,nO*nV*nOrb,n_diis,r_diis,t_diis,res,amp)
     
      end if

      ! Compute quasiparticle energy
 
      Sig(p) = 0d0

      do i=nC+1,nO
        do a=1,nV-nR
          do q=nC+1,nOrb-nR
 
            Sig(p) = Sig(p) + vec(i,a,q)*amp(i,a,q)
  
          end do
        end do
      end do

      ! Dump results
 
      write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.10,1X,A1,1X,F15.10,1X,A1,1X,F15.10,1X,A1,1X)') &
        '|',nSCF,'|',Sig(p)*HaToeV,'|',(eHF(p)+Sig(p))*HaToeV,'|',Conv,'|'
 
    end do

    write(*,*)'-------------------------------------------------------------'
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
      write(*,*)
 
    end if

    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)'| CC-based G0W0 calculation                                                   |'
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X)') &
              '|','Orb','|','e_HF (eV)','|','Sig_c (eV)','|','Z','|','e_QP (eV)','|'
    write(*,*)'-------------------------------------------------------------------------------'
 
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.10,1X,A1,1X,F15.10,1X,A1,1X,F15.10,1X,A1,1X,F15.10,1X,A1,1X)') &
    '|',p,'|',eHF(p)*HaToeV,'|',Sig(p)*HaToeV,'|',Z(p),'|',(eHF(p)+Sig(p))*HaToeV,'|'
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)

  end do

end subroutine 
