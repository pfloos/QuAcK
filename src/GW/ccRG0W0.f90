subroutine ccRG0W0(maxSCF,thresh,max_diis,nBas,nOrb,nC,nO,nV,nR,nS,ERI,ENuc,ERHF,eHF)

! CC-based GW module

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
  integer,intent(in)            :: nS
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: ERHF
  double precision,intent(in)   :: eHF(nOrb)
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)

! Local variables

  integer                       :: p,q,r,s
  integer                       :: i,j,k,l
  integer                       :: a,b,c,d
  integer                       :: m

  integer                       :: isp_W
  logical                       :: TDA_W
  logical                       :: dRPA
  integer                       :: nSCF
  double precision              :: Conv
  double precision              :: EcRPA

  double precision,allocatable  :: Sig(:)
  double precision,allocatable  :: Z(:)

  double precision,allocatable  :: del(:,:)
  double precision,allocatable  :: vec(:,:)
  double precision,allocatable  :: res(:,:)
  double precision,allocatable  :: amp(:,:)

  double precision,allocatable  :: Aph(:,:)
  double precision,allocatable  :: Bph(:,:)
  double precision,allocatable  :: Om(:)
  double precision,allocatable  :: XpY(:,:)
  double precision,allocatable  :: XmY(:,:)
  double precision,allocatable  :: rho(:,:,:)

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

  allocate(del(nS,nOrb))
  allocate(vec(nS,nOrb))
  allocate(res(nS,nOrb))
  allocate(amp(nS,nOrb))

  allocate(Sig(nOrb))
  allocate(Z(nOrb))

  allocate(r_diis(nS*nOrb,max_diis))
  allocate(t_diis(nS*nOrb,max_diis))

!-------------------!
! Compute screening !
!-------------------!

  ! Spin manifold 

  isp_W = 1
  TDA_W = .false.
  dRPA  = .true.

  ! Memory allocation

  allocate(Om(nS),Aph(nS,nS),Bph(nS,nS),XpY(nS,nS),XmY(nS,nS),rho(nOrb,nOrb,nS))

  call phLR_A(isp_W,dRPA,nOrb,nC,nO,nV,nR,nS,1d0,eHF,ERI,Aph)
  call phLR_B(isp_W,dRPA,nOrb,nC,nO,nV,nR,nS,1d0,ERI,Bph)

  call phLR(TDA_W,nS,Aph,Bph,EcRPA,Om,XpY,XmY)

  call RGW_excitation_density(nOrb,nC,nO,nR,nS,ERI,XpY,rho)

  deallocate(Aph,Bph,XpY,XmY)

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

    amp(:,:) = 0d0
    res(:,:) = 0d0
 
    ! Compute approximate hessians and coupling blocks
 
    do m=1,nS
      do j=nC+1,nO
        del(m,j) = Om(m) + eHF(j) - eHF(p)
        vec(m,j) = sqrt(2d0)*rho(p,j,m)
      end do
    end do
    
    do m=1,nS
      do b=1,nV-nR
        del(m,nO+b) = Om(m) + eHF(nO+b) - eHF(p)
        vec(m,nO+b) = sqrt(2d0)*rho(p,nO+b,m)
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
 
!     res(:,:) = vec(:,:) + (del(:,:) - Sig(p))*amp(:,:)
 
      do m=1,nS
        do j=nC+1,nO
          res(m,j) = vec(m,j) + (eHF(j) - Om(m) - eHF(p) - Sig(p))*amp(m,j) 
        end do
      end do
 
      ! Compute residual for 2p1h sector
 
      do m=nC+1,nO
        do b=1,nV-nR
          res(m,nO+b) = vec(m,nO+b) + (eHF(nO+b) + Om(m) - eHF(p) - Sig(p))*amp(m,nO+b) 
        end do
      end do
  
      ! Check convergence 
 
      Conv = maxval(abs(res))
    
      ! Update amplitudes

      amp(:,:) = amp(:,:) - res(:,:)/del(:,:)
 
      ! DIIS extrapolation

      if(max_diis > 1) then
        n_diis = min(n_diis+1,max_diis)
        call DIIS_extrapolation(rcond,nS*nOrb,nS*nOrb,n_diis,r_diis,t_diis,res,amp)
      end if

      ! Compute quasiparticle energy
 
      Sig(p) = 0d0

      do m=1,nS
        do q=nC+1,nOrb-nR
          Sig(p) = Sig(p) + vec(m,q)*amp(m,q)
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
