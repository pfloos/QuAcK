subroutine R_IP_ADC_2SOSEX(dotest,sig_inf,TDA_W,flow,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,eHF)

! Non-Dyson version of ADC-2SOSEX for IPs

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dotest

  logical,intent(in)            :: sig_inf
  logical,intent(in)            :: TDA_W
  double precision,intent(in)   :: flow
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: ERHF
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: eHF(nOrb)

! Local variables

  integer                       :: p,q,r,s
  integer                       :: i,j,k,l
  integer                       :: a,b,c,d
  integer                       :: mu,nu
  integer                       :: klc,kcd,ija,iab
  double precision              :: num,dem,reg
  double precision              :: omega

  logical                       :: print_W = .false.
  logical                       :: dRPA
  integer                       :: ispin
  double precision              :: EcRPA
  integer                       :: n2h1p,n2p1h,nH
  double precision,external     :: Kronecker_delta
  double precision,allocatable  :: H(:,:)
  double precision,allocatable  :: eGW(:)
  double precision,allocatable  :: Z(:)
  double precision,allocatable  :: Aph(:,:)
  double precision,allocatable  :: Bph(:,:)
  double precision,allocatable  :: Om(:)
  double precision,allocatable  :: XpY(:,:)
  double precision,allocatable  :: XmY(:,:)
  double precision,allocatable  :: rho(:,:,:)

  double precision,allocatable  :: F(:,:)
  double precision,allocatable  :: Vh(:,:)
  double precision,allocatable  :: Vx(:,:)
  double precision,allocatable  :: DM(:,:)
  double precision,allocatable  :: w(:,:,:)

  double precision,allocatable  :: U_2p1h(:,:)
  double precision,allocatable  :: K_2p1h(:,:)

  logical,parameter             :: verbose = .false.
  double precision,parameter    :: cutoff1 = 0.1d0
  double precision,parameter    :: cutoff2 = 0.01d0
  double precision              :: eF
  double precision,parameter    :: window = 1.5d0

  double precision              :: start_timing,end_timing,timing

! Output variables

! Hello world

  write(*,*)
  write(*,*)'****************************************'
  write(*,*)'* Restricted IP-ADC-2SOSEX Calculation *'
  write(*,*)'****************************************'
  write(*,*)

! Dimension of the supermatrix

  n2h1p = nO*nO*nV
  n2p1h = nV*nV*nO
  nH = nO + n2h1p

! Memory allocation

  allocate(H(nH,nH),eGW(nH),Z(nH))
  allocate(U_2p1h(n2p1h,nO),K_2p1h(n2p1h,n2p1h))

! Initialization

  dRPA = .true.
  EcRPA = 0d0

  eF = 0.5d0*(eHF(nO+1) + eHF(nO))

!-------------------!
! Compute screening !
!-------------------!

  ! Memory allocation 

  allocate(Om(nS),Aph(nS,nS),Bph(nS,nS),XpY(nS,nS),XmY(nS,nS),rho(nOrb,nOrb,nS))

  ! Spin manifold 

  ispin = 1

  call phRLR_A(ispin,dRPA,nOrb,nC,nO,nV,nR,nS,1d0,eHF,ERI,Aph)
  call phRLR_B(ispin,dRPA,nOrb,nC,nO,nV,nR,nS,1d0,ERI,Bph)

  call phRLR(TDA_W,nS,Aph,Bph,EcRPA,Om,XpY,XmY)

  ! Small shift to avoid hard zeros in amplitudes

  Om(:) = Om(:) + 1d-12

  if(print_W) call print_excitation_energies('phRPA@RHF','singlet',nS,Om)

  !--------------------------!
  ! Compute spectral weights !
  !--------------------------!

  call RGW_excitation_density(nOrb,nC,nO,nR,nS,ERI,XpY,rho)

  deallocate(Aph,Bph)

  !-------------------!
  ! Compute Sigma(oo) !
  !-------------------!

! allocate(F(nOrb,nOrb))
! F(:,:) = 0d0

  if(sig_inf) then

!   allocate(DM(nOrb,nOrb),Vh(nOrb,nOrb),Vx(nOrb,nOrb),w(nOrb,nOrb,nS))

!   ! call R_linDM_GW(nOrb,nC,nO,nV,nR,nS,eHF,Om,rho,0d0,DM)
!   call R_linDM_2SOSEX(nOrb,nC,nO,nV,nR,nS,eHF,Om,rho,ERI,0d0,DM)
!   call Hartree_matrix_AO_basis(nOrb,DM,ERI,Vh)
!   call exchange_matrix_AO_basis(nOrb,DM,ERI,Vx)
!
!   F(:,:) = Vh(:,:) + 0.5d0*Vx(:,:)
!
!   deallocate(Vh,Vx,DM,w,XpY,XmY)

  end if

! Initialization

  H(:,:) = 0d0

  !------------------------------!
  ! Compute IP-ADC-2SOSEX matrix !
  !------------------------------!
  !                              !
  !     | F      U_2h1p     |    ! 
  ! H = |                   |    ! 
  !     | U_2h1p (K+C)_2h1p |    ! 
  !                              !
  !------------------------------!

  call wall_time(start_timing)

  !---------!
  ! Block F !
  !---------!

  do i=nC+1,nO

    H(i,i) = eHF(i)

  end do

  !--------------!
  ! Block U_2p1h !
  !--------------!
  
  do j=nC+1,nO
  
    iab = 0
    do a=nO+1,nOrb-nR
      do mu=1,nS
        iab = iab + 1
  
        U_2p1h(iab,j) = sqrt(2d0)*rho(j,a,mu)
  
        do k=nC+1,nO
          do c=nO+1,nOrb-nR
  
            num = sqrt(2d0)*ERI(j,k,c,a)*rho(c,k,mu)
            dem = eHF(c) - eHF(k) - Om(mu)
            reg = (1d0 - exp(-2d0*flow*dem*dem))/dem
  
            U_2p1h(iab,j) = U_2p1h(iab,j) + num*reg
  
            num = sqrt(2d0)*ERI(j,c,k,a)*rho(k,c,mu)
            dem = eHF(c) - eHF(k) + Om(mu)
            reg = (1d0 - exp(-2d0*flow*dem*dem))/dem
  
            U_2p1h(iab,j) = U_2p1h(iab,j) + num*reg
  
          end do
        end do
  
      end do
    end do
  
  end do
  
  !------------------!
  ! Block (K+C)_2p1h !
  !------------------!

  do i=nC+1,nO
    do j=nC+1,nO

      omega = 0.5d0*(eHF(i) + eHF(j))
  
      K_2p1h(:,:) = 0d0
  
      iab = 0
      do a=nO+1,nOrb-nR
        do mu=1,nS
          iab = iab + 1
  
          K_2p1h(iab,iab) = 1d0/(omega - eHF(a) - Om(mu))
  
        end do
      end do

      H(i,j) = H(i,j) + dot_product(U_2p1h(:,i),matmul(K_2p1h,U_2p1h(:,j)))

    end do
  end do

  !--------------!
  ! Block U_2h1p !
  !--------------!

  do j=nC+1,nO

    ija = 0
    do i=nC+1,nO
      do mu=1,nS
        ija = ija + 1

        H(j     ,nO+ija) = sqrt(2d0)*rho(j,i,mu)
        H(nO+ija,j     ) = sqrt(2d0)*rho(j,i,mu)

        do k=nC+1,nO
          do c=nO+1,nO

            num = sqrt(2d0)*ERI(j,c,k,i)*rho(k,c,mu)
            dem = eHF(c) - eHF(k) - Om(mu)
            reg = (1d0 - exp(-2d0*flow*dem*dem))/dem

            H(j     ,nO+ija) = H(j     ,nO+ija) + num*reg
            H(nO+ija,j     ) = H(nO+ija,j     ) + num*reg

            num = sqrt(2d0)*ERI(j,k,c,i)*rho(c,k,mu)
            dem = eHF(c) - eHF(k) + Om(mu)
            reg = (1d0 - exp(-2d0*flow*dem*dem))/dem

            H(j     ,nO+ija) = H(j     ,nO+ija) + num*reg
            H(nO+ija,j     ) = H(nO+ija,j     ) + num*reg

          end do
        end do

      end do
    end do

  end do

  !------------------!
  ! Block (K+C)_2h1p !
  !------------------!

  ija = 0
  do i=nC+1,nO
    do mu=1,nS
      ija = ija + 1

      H(nO+ija,nO+ija) = eHF(i) - Om(mu) 

    end do
  end do

  call wall_time(end_timing)

  timing = end_timing - start_timing
  write(*,*)
  write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for construction of supermatrix = ',timing,' seconds'
  write(*,*)

!-------------------------!
! Diagonalize supermatrix !
!-------------------------!

  call wall_time(start_timing)

  call diagonalize_matrix(nH,H,eGW)

  call wall_time(end_timing)

  timing = end_timing - start_timing
  write(*,*)
  write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for diagonalization of supermatrix = ',timing,' seconds'
  write(*,*)

!-----------------!
! Compute weights !
!-----------------!

  Z(:) = 0d0
  do s=1,nH
    do i=nC+1,nO
      Z(s) = Z(s) + H(i,s)**2
    end do
  end do

!--------------!
! Dump results !
!--------------!

  write(*,*)'---------------------------------------------'
  write(*,'(1X,A45)')'| ADC-2SOSEX energies for all orbitals      |'
  write(*,*)'---------------------------------------------'
  write(*,'(1X,A1,1X,A5,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X)') &
            '|','#','|','e_QP (eV)','|','Z','|'
  write(*,*)'---------------------------------------------'

  do s=1,nH
!   if(eGW(s) < eF .and. eGW(s) > eF - window) then
    if(Z(s) > cutoff1) then
      write(*,'(1X,A1,1X,I5,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
      '|',s,'|',eGW(s)*HaToeV,'|',Z(s),'|'
    end if
  end do

  write(*,*)'---------------------------------------------'
  write(*,*)

  if(verbose) then 
  
    do s=1,nH

      if(eGW(s) < eF .and. eGW(s) > eF - window) then

        write(*,*)'-------------------------------------------------------------'
        write(*,'(1X,A12,1X,I3,A1,1X,A7,F12.6,A13,F6.4,1X)') & 
         'Eigenvalue #',s,':','e_QP = ',eGW(s)*HaToeV,' eV and Z = ',Z(s)
        write(*,*)'-------------------------------------------------------------'
        write(*,'(1X,A20,1X,A20,1X,A15,1X)') &
                  ' Conf. (p,mu)  ',' Coefficient ',' Weight ' 
        write(*,*)'-------------------------------------------------------------'
      
        do p=nC+1,nO 
          if(abs(H(p,s)) > cutoff2)                     &
            write(*,'(1X,A7,I3,A16,1X,F15.6,1X,F15.6)') &
            '      (',p,'    )           ',H(p,s),H(p,s)**2
        end do

        ija = 0
        do i=nC+1,nO
          do mu=1,nS
            ija = ija + 1
 
            if(abs(H(nO+ija,s)) > cutoff2)                  &
            write(*,'(1X,A7,I3,A1,I3,A12,1X,F15.6,1X,F15.6)') &
            '      (',i,',',mu,')           ',H(nO+ija,s),H(nO+ija,s)**2
       
          end do
        end do
       
        write(*,*)'-------------------------------------------------------------'
        write(*,*)

      end if

    end do

  end if

end subroutine 
