subroutine R_diag_half_GW(dotest,sig_inf,TDA_W,flow,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,eHF)

! Half-and-half version of GW in the diagonal approximation

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
  integer                       :: mu
  integer                       :: klc,kcd,ija,iab
  double precision              :: num,dem

  logical                       :: print_W = .false.
  logical                       :: dRPA
  integer                       :: ispin
  double precision              :: EcRPA
  integer                       :: n2h1p,nH
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

  logical,parameter             :: verbose = .false.
  double precision,parameter    :: cutoff1 = 0.1d0
  double precision,parameter    :: cutoff2 = 0.01d0
  double precision              :: eF
  double precision,parameter    :: window = 1.5d0

  double precision              :: start_timing,end_timing,timing

! Output variables

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'* Restricted diag-half-and-half-GW Calculation *'
  write(*,*)'************************************************'
  write(*,*)

! Diagonal approximation

  write(*,*)' Diagonal approximation enforced!'
  write(*,*)

! Dimension of the supermatrix

  n2h1p = nO*nO*nV
  nH = 1 + n2h1p

! Memory allocation

  allocate(H(nH,nH),eGW(nH),Z(nH))

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

  if(print_W) call print_excitation_energies('phRPA@RHF','singlet',nS,Om)

  !--------------------------!
  ! Compute spectral weights !
  !--------------------------!

  call RGW_excitation_density(nOrb,nC,nO,nR,nS,ERI,XpY,rho)

  deallocate(Aph,Bph,XpY,XmY)

  !-------------------!
  ! Compute Sigma(oo) !
  !-------------------!

  allocate(F(nOrb,nOrb))
  F(:,:) = 0d0

  if(sig_inf) then

    allocate(DM(nOrb,nOrb),Vh(nOrb,nOrb),Vx(nOrb,nOrb))
 
    call R_linDM_GW(nOrb,nC,nO,nV,nR,nS,eHF,Om,rho,0d0,DM)
    call Hartree_matrix_AO_basis(nOrb,DM,ERI,Vh)
    call exchange_matrix_AO_basis(nOrb,DM,ERI,Vx)
 
    F(:,:) = Vh(:,:) + 0.5d0*Vx(:,:)
 
    deallocate(Vh,Vx,DM)

  end if

!-------------------------!
! Main loop over orbitals !
!-------------------------!

  do i=nC+1,nO

    H(:,:) = 0d0

    !------------------------------!
    ! Compute IP-ADC-GW matrix     !
    !------------------------------!
    !                              !
    !     | F      U_2h1p        | !
    ! H = |                      | ! 
    !     | U_2h1p (K+C)_2h1p    | !
    !                              !
    !------------------------------!

    call wall_time(start_timing)
 
    !---------!
    ! Block F !
    !---------!
 
    H(1,1) = eHF(i) + F(i,i)
 
    !-------------------!
    ! Block static 2p1h !
    !-------------------!
 
    do mu=1,nS
      do a=nO+1,nOrb-nR
 
        num = 2d0*rho(i,a,mu)*rho(i,a,mu)
        dem = eHF(i) - eHF(a) - Om(mu)
        H(1,1) = H(1,1) + num/dem
 
      end do
    end do
 
    !--------------!
    ! Block U_2h1p !
    !--------------!
 
    ija = 0
    do k=nC+1,nO
      do mu=1,nS
        ija = ija + 1
 
        H(1       ,1+ija) = sqrt(2d0)*rho(i,k,mu)
        H(1+ija,1       ) = sqrt(2d0)*rho(i,k,mu)
 
      end do
    end do
 
    !------------------!
    ! Block (K+C)_2h1p !
    !------------------!
 
    ija = 0
    do k=nC+1,nO
      do mu=1,nS
        ija = ija + 1
 
        H(1+ija,1+ija) = eHF(k) - Om(mu) 
 
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
      Z(s) = Z(s) + H(1,s)**2
    end do

  !--------------!
  ! Dump results !
  !--------------!

    write(*,*)'-------------------------------------------'
    write(*,'(1X,A29,I3,A11)')'| ADC-GW energies for orbital',i,'|'
    write(*,*)'-------------------------------------------'
    write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X)') &
              '|','#','|','e_QP (eV)','|','Z','|'
    write(*,*)'-------------------------------------------'

    do s=1,nH
!     if(eGW(s) < eF .and. eGW(s) > eF - window) then
      if(Z(s) > cutoff1) then
        write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
        '|',s,'|',eGW(s)*HaToeV,'|',Z(s),'|'
      end if
    end do

    write(*,*)'-------------------------------------------'
    write(*,*)

    if(verbose) then

      do s=1,nH

        if(eGW(s) < eF .and. eGW(s) > eF - window) then

          write(*,*)'------------------------------------------------------------------------------'
          write(*,'(1X,A7,1X,I3,A6,I3,A1,1X,A7,F12.6,A13,F6.4,1X)') &
           'Orbital',i,' and #',s,':','e_QP = ',eGW(s)*HaToeV,' eV and Z = ',Z(s)
          write(*,*)'------------------------------------------------------------------------------'
          write(*,'(1X,A20,1X,A20,1X,A15,1X)') &
                    ' Conf. (i,ia)  ',' Coefficient ',' Weight '
          write(*,*)'------------------------------------------------------------------------------'

          if(i <= nO) &
            write(*,'(1X,A7,I3,A16,1X,F15.6,1X,F15.6,1X,F12.6)') &
            '      (',i,')               ',H(1,s),H(1,s)**2,-eHF(i)*HaToeV

          do k=nC+1,nO
            do ija=1,nS

              if(abs(H(1+ija,s)) > cutoff2)                     &
              write(*,'(1X,A7,I3,A1,I3,A12,1X,F15.6,1X,F15.6,1X,F12.6)') &
              '      (',k,',',ija,')           ',H(1+ija,s),H(1+ija,s)**2,(eHF(i) - Om(ija))*HaToeV

            end do
          end do

          write(*,*)'------------------------------------------------------------------------------'
          write(*,*)

        end if ! If state s should be print

      end do ! Loop on s

    end if ! If verbose

  end do

end subroutine 
