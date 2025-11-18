subroutine R_ADC_2SOSEX_diag(dotest,TDA_W,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,eHF)

! ADC version of 2SOSEX within the diagonal approximation

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dotest

  logical,intent(in)            :: TDA_W
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

  integer                       :: p
  integer                       :: s
  integer                       :: i,j,k,l
  integer                       :: a,b,c,d
  integer                       :: jb,kc,ia,ja
  integer                       :: mu,nu
  integer                       :: klc,kcd,ija,ijb,iab,jab

  logical                       :: print_W = .false.
  logical                       :: dRPA = .true.
  integer                       :: isp_W
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

  logical                       :: verbose = .false.
  double precision,parameter    :: cutoff1 = 0.01d0
  double precision,parameter    :: cutoff2 = 0.01d0
  double precision              :: eF
  double precision,parameter    :: window = 2.5d0

  double precision              :: start_timing,end_timing,timing

! Output variables

! Hello world

  write(*,*)
  write(*,*)'*************************************'
  write(*,*)'* Restricted ADC-2SOSEX Calculation *'
  write(*,*)'*************************************'
  write(*,*)

! Diagonal approximation

  write(*,*)' Diagonal approximation enforced!'
  write(*,*)

! Dimension of the supermatrix

  n2h1p = nO*nO*nV
  n2p1h = nV*nV*nO
  nH = 1 + n2h1p + n2p1h

! Memory allocation

  allocate(H(nH,nH),eGW(nH),Z(nH))
  
! Initialization

  EcRPA = 0d0

  eF = 0.5d0*(eHF(nO+1) + eHF(nO))

!-------------------!
! Compute screening !
!-------------------!

  ! Spin manifold 
 
  isp_W = 1

  ! Memory allocation

  allocate(Om(nS),Aph(nS,nS),Bph(nS,nS),XpY(nS,nS),XmY(nS,nS),rho(nOrb,nOrb,nS))
 
  call phRLR_A(isp_W,dRPA,nOrb,nC,nO,nV,nR,nS,1d0,eHF,ERI,Aph)
  call phRLR_B(isp_W,dRPA,nOrb,nC,nO,nV,nR,nS,1d0,ERI,Bph)
 
  call phRLR(TDA_W,nS,Aph,Bph,EcRPA,Om,XpY,XmY)

  if(print_W) call print_excitation_energies('phRPA@RHF','singlet',nS,Om)
 
  !--------------------------!
  ! Compute spectral weights !
  !--------------------------!
 
  call RGW_excitation_density(nOrb,nC,nO,nR,nS,ERI,XpY,rho)

  deallocate(Aph,Bph,XpY,XmY)

!-------------------------!
! Main loop over orbitals !
!-------------------------!

  do p=nO,nO

    H(:,:) = 0d0

    !--------------------------------------!
    !     Compute ADC-2SOSEX matrix        !
    !--------------------------------------!
    !                                      !
    !     | F      U_2h1p     U_2p1h     | ! 
    !     |                              | ! 
    ! H = | U_2h1p (K+C)_2h1p 0          | ! 
    !     |                              | ! 
    !     | U_2p1h 0          (K+C)_2p1h | ! 
    !                                      !
    !--------------------------------------!

    call wall_time(start_timing)

    !---------!
    ! Block F !
    !---------!

    H(1,1) = eHF(p)

    !--------------!
    ! Block U_2h1p !
    !--------------!
 
    ija = 0
    do i=nC+1,nO
      do mu=1,nS
        ija = ija + 1
 
        H(1    ,1+ija) = sqrt(2d0)*rho(p,i,mu)
        H(1+ija,1    ) = sqrt(2d0)*rho(p,i,mu)

        do k=nC+1,nO
          do c=nO+1,nOrb-nR
            H(1    ,1+ija) = H(1    ,1+ija) &
                           + sqrt(2d0)*ERI(p,c,k,i)*rho(k,c,mu)/(eHF(c) - eHF(k) - Om(mu)) &
                           + sqrt(2d0)*ERI(p,k,c,i)*rho(c,k,mu)/(eHF(c) - eHF(k) + Om(mu))
            H(1+ija,1    ) = H(1+ija,1    ) & 
                           + sqrt(2d0)*ERI(p,c,k,i)*rho(k,c,mu)/(eHF(c) - eHF(k) - Om(mu)) &
                           + sqrt(2d0)*ERI(p,k,c,i)*rho(c,k,mu)/(eHF(c) - eHF(k) + Om(mu))
          end do
        end do
 
      end do
    end do
 
    !--------------!
    ! Block U_2p1h !
    !--------------!
 
    iab = 0
    do a=nO+1,nOrb-nR
      do mu=1,nS
        iab = iab + 1
 
        H(1          ,1+n2h1p+iab) = sqrt(2d0)*rho(p,a,mu)
        H(1+n2h1p+iab,1          ) = sqrt(2d0)*rho(p,a,mu)
 
        do k=nC+1,nO
          do c=nO+1,nOrb-nR
            H(1    ,1+n2h1p+iab) = H(1    ,1+n2h1p+iab) &
                                 + sqrt(2d0)*ERI(p,k,c,a)*rho(c,k,mu)/(eHF(c) - eHF(k) - Om(mu)) &
                                 + sqrt(2d0)*ERI(p,c,k,a)*rho(k,c,mu)/(eHF(c) - eHF(k) + Om(mu))
            H(1+n2h1p+iab,1    ) = H(1+n2h1p+iab,1    ) &
                                 + sqrt(2d0)*ERI(p,k,c,a)*rho(c,k,mu)/(eHF(c) - eHF(k) - Om(mu)) &
                                 + sqrt(2d0)*ERI(p,c,k,a)*rho(k,c,mu)/(eHF(c) - eHF(k) + Om(mu))
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
 
        H(1+ija,1+ija) = eHF(i) - Om(mu) 

      end do
    end do

    !------------------!
    ! Block (K+C)_2p1h !
    !------------------!
 
    iab = 0
    do a=nO+1,nOrb-nR
      do mu=1,nS
        iab = iab + 1
 
        H(1+n2h1p+iab,1+n2h1p+iab) = eHF(a) + Om(mu)

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

    do s=1,nH
      Z(s) = H(1,s)**2
    end do

  !--------------!
  ! Dump results !
  !--------------!

    write(*,*)'-------------------------------------------'
    write(*,'(1X,A38,I3,A2)')'| ADC-2SOSEX energies (eV) for orbital',p,' |'
    write(*,*)'-------------------------------------------'
    write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X)') &
              '|','#','|','e_QP','|','Z','|'
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
           'Orbital',p,' and #',s,':','e_QP = ',eGW(s)*HaToeV,' eV and Z = ',Z(s)
          write(*,*)'------------------------------------------------------------------------------'
          write(*,'(1X,A20,1X,A20,1X,A15,1X)') &
                    ' Conf. (p,ia)  ',' Coefficient ',' Weight ' 
          write(*,*)'------------------------------------------------------------------------------'
         
          if(p <= nO) & 
            write(*,'(1X,A7,I3,A16,1X,F15.6,1X,F15.6,1X,F12.6)') &
            '      (',p,')               ',H(1,s),H(1,s)**2,-eHF(p)*HaToeV
          if(p > nO) & 
            write(*,'(1X,A16,I3,A7,1X,F15.6,1X,F15.6,1X,F12.6)') &
            '               (',p,')      ',H(1,s),H(1,s)**2,-eHF(p)*HaToeV
  
          ija = 0
          do i=nC+1,nO
            do ja=1,nS
              ija = ija + 1

              if(abs(H(1+ija,s)) > cutoff2)                     &
              write(*,'(1X,A7,I3,A1,I3,A12,1X,F15.6,1X,F15.6,1X,F12.6)') &
              '      (',i,',',ja,')           ',H(1+ija,s),H(1+ija,s)**2,(eHF(i) - Om(ja))*HaToeV
         
            end do
          end do
         
          iab = 0
          do ia=1,nS
            do b=nO+1,nOrb-nR
              iab = iab + 1

                if(abs(H(1+n2h1p+iab,s)) > cutoff2)                 &
                  write(*,'(1X,A7,I3,A1,I3,A12,1X,F15.6,1X,F15.6,1X,F12.6)') &
                  '      (',ia,',',b,')           ',H(1+n2h1p+iab,s),H(1+n2h1p+iab,s)**2,(eHF(b) + Om(ia))*HaToeV
                
            end do
          end do

          write(*,*)'------------------------------------------------------------------------------'
          write(*,*)

        end if ! If state s should be print

      end do ! Loop on s
 
    end if ! If verbose

  end do ! Loop on the orbital in the e block

end subroutine 
