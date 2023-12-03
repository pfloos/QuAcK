subroutine ufG0W0(dotest,TDA_W,nBas,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,eHF)

! Unfold G0W0 equations

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dotest

  logical,intent(in)            :: TDA_W
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: ERHF
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: eHF(nBas)

! Local variables

  integer                       :: p
  integer                       :: s
  integer                       :: i,j,k,l
  integer                       :: a,b,c,d
  integer                       :: jb,kc,ia,ja
  integer                       :: klc,kcd,ija,ijb,iab,jab

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

  logical                       :: verbose = .true.
  double precision,parameter    :: cutoff1 = 0.01d0
  double precision,parameter    :: cutoff2 = 0.01d0
  double precision              :: eF
  double precision,parameter    :: window = 2d0

  double precision              :: start_timing,end_timing,timing

! Output variables

! Hello world

  write(*,*)
  write(*,*)'****************************************'
  write(*,*)'* Restricted Upfolded G0W0 Calculation *'
  write(*,*)'****************************************'
  write(*,*)

! Dimension of the supermatrix

  n2h1p = nO*nO*nV
  n2p1h = nV*nV*nO
  nH = 1 + n2h1p + n2p1h

! Memory allocation

  allocate(H(nH,nH),eGW(nH),Z(nH))
  
! Initialization

  dRPA = .true.
  EcRPA = 0d0

  eF = 0.5d0*(eHF(nO+1) + eHF(nO))

!-------------------!
! Compute screening !
!-------------------!

  if(.not. TDA_W) then

    ! Spin manifold 
 
    ispin = 1

    ! Memory allocation

    allocate(Om(nS),Aph(nS,nS),Bph(nS,nS),XpY(nS,nS),XmY(nS,nS),rho(nBas,nBas,nS))
 
    call phLR_A(ispin,dRPA,nBas,nC,nO,nV,nR,nS,1d0,eHF,ERI,Aph)
    call phLR_B(ispin,dRPA,nBas,nC,nO,nV,nR,nS,1d0,ERI,Bph)
 
    call phLR(TDA_W,nS,Aph,Bph,EcRPA,Om,XpY,XmY)

    if(print_W) call print_excitation_energies('phRPA@RHF','singlet',nS,Om)
 
    !--------------------------!
    ! Compute spectral weights !
    !--------------------------!
 
    call GW_excitation_density(nBas,nC,nO,nR,nS,ERI,XpY,rho)

    deallocate(Aph,Bph,XpY,XmY)

  end if

!-------------------------!
! Main loop over orbitals !
!-------------------------!

  do p=nO-1,nO

    H(:,:) = 0d0

    if(TDA_W) then
 
      ! TDA for W
 
      write(*,*) 'Tamm-Dancoff approximation actived!'
      write(*,*)
 
      !---------------------------!
      !  Compute GW supermatrix   !
      !---------------------------!
      !                           !
      !     |   F   V2h1p V2p1h | ! 
      !     |                   | ! 
      ! H = | V2h1p C2h1p     0 | ! 
      !     |                   | ! 
      !     | V2p1h   0   C2p1h | ! 
      !                           !
      !---------------------------!
 
      call wall_time(start_timing)

      !---------!
      ! Block F !
      !---------!
      
      H(1,1) = eHF(p)
 
      !-------------!
      ! Block V2h1p !
      !-------------!

      ija = 0
      do i=nC+1,nO
        do j=nC+1,nO
          do a=nO+1,nBas-nR
            ija = ija + 1
               
            H(1    ,1+ija) = sqrt(2d0)*ERI(p,a,i,j)
            H(1+ija,1    ) = sqrt(2d0)*ERI(p,a,i,j)
 
          end do
        end do
      end do

      !-------------!
      ! Block V2p1h !
      !-------------!     
 
      iab = 0
      do i=nC+1,nO
        do a=nO+1,nBas-nR
          do b=nO+1,nBas-nR
            iab = iab + 1   
 
            H(1          ,1+n2h1p+iab) = sqrt(2d0)*ERI(p,i,b,a)
            H(1+n2h1p+iab,1          ) = sqrt(2d0)*ERI(p,i,b,a)
               
          end do
        end do
      end do
 
      !-------------!
      ! Block C2h1p !
      !-------------!
 
      ija = 0
      do i=nC+1,nO
        do j=nC+1,nO
          do a=nO+1,nBas-nR
            ija = ija + 1
               
            klc = 0
            do k=nC+1,nO
              do l=nC+1,nO
                do c=nO+1,nBas-nR
                  klc = klc + 1
                     
                  H(1+ija,1+klc) & 
                       = ((eHF(i) + eHF(j) - eHF(a))*Kronecker_delta(j,l)*Kronecker_delta(a,c) & 
                       - 2d0*ERI(j,c,a,l))*Kronecker_delta(i,k)
                     
                end do
              end do
            end do
 
          end do
        end do
      end do
 
      !-------------!
      ! Block C2p1h !
      !-------------!
      
      iab = 0
      do i=nC+1,nO
        do a=nO+1,nBas-nR
          do b=nO+1,nBas-nR
            iab = iab + 1
               
            kcd = 0
            do k=nC+1,nO
              do c=nO+1,nBas-nR
                do d=nO+1,nBas-nR
                  kcd = kcd + 1
                     
                  H(1+n2h1p+iab,1+n2h1p+kcd) &
                       = ((eHF(a) + eHF(b) - eHF(i))*Kronecker_delta(i,k)*Kronecker_delta(a,c) & 
                       + 2d0*ERI(a,k,i,c))*Kronecker_delta(b,d)
                     
                end do
              end do
            end do
        
          end do
        end do
      end do
       
      call wall_time(end_timing)
 
      timing = end_timing - start_timing
      write(*,*)
      write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for construction of supermatrix = ',timing,' seconds'
      write(*,*)

    else
 
      write(*,*) 'Tamm-Dancoff approximation deactivated!'
      write(*,*)
 
      !---------------------------!
      !  Compute GW supermatrix   !
      !---------------------------!
      !                           !
      !     |   F   W2h1p W2p1h | ! 
      !     |                   | ! 
      ! H = | W2h1p D2h1p     0 | ! 
      !     |                   | ! 
      !     | W2p1h   0   D2p1h | ! 
      !                           !
      !---------------------------!

      call wall_time(start_timing)
 
      !---------!
      ! Block F !
      !---------!
 
      H(1,1) = eHF(p)
 
      !-------------!
      ! Block D2h1p !
      !-------------!
 
      ija = 0
      do i=nC+1,nO
        do ja=1,nS
          ija = ija + 1
 
          H(1+ija,1+ija) = eHF(i) - Om(ja) 
 
        end do
      end do
 
      !-------------!
      ! Block W2h1p !
      !-------------!
 
      ija = 0
      do i=nC+1,nO
        do ja=1,nS
          ija = ija + 1
 
          H(1    ,1+ija) = sqrt(2d0)*rho(p,i,ja)
          H(1+ija,1    ) = sqrt(2d0)*rho(p,i,ja)
 
        end do
      end do
 
      !-------------!
      ! Block D2p1h !
      !-------------!
 
      iab = 0
      do ia=1,nS
        do b=nO+1,nBas-nR
          iab = iab + 1
 
          H(1+n2h1p+iab,1+n2h1p+iab) = eHF(b) + Om(ia)
 
        end do
      end do
 
      !-------------!
      ! Block W2p1h !
      !-------------!
 
      iab = 0
      do ia=1,nS
        do b=nO+1,nBas-nR
          iab = iab + 1
 
          H(1          ,1+n2h1p+iab) = sqrt(2d0)*rho(p,b,ia)
          H(1+n2h1p+iab,1          ) = sqrt(2d0)*rho(p,b,ia)
 
        end do
      end do

      call wall_time(end_timing)
 
      timing = end_timing - start_timing
      write(*,*)
      write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for construction of supermatrix = ',timing,' seconds'
      write(*,*)
       
    end if
    
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
    write(*,'(1X,A32,I3,A8)')'| G0W0 energies (eV) for orbital',p,'      |'
    write(*,*)'-------------------------------------------'
    write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X)') &
              '|','#','|','e_QP','|','Z','|'
    write(*,*)'-------------------------------------------'
  
    do s=1,nH
      if(eGW(s) < eF .and. eGW(s) > eF - window) then
!     if(Z(s) > cutoff1) then
        write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
        '|',s,'|',eGW(s)*HaToeV,'|',Z(s),'|'
      end if
    end do
  
    write(*,*)'-------------------------------------------'
    write(*,*)
 
    if(verbose) then 
 
      if(TDA_W) then  
 
        ! TDA printing format
 
        do s=1,nH
       
          if(eGW(s) < eF .and. eGW(s) > eF - window) then
       
            write(*,*)'-------------------------------------------------------------'
            write(*,'(1X,A7,1X,I3,A6,I3,A1,1X,A7,F12.6,A13,F6.4,1X)') & 
             'Orbital',p,' and #',s,':','e_QP = ',eGW(s)*HaToeV,' eV and Z = ',Z(s)
            write(*,*)'-------------------------------------------------------------'
            write(*,'(1X,A20,1X,A20,1X,A15,1X)') &
                      ' Configuration ',' Coefficient ',' Weight ' 
            write(*,*)'-------------------------------------------------------------'
           
            if(p <= nO) & 
              write(*,'(1X,A7,I3,A16,1X,F15.6,1X,F15.6)') &
              '      (',p,')               ',H(1,s),H(1,s)**2
            if(p > nO) & 
              write(*,'(1X,A16,I3,A7,1X,F15.6,1X,F15.6)') &
              '               (',p,')      ',H(1,s),H(1,s)**2
    
            ija = 0
            do i=nC+1,nO
              do j=nC+1,nO
                do a=nO+1,nBas-nR
                  ija = ija + 1
  
                  if(abs(H(1+ija,s)) > cutoff2)               &
                  write(*,'(1X,A3,I3,A1,I3,A6,I3,A7,1X,F15.6,1X,F15.6)') &
                  '  (',i,',',j,') -> (',a,')      ',H(1+ija,s),H(1+ija,s)**2
           
                end do
              end do
            end do
           
            iab = 0
            do i=nC+1,nO
              do a=nO+1,nBas-nR
                do b=nO+1,nBas-nR
                  iab = iab + 1
 
                  if(abs(H(1+n2h1p+iab,s)) > cutoff2)           &
                    write(*,'(1X,A7,I3,A6,I3,A1,I3,A3,1X,F15.6,1X,F15.6)') &
                    '      (',i,') -> (',a,',',b,')  ',H(1+n2h1p+iab,s),H(1+n2h1p+iab,s)**2
                  
                end do
              end do
            end do

            write(*,*)'-------------------------------------------------------------'
            write(*,*)

          end if

        end do
 
      else 
  
        ! non-TDA printing format
 
        do s=1,nH
        
          if(eGW(s) < eF .and. eGW(s) > eF - window) then
        
            write(*,*)'-------------------------------------------------------------'
            write(*,'(1X,A7,1X,I3,A6,I3,A1,1X,A7,F12.6,A13,F6.4,1X)') & 
             'Orbital',p,' and #',s,':','e_QP = ',eGW(s)*HaToeV,' eV and Z = ',Z(s)
            write(*,*)'-------------------------------------------------------------'
            write(*,'(1X,A20,1X,A20,1X,A15,1X)') &
                      ' Conf. (p,ia)  ',' Coefficient ',' Weight ' 
            write(*,*)'-------------------------------------------------------------'
           
            if(p <= nO) & 
              write(*,'(1X,A7,I3,A16,1X,F15.6,1X,F15.6)') &
              '      (',p,')               ',H(1,s),H(1,s)**2
            if(p > nO) & 
              write(*,'(1X,A16,I3,A7,1X,F15.6,1X,F15.6)') &
              '               (',p,')      ',H(1,s),H(1,s)**2
    
            ija = 0
            do i=nC+1,nO
              do ja=1,nS
                ija = ija + 1
  
                if(abs(H(1+ija,s)) > cutoff2)                     &
                write(*,'(1X,A7,I3,A1,I3,A12,1X,F15.6,1X,F15.6)') &
                '      (',i,',',ja,')           ',H(1+ija,s),H(1+ija,s)**2
           
              end do
            end do
           
            iab = 0
            do ia=1,nS
              do b=nO+1,nBas-nR
                iab = iab + 1
 
                  if(abs(H(1+n2h1p+iab,s)) > cutoff2)                 &
                    write(*,'(1X,A7,I3,A1,I3,A12,1X,F15.6,1X,F15.6)') &
                    '      (',ia,',',b,')           ',H(1+n2h1p+iab,s),H(1+n2h1p+iab,s)**2
                  
              end do
            end do
 
            write(*,*)'-------------------------------------------------------------'
            write(*,*)

          end if ! If state s should be print
 
        end do ! Loop on s

      end if ! If TDA
 
    end if ! If verbose

  end do ! Loop on the orbital in the e block

end subroutine 
