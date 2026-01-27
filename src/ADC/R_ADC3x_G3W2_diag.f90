subroutine R_ADC3x_G3W2_diag(dotest,sig_inf,TDA_W,flow,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,eHF)

! ADC(3)-x version of G3W2 within the diagonal approximation

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

  integer                       :: p,r,s
  integer                       :: i,j,k,l
  integer                       :: a,b,c,d
  integer                       :: jb,kc,ia,ja
  integer                       :: mu,nu
  integer                       :: klc,kcd,ija,ijb,iab,jab
  double precision              :: num,num1,num2
  double precision              :: dem,dem1,dem2,dem3
  double precision              :: reg,reg1,reg2,reg3

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

  double precision,allocatable  :: Reigv(:,:) 

  double precision              :: start_timing,end_timing,t_build,t_diag

  logical                       :: add_U1_2h1p
  logical                       :: add_U2_2h1p
  logical                       :: add_U3_2h1p

  logical                       :: add_U1_2p1h
  logical                       :: add_U2_2p1h
  logical                       :: add_U3_2p1h

  logical                       :: add_K_2h1p
  logical                       :: add_K_2p1h

  logical                       :: add_C1_2h1p
  logical                       :: add_C1_2p1h

! Output variables

! Hello world

  write(*,*)
  write(*,*)'***************************************'
  write(*,*)'* Restricted ADC(3x)-G3W2 Calculation *'
  write(*,*)'***************************************'
  write(*,*)

! Diagonal approximation

  write(*,*)' Diagonal approximation enforced!'
  write(*,*)

! Dimension of the supermatrix

  n2h1p = nO*nO*nV
  n2p1h = nV*nV*nO
  nH = 1 + n2h1p + n2p1h

! Select matrix components

! ADC(3x)-G3W2

  add_U3_2h1p = .true.
  add_U3_2p1h = .true.

! ADC(3)-G3W2

  add_C1_2h1p = .true.
  add_C1_2p1h = .true.

! ADC-SOSEX

  add_U2_2h1p = .true.
  add_U2_2p1h = .true.

! ADC-GW

  add_K_2h1p  = .true.
  add_K_2p1h  = .true.

  add_U1_2h1p = .true.
  add_U1_2p1h = .true.

! Memory allocation

  allocate(H(nH,nH),eGW(nH),Z(nH),Reigv(nH,nH))
  
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

  ! Small shift to avoid hard zeros in amplitudes

  Om(:) = Om(:) + 1d-12

  if(print_W) call print_excitation_energies('phRPA@RHF','singlet',nS,Om)
 
  !--------------------------!
  ! Compute spectral weights !
  !--------------------------!
 
  call RGW_excitation_density(nOrb,nC,nO,nR,nS,ERI,XpY,rho)

  deallocate(Aph,Bph,XpY,XmY)

!-------------------------!
! Main loop over orbitals !
!-------------------------!

  do p=nC+1,nO

    H(:,:) = 0d0
 
    !--------------------------------------!
    !  Compute ADC(3x)-G3W2 matrix         !
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

    ! First-order terms

    if(add_U1_2h1p) then

      ija = 0
      do i=nC+1,nO
        do mu=1,nS
          ija = ija + 1
  

          H(1    ,1+ija) = sqrt(2d0)*rho(p,i,mu)

          H(1+ija,1    ) = sqrt(2d0)*rho(p,i,mu)

        end do 
      end do 

    end if

    ! Second-order terms

    if(add_U2_2h1p) then

      ija = 0
      do i=nC+1,nO
        do mu=1,nS
          ija = ija + 1

          do k=nC+1,nO
            do c=nO+1,nOrb-nR

            num = sqrt(2d0)*rho(k,c,mu)*ERI(i,k,c,p)
            dem = eHF(c) - eHF(k) - Om(mu)
            reg = (1d0 - exp(-2d0*flow*dem*dem))/dem

            H(1    ,1+ija) = H(1    ,1+ija) + num*reg
            H(1+ija,1    ) = H(1+ija,1    ) + num*reg

            num = sqrt(2d0)*rho(c,k,mu)*ERI(i,c,k,p)
            dem = eHF(c) - eHF(k) + Om(mu)
            reg = (1d0 - exp(-2d0*flow*dem*dem))/dem

            H(1    ,1+ija) = H(1    ,1+ija) + num*reg
            H(1+ija,1    ) = H(1+ija,1    ) + num*reg

            end do
          end do

        end do
      end do

    end if

    ! Third-order terms

    if(add_U3_2h1p) then

      ija = 0
      do i=nC+1,nO
        do mu=1,nS
          ija = ija + 1

          do k=nC+1,nO
            do c=nO+1,nOrb-nR
              do nu=1,nS

                num = 2d0*sqrt(2d0)*rho(c,k,mu)*rho(i,k,nu)*rho(p,c,nu)
                dem1 = eHF(c) - eHF(k) + Om(mu)
                dem2 = eHF(k) - eHF(i) - Om(nu)

                reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))/dem1
                reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))/dem2

                H(1    ,1+ija) = H(1    ,1+ija) + num*reg1*reg2
                H(1+ija,1    ) = H(1+ija,1    ) + num*reg1*reg2

                num = 2d0*sqrt(2d0)*rho(k,c,mu)*rho(c,i,nu)*rho(k,p,nu)
                dem1 = eHF(k) - eHF(c) + Om(mu)
                dem2 = eHF(i) - eHF(c) - Om(nu)

                reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))/dem1
                reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))/dem2

                H(1    ,1+ija) = H(1    ,1+ija) - num*reg1*reg2
                H(1+ija,1    ) = H(1+ija,1    ) - num*reg1*reg2

                num = 2d0*sqrt(2d0)*rho(k,c,mu)*rho(i,c,nu)*rho(p,k,nu)
                dem1 = eHF(k) - eHF(c) + Om(mu)
                dem2 = eHF(c) - eHF(i) - Om(nu)

                reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))/dem1
                reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))/dem2

                H(1    ,1+ija) = H(1    ,1+ija) - 0.5d0*num*reg1*reg2
                H(1+ija,1    ) = H(1+ija,1    ) - 0.5d0*num*reg1*reg2

                num = 2d0*sqrt(2d0)*rho(k,i,nu)*rho(c,k,mu)*rho(c,p,nu)
                dem1 = eHF(i) - eHF(c) - Om(nu) - Om(mu)
                dem2 = eHF(c) - eHF(k) + Om(mu)

                reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))/dem1
                reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))/dem2

                H(1    ,1+ija) = H(1    ,1+ija) + num*reg1*reg2
                H(1+ija,1    ) = H(1+ija,1    ) + num*reg1*reg2

              end do
            end do
          end do

          do j=nC+1,nO
            do k=nC+1,nO
              do nu=1,nS

                num = 2d0*sqrt(2d0)*rho(k,j,mu)*rho(i,j,nu)*rho(p,k,nu)
                dem1 = eHF(k) - eHF(j) + Om(mu)
                dem2 = eHF(j) - eHF(i) - Om(nu)

                reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))/dem1
                reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))/dem2

                H(1    ,1+ija) = H(1    ,1+ija) + 0.5d0*num*reg1*reg2
                H(1+ija,1    ) = H(1+ija,1    ) + 0.5d0*num*reg1*reg2

              end do
            end do
          end do

          do a=nO+1,nOrb-nR
            do b=nO+1,nOrb-nR
              do nu=1,nS

                num = 2d0*sqrt(2d0)*rho(a,i,nu)*rho(b,a,mu)*rho(b,p,nu)
                dem1 = eHF(i) - eHF(b) - Om(nu) - Om(mu)
                dem2 = eHF(i) - eHF(a) - Om(nu)

                reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))/dem1
                reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))/dem2

                H(1    ,1+ija) = H(1    ,1+ija) + num*reg1*reg2
                H(1+ija,1    ) = H(1+ija,1    ) + num*reg1*reg2

              end do
            end do
          end do
  
        end do
      end do

    end if

    !--------------!
    ! Block U_2p1h !
    !--------------!

    ! First-order terms

    if(add_U1_2p1h) then

      iab = 0
      do a=nO+1,nOrb-nR
        do mu=1,nS
          iab = iab + 1
  

          H(1          ,1+n2h1p+iab) = sqrt(2d0)*rho(a,p,mu)

          H(1+n2h1p+iab,1          ) = sqrt(2d0)*rho(a,p,mu)

        end do
      end do

    end if
  
    ! Second-order terms

    if(add_U2_2p1h) then

      iab = 0
      do a=nO+1,nOrb-nR
        do mu=1,nS
          iab = iab + 1
          do k=nC+1,nO
            do c=nO+1,nOrb-nR

            num = sqrt(2d0)*rho(k,c,mu)*ERI(a,c,k,p)
            dem = eHF(c) - eHF(k) - Om(mu)
            reg = (1d0 - exp(-2d0*flow*dem*dem))/dem

            H(1    ,1+n2h1p+iab) = H(1    ,1+n2h1p+iab) + num*reg
            H(1+n2h1p+iab,1    ) = H(1+n2h1p+iab,1    ) + num*reg

            num = sqrt(2d0)*rho(c,k,mu)*ERI(a,k,c,p)
            dem = eHF(c) - eHF(k) + Om(mu)
            reg = (1d0 - exp(-2d0*flow*dem*dem))/dem

            H(1    ,1+n2h1p+iab) = H(1    ,1+n2h1p+iab) + num*reg
            H(1+n2h1p+iab,1    ) = H(1+n2h1p+iab,1    ) + num*reg

            end do
          end do

        end do
      end do

    end if

    ! Third-order terms

    if(add_U3_2p1h) then

      iab = 0
      do a=nO+1,nOrb-nR
        do mu=1,nS
          iab = iab + 1

          do k=nC+1,nO
            do c=nO+1,nOrb-nR
              do nu=1,nS

                num = 2d0*sqrt(2d0)*rho(c,k,mu)*rho(c,a,nu)*rho(k,p,nu)
                dem1 = eHF(c) - eHF(k) + Om(mu)
                dem2 = eHF(a) - eHF(c) - Om(nu)

                reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))/dem1
                reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))/dem2

                H(1          ,1+n2h1p+iab) = H(1          ,1+n2h1p+iab) + num*reg1*reg2
                H(1+n2h1p+iab,1          ) = H(1+n2h1p+iab,1          ) + num*reg1*reg2

                num = 2d0*sqrt(2d0)*rho(k,c,mu)*rho(a,k,nu)*rho(p,c,nu)
                dem1 = eHF(k) - eHF(c) + Om(mu)
                dem2 = eHF(k) - eHF(a) - Om(nu)

                reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))/dem1
                reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))/dem2

                H(1          ,1+n2h1p+iab) = H(1          ,1+n2h1p+iab) - num*reg1*reg2
                H(1+n2h1p+iab,1          ) = H(1+n2h1p+iab,1          ) - num*reg1*reg2

                num = 2d0*sqrt(2d0)*rho(k,c,mu)*rho(k,a,nu)*rho(c,p,nu)
                dem1 = eHF(k) - eHF(c) + Om(mu)
                dem2 = eHF(a) - eHF(k) - Om(nu)

                reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))/dem1
                reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))/dem2

                H(1          ,1+n2h1p+iab) = H(1          ,1+n2h1p+iab) - 0.5d0*num*reg1*reg2
                H(1+n2h1p+iab,1          ) = H(1+n2h1p+iab,1          ) - 0.5d0*num*reg1*reg2

                num = 2d0*sqrt(2d0)*rho(a,c,nu)*rho(c,k,mu)*rho(p,k,nu)
                dem1 = eHF(a) - eHF(k) + Om(nu) + Om(mu)
                dem2 = eHF(k) - eHF(c) - Om(mu)
                
                reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))/dem1
                reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))/dem2
                
                H(1          ,1+n2h1p+iab) = H(1          ,1+n2h1p+iab) + num*reg1*reg2
                H(1+n2h1p+iab,1          ) = H(1+n2h1p+iab,p          ) + num*reg1*reg2

              end do
            end do
          end do

          do b=nO+1,nOrb-nR
            do c=nO+1,nOrb-nR
              do nu=1,nS

                num = 2d0*sqrt(2d0)*rho(b,c,mu)*rho(b,a,nu)*rho(c,p,nu)
                dem1 = eHF(b) - eHF(c) + Om(mu)
                dem2 = eHF(a) - eHF(b) - Om(nu)

                reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))/dem1
                reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))/dem2

                H(1          ,1+n2h1p+iab) = H(1          ,1+n2h1p+iab) + 0.5d0*num*reg1*reg2
                H(1+n2h1p+iab,1    )       = H(1+n2h1p+iab,1          ) + 0.5d0*num*reg1*reg2

              end do
            end do
          end do

          do i=nC+1,nO
            do j=nC+1,nO
              do nu=1,nS

                num = 2d0*sqrt(2d0)*rho(a,j,nu)*rho(j,i,mu)*rho(p,i,nu)
                dem1 = eHF(a) - eHF(i) + Om(nu) + Om(mu)
                dem2 = eHF(a) - eHF(j) + Om(nu)

                reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))/dem1
                reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))/dem2

                H(1          ,1+n2h1p+iab) = H(1          ,1+n2h1p+iab) + num*reg1*reg2
                H(1+n2h1p+iab,1          ) = H(1+n2h1p+iab,1          ) + num*reg1*reg2

              end do
            end do
          end do

        end do
      end do

    end if

    !-----------------------!
    ! Block (K+C)_2h1p-2h1p !
    !-----------------------!

    ! Zeroth-order terms

    if(add_K_2h1p) then

      ija = 0
      do i=nC+1,nO
        do mu=1,nS
          ija = ija + 1
  
     
          H(1+ija,1+ija) = eHF(i) - Om(mu) 

        end do
      end do

    end if

    ! First-order terms

    if(add_C1_2h1p) then

      ija = 0
      do i=nC+1,nO
        do mu=1,nS
          ija = ija + 1

          klc = 0
          do k=nC+1,nO
            do nu=1,nS
              klc = klc + 1
          
              do r=nC+1,nOrb-nR
       
                num = rho(k,r,mu)*rho(i,r,nu)
                dem = eHF(i) - eHF(r) + Om(nu)
                reg = (1d0 - exp(-2d0*flow*dem*dem))/dem
               
                H(1+ija,1+klc) = H(1+ija,1+klc) + num*reg
               
                num = rho(k,r,mu)*rho(i,r,nu)
                dem = eHF(k) - eHF(r) + Om(mu)
                reg = (1d0 - exp(-2d0*flow*dem*dem))/dem
               
                H(1+ija,1+klc) = H(1+ija,1+klc) + num*reg
       
              end do
   
            end do
          end do

        end do
      end do

    end if

    !-----------------------!
    ! Block (K+C)_2p1h-2p1h !
    !-----------------------!

    ! Zeroth-order terms

    if(add_K_2p1h) then

      iab = 0
      do a=nO+1,nOrb-nR
        do mu=1,nS
          iab = iab + 1
     
     
          H(1+n2h1p+iab,1+n2h1p+iab) = eHF(a) + Om(mu)
     
        end do
      end do

    end if

    ! First-order terms

    if(add_C1_2p1h) then

      iab = 0
      do a=nO+1,nOrb-nR
        do mu=1,nS
          iab = iab + 1

          kcd = 0
          do c=nO+1,nOrb-nR
            do nu=1,nS
              kcd = kcd + 1
         
              do r=nC+1,nOrb-nR

                num = rho(r,c,mu)*rho(r,a,nu)
                dem = eHF(c) - eHF(r) - Om(mu)
                reg = (1d0 - exp(-2d0*flow*dem*dem))/dem
               
                H(1+n2h1p+iab,1+n2h1p+kcd) = H(1+n2h1p+iab,1+n2h1p+kcd) + num*reg
               
                num = rho(r,c,mu)*rho(r,a,nu)
                dem = eHF(a) - eHF(r) - Om(nu)
                reg = (1d0 - exp(-2d0*flow*dem*dem))/dem
               
                H(1+n2h1p+iab,1+n2h1p+kcd) = H(1+n2h1p+iab,1+n2h1p+kcd) + num*reg

              end do
  
            end do
          end do
  
        end do
      end do

    end if

    call wall_time(end_timing)

    t_build = end_timing - start_timing
    write(*,*)
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for construction of supermatrix = ',t_build,' seconds'
    write(*,*)
  
  !-------------------------!
  ! Diagonalize supermatrix !
  !-------------------------!

  call wall_time(start_timing)

  call diagonalize_general_matrix(nH,H,eGW,Reigv)

  call wall_time(end_timing)

  t_diag = end_timing - start_timing
  write(*,*)
  write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for diagonalization of supermatrix = ',t_diag,' seconds'
  write(*,*)

  !-----------------!
  ! Compute weights !
  !-----------------!

    do s=1,nH
      Z(s) = Reigv(1,s)**2
    end do

  !--------------!
  ! Dump results !
  !--------------!

    write(*,*)'-------------------------------------------'
    write(*,'(1X,A35,I3,A5)')'| ADC(3x)-G3W2 energies for orbital',p,'    |'
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
           'Orbital',p,' and #',s,':','e_QP = ',eGW(s)*HaToeV,' eV and Z = ',Z(s)
          write(*,*)'------------------------------------------------------------------------------'
          write(*,'(1X,A20,1X,A20,1X,A15,1X)') &
                    ' Conf. (p,ia)  ',' Coefficient ',' Weight ' 
          write(*,*)'------------------------------------------------------------------------------'
         
          if(p <= nO) & 
            write(*,'(1X,A7,I3,A16,1X,F15.6,1X,F15.6,1X,F12.6)') &
            '      (',p,')               ',Reigv(1,s),Reigv(1,s)**2,-eHF(p)*HaToeV
          if(p > nO) & 
            write(*,'(1X,A16,I3,A7,1X,F15.6,1X,F15.6,1X,F12.6)') &
            '               (',p,')      ',Reigv(1,s),Reigv(1,s)**2,-eHF(p)*HaToeV
  
          ija = 0
          do i=nC+1,nO
            do ja=1,nS
              ija = ija + 1

              if(abs(Reigv(1+ija,s)) > cutoff2)                     &
              write(*,'(1X,A7,I3,A1,I3,A12,1X,F15.6,1X,F15.6,1X,F12.6)') &
              '      (',i,',',ja,')           ',Reigv(1+ija,s),Reigv(1+ija,s)**2,(eHF(i) - Om(ja))*HaToeV
         
            end do
          end do
         
          iab = 0
          do ia=1,nS
            do b=nO+1,nOrb-nR
              iab = iab + 1

                if(abs(Reigv(1+n2h1p+iab,s)) > cutoff2)                 &
                  write(*,'(1X,A7,I3,A1,I3,A12,1X,F15.6,1X,F15.6,1X,F12.6)') &
                  '      (',ia,',',b,')           ',Reigv(1+n2h1p+iab,s),Reigv(1+n2h1p+iab,s)**2,(eHF(b) + Om(ia))*HaToeV
                
            end do
          end do

          write(*,*)'------------------------------------------------------------------------------'
          write(*,*)

        end if ! If state s should be print

      end do ! Loop on s
 
    end if ! If verbose

  end do ! Loop on the orbital in the e block

end subroutine 

