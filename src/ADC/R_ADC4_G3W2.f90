! subroutine R_ADC4_G3W2(dotest,sig_inf,TDA_W,flow,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,eHF)

! ! ADC(4) version of G3W2

!   implicit none
!   include 'parameters.h'

! ! Input variables

!   logical,intent(in)            :: dotest

!   logical,intent(in)            :: sig_inf
!   logical,intent(in)            :: TDA_W
!   double precision,intent(in)   :: flow
!   integer,intent(in)            :: nBas
!   integer,intent(in)            :: nOrb
!   integer,intent(in)            :: nC
!   integer,intent(in)            :: nO
!   integer,intent(in)            :: nV
!   integer,intent(in)            :: nR
!   integer,intent(in)            :: nS
!   double precision,intent(in)   :: ENuc
!   double precision,intent(in)   :: ERHF
!   double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)
!   double precision,intent(in)   :: eHF(nOrb)

! ! Local variables

!   integer                       :: p,q,r,s
!   integer                       :: i,j,k,l
!   integer                       :: a,b,c,d
!   integer                       :: mu,nu
!   integer                       :: klc,kcd,ija,iab
!   double precision              :: num,num1,num2
!   double precision              :: dem1,dem2,dem3
!   double precision              :: reg1,reg2,reg3

!   logical                       :: print_W = .false.
!   logical                       :: dRPA
!   integer                       :: ispin
!   double precision              :: EcRPA
!   integer                       :: n2h1p,n2p1h,nH
!   double precision,external     :: Kronecker_delta
!   double precision,allocatable  :: H(:,:)
!   double precision,allocatable  :: eGW(:)
!   double precision,allocatable  :: Z(:)
!   double precision,allocatable  :: Aph(:,:)
!   double precision,allocatable  :: Bph(:,:)
!   double precision,allocatable  :: Om(:)
!   double precision,allocatable  :: XpY(:,:)
!   double precision,allocatable  :: XmY(:,:)
!   double precision,allocatable  :: rho(:,:,:)

!   logical                       :: verbose = .false.
!   double precision,parameter    :: cutoff1 = 0.1d0
!   double precision,parameter    :: cutoff2 = 0.01d0
!   double precision              :: eF
!   double precision,parameter    :: window = 1.5d0
!   double precision,external     :: SRG_reg

!   double precision              :: start_timing,end_timing,timing

!   double precision,allocatable  :: Reigv(:,:)

!   integer                       :: nIt,maxIt
!   integer,allocatable           :: idx(:)
!   integer,allocatable           :: order(:)
!   double precision,allocatable  :: err(:),eOld(:)
!   double precision              :: w,thresh,Conv

! ! Output variables

! ! Hello world

!   write(*,*)
!   write(*,*)'**************************************'
!   write(*,*)'* Restricted ADC(4)-G3W2 Calculation *'
!   write(*,*)'**************************************'
!   write(*,*)

! ! Dimension of the supermatrix

!   n2h1p = nO*nO*nV
!   n2p1h = nV*nV*nO
!   nH = nOrb + n2h1p + n2p1h

! ! Memory allocation

!   allocate(H(nH,nH),eGW(nH),Z(nH),Reigv(nH,nH))
!   allocate(idx(nOrb),err(nOrb),eOld(nOrb),order(nH))

! ! Initialization

!   dRPA = .true.
!   EcRPA = 0d0

!   eF = 0.5d0*(eHF(nO+1) + eHF(nO))

! !-------------------!
! ! Compute screening !
! !-------------------!

!   ! Memory allocation 

!   allocate(Om(nS),Aph(nS,nS),Bph(nS,nS),XpY(nS,nS),XmY(nS,nS),rho(nOrb,nOrb,nS))

!   ! Spin manifold 

!   ispin = 1

!   call phRLR_A(ispin,dRPA,nOrb,nC,nO,nV,nR,nS,1d0,eHF,ERI,Aph)
!   call phRLR_B(ispin,dRPA,nOrb,nC,nO,nV,nR,nS,1d0,ERI,Bph)

!   call phRLR(TDA_W,nS,Aph,Bph,EcRPA,Om,XpY,XmY)

!   if(print_W) call print_excitation_energies('phRPA@RHF','singlet',nS,Om)

!   !--------------------------!
!   ! Compute spectral weights !
!   !--------------------------!

!   call RGW_excitation_density(nOrb,nC,nO,nR,nS,ERI,XpY,rho)

!   deallocate(Aph,Bph,XpY,XmY)

! ! Initialization

!   H(:,:) = 0d0

!   !--------------------------------------------------------------!
!   !     Compute ADC-G3W2 matrix up to 3h2p/3p2h                  !
!   !--------------------------------------------------------------!
!   !                                                              !
!   !     | F      U_2h1p          U_2p1h          U_3h2p U_3p2h | ! 
!   !     |                                                      | ! 
!   ! H = | U_2h1p (K+C)_2h1p-2h1p C_2p1h-2h1p     0      0      | ! 
!   !     |                                                      | ! 
!   !     | U_2p1h C_2h1p-2p1h     (K+C)_2p1h-2p1h 0      0      | ! 
!   !     |                                                      | ! 
!   !     | U_3h2p 0               0               C_3h2p 0      | !    
!   !     |                                                      | ! 
!   !     | U_3p2h 0               0               0      C_3p2h | !
!   !                                                              !
!   !--------------------------------------------------------------!

!   call wall_time(start_timing)

!   ! Starting loop for non-linear 3h2p/3p2h part

!   thresh = 1d-5
!   Conv   = 1d0
!   nIt    = 0
!   maxIt  = 40

!   w = eHF(nO)

!   write(*,*)'--------------------------------------------------------'
!   write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A10,1X,A1,1X)')'|','It.','|','e_QP (eV)','|','Z','|','Conv.','|'
!   write(*,*)'--------------------------------------------------------'

!   do while(Conv > thresh .and. nIt < maxIt)

!     nIt = nIt + 1

!     !---------!
!     ! Block F !
!     !---------!
 
!     do p=nC+1,nOrb-nR
!       H(p,p) = eHF(p)
!     end do
 
!     do p=nC+1,nOrb-nR
!       do q=nC+1,nOrb-nR
 
!         ! Downfolding the 3h2p configurations
 
!         do i=nC+1,nO
!           do mu=1,nS
!           do nu=1,nS
!               do r=nC+1,nOrb-nR
!               do s=nC+1,nOrb-nR
 
!                 num1 = 2d0*rho(p,r,mu)*rho(r,i,nu)
!                 num2 = 2d0*rho(s,i,mu)*rho(q,s,nu)
!                 dem1 = eHF(i) - eHF(r) - Om(nu)
!                 dem2 = w - eHF(i) + Om(nu) + Om(mu)
!                 dem3 = eHF(i) - eHF(s) - Om(mu)

!                 reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))/dem1
!                 reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))/dem2
!                 reg3 = (1d0 - exp(-2d0*flow*dem3*dem3))/dem3

!                 H(p,q) = H(p,q) + num1*num2*reg1*reg2*reg3
 
!              end do
!              end do
!            end do
!            end do
!          end do
 
!         ! Downfolding the 3p2h configurations
 
!         do a=nO+1,nOrb-nR
!           do mu=1,nS
!           do nu=1,nS
!               do r=nC+1,nOrb-nR
!               do s=nC+1,nOrb-nR
 
!                 num1 = 2d0*rho(r,p,mu)*rho(a,r,nu)
!                 num2 = 2d0*rho(a,s,mu)*rho(s,q,nu)
!                 dem1 = eHF(r) - eHF(a) - Om(nu)
!                 dem2 = w - eHF(a) - Om(nu) - Om(mu)
!                 dem3 = eHF(s) - eHF(a) - Om(mu)
 
!                 reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))/dem1
!                 reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))/dem2
!                 reg3 = (1d0 - exp(-2d0*flow*dem3*dem3))/dem3
 
!                 H(p,q) = H(p,q) + num1*num2*reg1*reg2*reg3
 
!              end do
!              end do
!            end do
!            end do
!          end do
 
!       end do
!     end do
 
!     !--------------!
!     ! Block U_2h1p !
!     !--------------!
 
!     do p=nC+1,nOrb-nR
 
!       ija = 0
!       do i=nC+1,nO
!         do mu=1,nS
!           ija = ija + 1
 
!           H(p       ,nOrb+ija) = sqrt(2d0)*rho(p,i,mu)
!           H(nOrb+ija,p       ) = sqrt(2d0)*rho(p,i,mu)
 
!           do k=nC+1,nO
!             do c=nO+1,nOrb-nR
!               H(p    ,nOrb+ija) = H(p    ,nOrb+ija) &
!                                 + sqrt(2d0)*rho(k,c,mu)*ERI(i,k,c,p)/(eHF(c) - eHF(k) - Om(mu)) &
!                                 + sqrt(2d0)*rho(c,k,mu)*ERI(i,c,k,p)/(eHF(c) - eHF(k) + Om(mu))
!               H(nOrb+ija,p    ) = H(nOrb+ija,p    ) &
!                                 + sqrt(2d0)*rho(k,c,mu)*ERI(i,k,c,p)/(eHF(c) - eHF(k) - Om(mu)) &
!                                 + sqrt(2d0)*rho(c,k,mu)*ERI(i,c,k,p)/(eHF(c) - eHF(k) + Om(mu))
!             end do
!           end do
 
!         end do
!       end do
 
!     end do
 
!     !--------------!
!     ! Block U_2p1h !
!     !--------------!
 
!     do p=nC+1,nOrb-nR
 
!       iab = 0
!       do a=nO+1,nOrb-nR
!         do mu=1,nS
!           iab = iab + 1
 
!           H(p             ,nOrb+n2h1p+iab) = sqrt(2d0)*rho(p,a,mu)
!           H(nOrb+n2h1p+iab,p             ) = sqrt(2d0)*rho(p,a,mu)
 
!           do k=nC+1,nO
!             do c=nO+1,nOrb-nR
!               H(p    ,nOrb+n2h1p+iab) = H(p    ,nOrb+n2h1p+iab) &
!                                       + sqrt(2d0)*rho(k,c,mu)*ERI(a,c,k,p)/(eHF(c) - eHF(k) - Om(mu)) &
!                                       + sqrt(2d0)*rho(c,k,mu)*ERI(a,k,c,p)/(eHF(c) - eHF(k) + Om(mu))
!               H(nOrb+n2h1p+iab,p    ) = H(nOrb+n2h1p+iab,p    ) &
!                                       + sqrt(2d0)*rho(k,c,mu)*ERI(a,c,k,p)/(eHF(c) - eHF(k) - Om(mu)) &
!                                       + sqrt(2d0)*rho(c,k,mu)*ERI(a,k,c,p)/(eHF(c) - eHF(k) + Om(mu))
!             end do
!           end do
 
!         end do
!       end do
 
!     end do
 
!     !------------------!
!     ! Block (K+C)_2h1p !
!     !------------------!
 
!     ija = 0
!     do i=nC+1,nO
!       do mu=1,nS
!         ija = ija + 1
 
!         H(nOrb+ija,nOrb+ija) = eHF(i) - Om(mu) 
 
!        ! First-order terms
 
!         klc = 0
!         do k=nC+1,nO
!           do nu=1,nS
!             klc = klc + 1
 
!             do r=nC+1,nOrb-nR
!               H(nOrb+ija,nOrb+klc) = H(nOrb+ija,nOrb+klc) &
!                                    + 1d0*rho(k,r,mu)*rho(i,r,nu)/(eHF(i) - eHF(r) + Om(nu)) &
!                                    + 1d0*rho(k,r,mu)*rho(i,r,nu)/(eHF(k) - eHF(r) + Om(mu))
!             end do
 
!           end do
!         end do
 
!       end do
!     end do
 
!     !------------------!
!     ! Block (K+C)_2p1h !
!     !------------------!
 
!     iab = 0
!     do a=nO+1,nOrb-nR
!       do mu=1,nS
!         iab = iab + 1
 
!         H(nOrb+n2h1p+iab,nOrb+n2h1p+iab) = eHF(a) + Om(mu)
 
!         kcd = 0
!         do c=nO+1,nOrb-nR
!           do nu=1,nS
!             kcd = kcd + 1
 
!             do r=nC+1,nOrb-nR
!               H(nOrb+n2h1p+iab,nOrb+n2h1p+kcd) = H(nOrb+n2h1p+iab,nOrb+n2h1p+kcd) &
!                                                + 1d0*rho(r,c,mu)*rho(r,a,nu)/(eHF(c) - eHF(r) - Om(mu)) &
!                                                + 1d0*rho(r,c,mu)*rho(r,a,nu)/(eHF(a) - eHF(r) - Om(nu))
!             end do
 
!           end do
!         end do
 
!       end do
!     end do
 
!     !-------------------!
!     ! Block C_2h1p-2p1h !
!     !-------------------!
 
!     ija = 0
!     do i=nC+1,nO
!       do mu=1,nS
!         ija = ija + 1
 
!         kcd = 0
!         do a=nO+1,nOrb-nR
!           do nu=1,nS
!             kcd = kcd + 1
 
!             ! First-order terms
 
!             do k=nC+1,nO
 
!               H(nOrb+ija      ,nOrb+n2h1p+kcd) = H(nOrb+ija      ,nOrb+n2h1p+kcd) &
!                                                + 2d0*rho(k,i,nu)*rho(a,k,mu)/(eHF(a) - eHF(k) + Om(nu))
 
!               H(nOrb+n2h1p+kcd,nOrb+ija      ) = H(nOrb+n2h1p+kcd,nOrb+ija      ) &
!                                                + 2d0*rho(k,i,nu)*rho(a,k,mu)/(eHF(a) - eHF(k) + Om(nu))
 
!             end do
 
!             do c=nO+1,nOrb-nR
 
!               H(nOrb+ija      ,nOrb+n2h1p+kcd) = H(nOrb+ija      ,nOrb+n2h1p+kcd) &
!                                                + 2d0*rho(a,c,nu)*rho(c,i,mu)/(eHF(i) - eHF(c) - Om(mu))
 
!               H(nOrb+n2h1p+kcd,nOrb+ija      ) = H(nOrb+n2h1p+kcd,nOrb+ija      ) &
!                                                + 2d0*rho(a,c,nu)*rho(c,i,mu)/(eHF(i) - eHF(c) - Om(mu))
 
!             end do
 
!           end do
!         end do
 
!       end do
!     end do
 
!     call wall_time(end_timing)

!     timing = end_timing - start_timing
! !   write(*,*)
! !   write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for construction of supermatrix = ',timing,' seconds'
! !   write(*,*)

!   !-------------------------!
!   ! Diagonalize supermatrix !
!   !-------------------------!

!     call wall_time(start_timing)
 
!     call diagonalize_general_matrix(nH,H,eGW,Reigv)
 
!     call wall_time(end_timing)
 
!     timing = end_timing - start_timing
! !   write(*,*)
! !   write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for diagonalization of supermatrix = ',timing,' seconds'
! !   write(*,*)

!   !-----------------!
!   ! Compute weights !
!   !-----------------!

!     Z(:) = 0d0
!     do s=1,nH
!       do p=nC+1,nOrb-nR
!         Z(s) = Z(s) + Reigv(p,s)**2
!       end do
!     end do

!     do s=1,nH
!       order(s) = s
!     end do                      

!     call quick_sort(Z,order,nH)

!     p = 0
!     do s=nH-nOrb+1,nH
!       p = p + 1
!       err(p) = eGW(order(s)) - eOld(p)
!       eOld(p) = eGW(order(s))
!     end do

!     call vecout(nOrb,err)

!   !-----------------------------!
!   ! Update quasiparticle energy !
!   !-----------------------------!

!     Conv = maxval(abs(err))

!     write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F10.6,1X,A1,1X)') '|',nIt,'|',eGW(nO)*HaToeV,'|',Z(nO),'|',Conv,'|'

!   end do

!   write(*,*)'--------------------------------------------------------'
!   write(*,*)

!   if(nIt == maxIt) then

!     write(*,*)
!     write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
!     write(*,*)'                 Convergence failed                 '
!     write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
!     write(*,*)

!   end if

! !--------------!
! ! Dump results !
! !--------------!

!   write(*,*)'---------------------------------------------'
!   write(*,'(1X,A45)')'| ADC(4)-G3W2 energies for all orbitals     |'
!   write(*,*)'---------------------------------------------'
!   write(*,'(1X,A1,1X,A5,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X)') &
!             '|','#','|','e_QP (eV)','|','Z','|'
!   write(*,*)'---------------------------------------------'

!   do s=1,nH
! !   if(eGW(s) < eF .and. eGW(s) > eF - window) then
!     if(Z(s) > cutoff1) then
!       write(*,'(1X,A1,1X,I5,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
!       '|',s,'|',eGW(s)*HaToeV,'|',Z(s),'|'
!     end if
!   end do

!   write(*,*)'---------------------------------------------'
!   write(*,*)

!   if(verbose) then 
  
!     do s=1,nH

!       if(eGW(s) < eF .and. eGW(s) > eF - window) then

!         write(*,*)'-------------------------------------------------------------'
!         write(*,'(1X,A12,1X,I3,A1,1X,A7,F12.6,A13,F6.4,1X)') & 
!          'Eigenvalue #',s,':','e_QP = ',eGW(s)*HaToeV,' eV and Z = ',Z(s)
!         write(*,*)'-------------------------------------------------------------'
!         write(*,'(1X,A20,1X,A20,1X,A15,1X)') &
!                   ' Conf. (p,mu)  ',' Coefficient ',' Weight ' 
!         write(*,*)'-------------------------------------------------------------'
      
!         do p=nC+1,nO 
!           if(abs(Reigv(p,s)) > cutoff2)                     &
!             write(*,'(1X,A7,I3,A16,1X,F15.6,1X,F15.6)') &
!             '      (',p,'    )           ',Reigv(p,s),Reigv(p,s)**2
!         end do
!         do p=nO+1,nOrb-nR
!           if(abs(Reigv(p,s)) > cutoff2)                     &
!             write(*,'(1X,A7,I3,A16,1X,F15.6,1X,F15.6)') &
!             '      (',p,'    )           ',Reigv(p,s),Reigv(p,s)**2
!         end do

!         ija = 0
!         do i=nC+1,nO
!           do mu=1,nS
!             ija = ija + 1
 
!             if(abs(Reigv(nOrb+ija,s)) > cutoff2)                  &
!             write(*,'(1X,A7,I3,A1,I3,A12,1X,F15.6,1X,F15.6)') &
!             '      (',i,',',mu,')           ',Reigv(nOrb+ija,s),Reigv(nOrb+ija,s)**2
       
!           end do
!         end do
       
!         iab = 0
!         do mu=1,nS
!           do b=nO+1,nOrb-nR
!             iab = iab + 1

!               if(abs(Reigv(nOrb+n2h1p+iab,s)) > cutoff2)              &
!                 write(*,'(1X,A7,I3,A1,I3,A12,1X,F15.6,1X,F15.6)') &
!                 '      (',mu,',',b,')           ',Reigv(nOrb+n2h1p+iab,s),Reigv(nOrb+n2h1p+iab,s)**2
              
!           end do
!         end do

!         write(*,*)'-------------------------------------------------------------'
!         write(*,*)

!       end if

!     end do

!   end if

! end subroutine R_ADC4_G3W2

subroutine R_ADC4_G3W2(dotest,sig_inf,TDA_W,flow,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,eHF)

! ADC(4) version of G3W2

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
  integer                       :: klc,kcd,ija,iab,ijkab,ijabc
  double precision              :: num,num1,num2
  double precision              :: dem,dem1,dem2,dem3
  double precision              :: reg,reg1,reg2,reg3

  logical                       :: print_W = .false.
  logical                       :: dRPA
  integer                       :: ispin
  double precision              :: EcRPA
  integer                       :: n2h1p,n2p1h,n3h2p,n3p2h,nH
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
  double precision,parameter    :: cutoff1 = 0.1d0
  double precision,parameter    :: cutoff2 = 0.01d0
  double precision              :: eF
  double precision,parameter    :: window = 1.5d0
  double precision,external     :: SRG_reg

  double precision              :: start_timing,end_timing,timing
  double precision              :: start_time,end_time,time

  logical                       :: add_U1_2h1p
  logical                       :: add_U2_2h1p
  logical                       :: add_U3_2h1p

  logical                       :: add_U1_2p1h
  logical                       :: add_U2_2p1h
  logical                       :: add_U3_2p1h

  logical                       :: add_K_2h1p
  logical                       :: add_K_2p1h

  logical                       :: add_C1_2h1p_2h1p
  logical                       :: add_C1_2p1h_2p1h
  logical                       :: add_C1_2h1p_2p1h
  
  logical                       :: add_K_3h2p
  logical                       :: add_K_3p2h
  
  logical                       :: add_U1_3h2p
  logical                       :: add_U1_3p2h

! Output variables

! Hello world

  write(*,*)
  write(*,*)'**************************************'
  write(*,*)'* Restricted ADC(4)-G3W2 Calculation *'
  write(*,*)'**************************************'
  write(*,*)

! Dimension of the supermatrix

  n2h1p = nO*nO*nV
  n2p1h = nV*nV*nO
  n3h2p = nO*nO*nO*nV*nV
  n3p2h = nV*nV*nV*nO*nO
  nH = nOrb + n2h1p + n2p1h + n3h2p + n3p2h

! Select matrix components

! ADC(4)-G3W2

  add_K_3h2p  = .true.
  add_K_3p2h  = .true.

  add_U1_3h2p = .true.
  add_U1_3p2h = .true.
  
  add_U3_2h1p = .true.
  add_U3_2p1h = .true.
  
! ADC(3)-G3W2

  add_C1_2h1p_2h1p = .true.
  add_C1_2p1h_2p1h = .true.
  add_C1_2h1p_2p1h = .true.

! ADC-SOSEX

  add_U2_2h1p = .true.
  add_U2_2p1h = .true.

  add_K_2h1p  = .true.
  add_K_2p1h  = .true.

! ADC-GW

  add_U1_2h1p = .true.
  add_U1_2p1h = .true.

  
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
  
  ! Small shift to avoid hard zeros in amplitudes

  Om(:) = Om(:) + 1d-12

  if(print_W) call print_excitation_energies('phRPA@RHF','singlet',nS,Om)

!--------------------------!
! Compute spectral weights !
!--------------------------!

  call RGW_excitation_density(nOrb,nC,nO,nR,nS,ERI,XpY,rho)

  deallocate(Aph,Bph,XpY,XmY)

! Initialization

  H(:,:) = 0d0

!--------------------------------------------------------------!
!     Compute ADC-G3W2 matrix up to 3h2p/3p2h                  !
!--------------------------------------------------------------!
!                                                              !
!     | F      U_2h1p          U_2p1h          U_3h2p U_3p2h | ! 
!     |                                                      | ! 
! H = | U_2h1p (K+C)_2h1p-2h1p C_2p1h-2h1p     0      0      | ! 
!     |                                                      | ! 
!     | U_2p1h C_2h1p-2p1h     (K+C)_2p1h-2p1h 0      0      | ! 
!     |                                                      | ! 
!     | U_3h2p 0               0               C_3h2p 0      | !    
!     |                                                      | ! 
!     | U_3p2h 0               0               0      C_3p2h | !
!                                                              !
!--------------------------------------------------------------!

  call wall_time(start_timing)

  write(*,*)'--------------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A10,1X,A1,1X)')'|','It.','|','e_QP (eV)','|','Z','|'
  write(*,*)'--------------------------------------------------------'

  call wall_time(start_time)
  
  !---------!
  ! Block F !
  !---------!
 
  do p=nC+1,nOrb-nR
    H(p,p) = eHF(p)
  end do
 
  !--------------!
  ! Block U_2h1p !
  !--------------!
  
  ! First-order terms

  if(add_U1_2h1p) then

    do p=nC+1,nOrb-nR

      ija = 0
      do i=nC+1,nO
        do mu=1,nS
          ija = ija + 1
 
          H(p       ,nOrb+ija) = sqrt(2d0)*rho(p,i,mu)
          H(nOrb+ija,p       ) = sqrt(2d0)*rho(p,i,mu)
 
        end do
      end do
     
    end do
  
  end if
  
     
  ! Second-order terms

  if(add_U2_2h1p) then

    do p=nC+1,nOrb-nR

      ija = 0
      do i=nC+1,nO
        do mu=1,nS
          ija = ija + 1
 
          do k=nC+1,nO
            do c=nO+1,nOrb-nR
 
              num = sqrt(2d0)*rho(k,c,mu)*ERI(i,k,c,p)
              dem = eHF(c) - eHF(k) - Om(mu)
              reg = (1d0 - exp(-2d0*flow*dem*dem))/dem
 
              H(p       ,nOrb+ija) = H(p       ,nOrb+ija) + num*reg
              H(nOrb+ija,p       ) = H(nOrb+ija,p       ) + num*reg
 
              num = sqrt(2d0)*rho(c,k,mu)*ERI(i,c,k,p)
              dem = eHF(c) - eHF(k) + Om(mu)
              reg = (1d0 - exp(-2d0*flow*dem*dem))/dem
 
              H(p       ,nOrb+ija) = H(p       ,nOrb+ija) + num*reg
              H(nOrb+ija,p       ) = H(nOrb+ija,p       ) + num*reg
 
            end do
          end do
 
        end do
      end do

    end do

  end if
 
  ! Third-order terms

  if(add_U3_2h1p) then
 
     do p=nC+1,nOrb-nR
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
     
                       H(p    ,nOrb+ija) = H(p    ,nOrb+ija) + num*reg1*reg2
                       H(nOrb+ija,p    ) = H(nOrb+ija,p    ) + num*reg1*reg2
     
                       num = 2d0*sqrt(2d0)*rho(k,c,mu)*rho(c,i,nu)*rho(k,p,nu)
                       dem1 = eHF(k) - eHF(c) + Om(mu)
                       dem2 = eHF(i) - eHF(c) - Om(nu)
     
                       reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))/dem1
                       reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))/dem2
     
                       H(p    ,nOrb+ija) = H(p    ,nOrb+ija) - num*reg1*reg2
                       H(nOrb+ija,p    ) = H(nOrb+ija,p    ) - num*reg1*reg2
     
                       num = 2d0*sqrt(2d0)*rho(k,c,mu)*rho(i,c,nu)*rho(p,k,nu)
                       dem1 = eHF(k) - eHF(c) + Om(mu)
                       dem2 = eHF(c) - eHF(i) - Om(nu)
     
                       reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))/dem1
                       reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))/dem2
 
                       H(p    ,nOrb+ija) = H(p    ,nOrb+ija) - 0.5d0*num*reg1*reg2
                       H(nOrb+ija,p    ) = H(nOrb+ija,p    ) - 0.5d0*num*reg1*reg2
 
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
     
                       H(p    ,nOrb+ija) = H(p    ,nOrb+ija) + 0.5d0*num*reg1*reg2
                       H(nOrb+ija,p    ) = H(nOrb+ija,p    ) + 0.5d0*num*reg1*reg2
 
                    end do
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

    do p=nC+1,nOrb-nR

      iab = 0
      do a=nO+1,nOrb-nR
        do mu=1,nS
          iab = iab + 1
 
          H(p             ,nOrb+n2h1p+iab) = sqrt(2d0)*rho(p,a,mu)
          H(nOrb+n2h1p+iab,p             ) = sqrt(2d0)*rho(p,a,mu)
 
        end do
      end do

    end do

  end if
     
  ! Second-order terms

  if(add_U2_2p1h) then

    do p=nC+1,nOrb-nR

      iab = 0
      do a=nO+1,nOrb-nR
        do mu=1,nS
          iab = iab + 1
 
          do k=nC+1,nO
            do c=nO+1,nOrb-nR
 
              num = sqrt(2d0)*rho(k,c,mu)*ERI(a,c,k,p)
              dem = eHF(c) - eHF(k) - Om(mu)
              reg = (1d0 - exp(-2d0*flow*dem*dem))/dem
 
              H(p             ,nOrb+n2h1p+iab) = H(p             ,nOrb+n2h1p+iab) + num*reg
              H(nOrb+n2h1p+iab,p             ) = H(nOrb+n2h1p+iab,p             ) + num*reg
 
              num = sqrt(2d0)*rho(c,k,mu)*ERI(a,k,c,p)
              dem = eHF(c) - eHF(k) + Om(mu)
              reg = (1d0 - exp(-2d0*flow*dem*dem))/dem
 
              H(p             ,nOrb+n2h1p+iab) = H(p             ,nOrb+n2h1p+iab) + num*reg
              H(nOrb+n2h1p+iab,p             ) = H(nOrb+n2h1p+iab,p             ) + num*reg
 
            end do
          end do
 
        end do
      end do

    end do

  end if
 
  ! Third-order terms

  if(add_U3_2p1h) then

     do p=nC+1,nOrb-nR
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
     
                       H(p          ,nOrb+n2h1p+iab) = H(p          ,nOrb+n2h1p+iab) + num*reg1*reg2
                       H(nOrb+n2h1p+iab,p          ) = H(nOrb+n2h1p+iab,p          ) + num*reg1*reg2
     
                       num = 2d0*sqrt(2d0)*rho(k,c,mu)*rho(a,k,nu)*rho(p,c,nu)
                       dem1 = eHF(k) - eHF(c) + Om(mu)
                       dem2 = eHF(k) - eHF(a) - Om(nu)
 
                       reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))/dem1
                       reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))/dem2
     
                       H(p          ,nOrb+n2h1p+iab) = H(p          ,nOrb+n2h1p+iab) - num*reg1*reg2
                       H(nOrb+n2h1p+iab,p          ) = H(nOrb+n2h1p+iab,p          ) - num*reg1*reg2
     
                       num = 2d0*sqrt(2d0)*rho(k,c,mu)*rho(k,a,nu)*rho(c,p,nu)
                       dem1 = eHF(k) - eHF(c) + Om(mu)
                       dem2 = eHF(a) - eHF(k) - Om(nu)
     
                       reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))/dem1
                       reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))/dem2
     
                       H(p          ,nOrb+n2h1p+iab) = H(p          ,nOrb+n2h1p+iab) - 0.5d0*num*reg1*reg2
                       H(nOrb+n2h1p+iab,p          ) = H(nOrb+n2h1p+iab,p          ) - 0.5d0*num*reg1*reg2
     
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
     
                       H(p          ,nOrb+n2h1p+iab) = H(p          ,nOrb+n2h1p+iab) + 0.5d0*num*reg1*reg2
                       H(nOrb+n2h1p+iab,p          ) = H(nOrb+n2h1p+iab,p          ) + 0.5d0*num*reg1*reg2
     
                    end do
                 end do
              end do
     
           end do
        end do
     
     end do
     
  end if
   
  call wall_time(end_time)

  time = end_time - start_time

  write(*,*)
  write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for construction of U 2h1p/2p1h = ',time,' seconds'
  write(*,*)

  call wall_time(start_time)
  
  !------------------!
  ! Block (K+C)_2h1p !
  !------------------!
 
  ! Zeroth-order terms

  if(add_K_2h1p) then

    ija = 0
    do i=nC+1,nO
      do mu=1,nS
        ija = ija + 1
 
        H(nOrb+ija,nOrb+ija) = eHF(i) - Om(mu) 
 
      end do
    end do

  end if

 ! First-order terms

  if(add_C1_2h1p_2h1p) then

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
     
              H(nOrb+ija,nOrb+klc) = H(nOrb+ija,nOrb+klc) + num*reg
     
              num = rho(k,r,mu)*rho(i,r,nu)
              dem = eHF(k) - eHF(r) + Om(mu)
              reg = (1d0 - exp(-2d0*flow*dem*dem))/dem
     
              H(nOrb+ija,nOrb+klc) = H(nOrb+ija,nOrb+klc) + num*reg
     
            end do
     
          end do
        end do

      end do
    end do

  end if

  !------------------!
  ! Block (K+C)_2p1h !
  !------------------!

  ! Zeroth-order terms

  if(add_K_2p1h) then

    iab = 0
    do a=nO+1,nOrb-nR
      do mu=1,nS
        iab = iab + 1
 
        H(nOrb+n2h1p+iab,nOrb+n2h1p+iab) = eHF(a) + Om(mu)

      end do
    end do

  end if

  ! First-order terms

  if(add_C1_2p1h_2p1h) then

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
 
              H(nOrb+n2h1p+iab,nOrb+n2h1p+kcd) = H(nOrb+n2h1p+iab,nOrb+n2h1p+kcd) + num*reg
 
              num = rho(r,c,mu)*rho(r,a,nu)
              dem = eHF(a) - eHF(r) - Om(nu)
              reg = (1d0 - exp(-2d0*flow*dem*dem))/dem
 
              H(nOrb+n2h1p+iab,nOrb+n2h1p+kcd) = H(nOrb+n2h1p+iab,nOrb+n2h1p+kcd) + num*reg
 
            end do
 
          end do
        end do
 
      end do
    end do

  end if

  !-------------------!
  ! Block C_2h1p-2p1h !
  !-------------------!

  ! First-order terms

  if(add_C1_2h1p_2p1h) then

    ija = 0
    do i=nC+1,nO
      do mu=1,nS
        ija = ija + 1
 
        kcd = 0
        do a=nO+1,nOrb-nR
          do nu=1,nS
            kcd = kcd + 1
 
              do k=nC+1,nO
 
                num = 2d0*rho(k,i,mu)*rho(a,k,nu)
                dem = eHF(a) - eHF(k) + Om(nu)
                reg = (1d0 - exp(-2d0*flow*dem*dem))/dem
 
                H(nOrb+ija      ,nOrb+n2h1p+kcd) = H(nOrb+ija      ,nOrb+n2h1p+kcd) + num*reg
 
                H(nOrb+n2h1p+kcd,nOrb+ija      ) = H(nOrb+n2h1p+kcd,nOrb+ija      ) + num*reg
 
              end do
 
              do c=nO+1,nOrb-nR
 
                num = 2d0*rho(c,i,mu)*rho(a,c,nu)
                dem = eHF(i) - eHF(c) - Om(mu)
                reg = (1d0 - exp(-2d0*flow*dem*dem))/dem
 
                H(nOrb+ija      ,nOrb+n2h1p+kcd) = H(nOrb+ija      ,nOrb+n2h1p+kcd) + num*reg
 
                H(nOrb+n2h1p+kcd,nOrb+ija      ) = H(nOrb+n2h1p+kcd,nOrb+ija      ) + num*reg
 
              end do
 
          end do
        end do
 
      end do
    end do

  end if
   
  call wall_time(end_time)

  time = end_time - start_time

  write(*,*)
  write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for construction of C 2h1p/2p1h = ',time,' seconds'
  write(*,*)

  call wall_time(start_time)
  
!------------------!
!   Block U_3h2p   !
!------------------!

  if(add_U1_3h2p) then

     do p=nC+1,nOrb-nR

        ijkab = 0
        do i=nC+1,nO
           do mu=1,nS
              do nu=1,nS
                 ijkab = ijkab + 1  

                 do r=nC+1,nOrb-nR

                    num = 2d0*rho(p,r,mu)*rho(r,i,nu)
                    dem1 = eHF(r) - eHF(i) + Om(nu)
                    reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))/dem1
                    
                    H(p, nOrb + n2h1p + n2p1h + ijkab) = H(p, nOrb + n2h1p + n2p1h + ijkab) + num*reg1
                    H(nOrb + n2h1p + n2p1h + ijkab, p) = H(nOrb + n2h1p + n2p1h + ijkab, p) + num*reg1

                 end do
              
              end do
           end do
        end do
       
    end do
     
  end if
  
!------------------!
!   Block U_3p2h   !
!------------------!
 
  if(add_U1_3p2h) then

     do p=nC+1,nOrb-nR

        ijabc = 0
        do a=nO+1,nOrb-nR
           do mu=1,nS
              do nu=1,nS
                 ijabc = ijabc + 1  

                 do r=nC+1,nOrb-nR

                    num = 2d0*rho(r,p,mu)*rho(a,r,nu)
                    dem1 = eHF(r) - eHF(a) - Om(nu)
                    reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))/dem1
                 
                    H(p, nOrb + n2h1p + n2p1h + n3h2p + ijabc) = H(p, nOrb + n2h1p + n2p1h + n3h2p + ijabc) + num*reg1
                    H(nOrb + n2h1p + n2p1h + n3h2p + ijabc, p) = H(nOrb + n2h1p + n2p1h + n3h2p + ijabc, p) + num*reg1

                 end do
                 
              end do
           end do
        end do

    end do
     
  end if
   
  call wall_time(end_time)

  time = end_time - start_time

  write(*,*)
  write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for construction of U 3h2p/3p2h = ',time,' seconds'
  write(*,*)

  call wall_time(start_time)
  
!------------------!
!   Block K_3h2p   !
!------------------!

  if(add_K_3h2p) then
     ijkab = 0
     do i=nC+1,nO
        do mu=1,nS
           do nu=1,nS
              ijkab = ijkab + 1  
              
              H(nOrb + n2h1p + n2p1h + ijkab, nOrb + n2h1p + n2p1h + ijkab) = eHF(i) - Om(mu) - Om(nu)
              
           end do
        end do
     end do
  end if
  
  
!------------------!
!   Block K_3p2h   !
!------------------!

  if(add_K_3h2p) then
     ijabc = 0
     do a=nO+1,nOrb-nR
        do mu=1,nS
           do nu=1,nS
              ijabc = ijabc + 1  
  
              H(nOrb + n2h1p + n2p1h + n3h2p + ijabc, nOrb + n2h1p + n2p1h + n3h2p + ijabc) = eHF(a) + Om(mu) + Om(nu)
              
           end do
        end do
     end do
  end if
   
  call wall_time(end_time)

  time = end_time - start_time

  write(*,*)
  write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for construction of K 3h2p/3p2h = ',time,' seconds'
  write(*,*)
 
  call wall_time(end_timing)

  timing = end_timing - start_timing
! write(*,*)
! write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for construction of supermatrix = ',timing,' seconds'
! write(*,*)

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
    do p=nC+1,nOrb-nR
      Z(s) = Z(s) + H(p,s)**2
    end do
  end do          

  write(*,*)'--------------------------------------------------------'
  write(*,*)

!--------------!
! Dump results !
!--------------!

  write(*,*)'---------------------------------------------'
  write(*,'(1X,A45)')'| ADC(4)-G3W2 energies for all orbitals     |'
  write(*,*)'---------------------------------------------'
  write(*,'(1X,A1,1X,A5,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X)') &
            '|','#','|','e_QP (eV)','|','Z','|'
  write(*,*)'---------------------------------------------'

  do s=1,nH
    if(eGW(s) < eF .and. eGW(s) > eF - window) then
    ! if(Z(s) > cutoff1) then
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
        do p=nO+1,nOrb-nR
          if(abs(H(p,s)) > cutoff2)                     &
            write(*,'(1X,A7,I3,A16,1X,F15.6,1X,F15.6)') &
            '      (',p,'    )           ',H(p,s),H(p,s)**2
        end do

        ija = 0
        do i=nC+1,nO
          do mu=1,nS
            ija = ija + 1
 
            if(abs(H(nOrb+ija,s)) > cutoff2)                  &
            write(*,'(1X,A7,I3,A1,I3,A12,1X,F15.6,1X,F15.6)') &
            '      (',i,',',mu,')           ',H(nOrb+ija,s),H(nOrb+ija,s)**2
       
          end do
        end do
       
        iab = 0
        do mu=1,nS
          do b=nO+1,nOrb-nR
            iab = iab + 1

              if(abs(H(nOrb+n2h1p+iab,s)) > cutoff2)              &
                write(*,'(1X,A7,I3,A1,I3,A12,1X,F15.6,1X,F15.6)') &
                '      (',mu,',',b,')           ',H(nOrb+n2h1p+iab,s),H(nOrb+n2h1p+iab,s)**2
              
          end do
        end do

        write(*,*)'-------------------------------------------------------------'
        write(*,*)

      end if

    end do

  end if

end subroutine R_ADC4_G3W2
