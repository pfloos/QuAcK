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
  integer                       :: mu,nu,la
  integer                       :: bra,ket
  double precision              :: num,num1,num2
  double precision              :: dem,dem1,dem2
  double precision              :: reg,reg1,reg2

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
  double precision,external     :: SRG_reg2

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

  logical                       :: add_C1_2h1p
  logical                       :: add_C1_2p1h

  logical                       :: add_C1_2h1p_3h2p
  logical                       :: add_C1_2p1h_3p2h
  
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

  add_C1_2h1p_3h2p  = .true.
  add_C1_2p1h_3p2h  = .true.

  add_U1_3h2p = .true.
  add_U1_3p2h = .true.

  add_K_3h2p  = .true.
  add_K_3p2h  = .true.

! ADC(3x)-G3W2

  add_U3_2h1p = .true.
  add_U3_2p1h = .true.
  
! ADC(3)-G3W2

  add_C1_2h1p = .true.
  add_C1_2p1h = .true.

! ADC(2x)-G3W2

  add_U2_2h1p = .true.
  add_U2_2p1h = .true.

  add_K_2h1p  = .true.
  add_K_2p1h  = .true.

! ADC(2)-G3W2

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

!------------------------------------------------------------------!
!     Compute ADC-G3W2 matrix up to 3h2p/3p2h                      !
!------------------------------------------------------------------!
!                                                                  !
!     | F      U_2h1p      U_2p1h      U_3h2p        U_3p2h      | ! 
!     |                                                          | ! 
! H = | U_2h1p (K+C)_2h1p  0           C_2h1p-3h2p   0           | ! 
!     |                                                          | ! 
!     | U_2p1h 0           (K+C)_2p1h  0             C_2h1p-3h2p | !
!     |                                                          | ! 
!     | U_3h2p C_2h1p-3h2p 0           C_3h2p 0                  | !    
!     |                                                          | ! 
!     | U_3p2h 0           C_2p1h-3p2h 0      C_3p2h 0           | !
!                                                                  !
!------------------------------------------------------------------!

  call wall_time(start_timing)

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

  call wall_time(start_time)

  if(add_U1_2h1p) then

    do p=nC+1,nOrb-nR

      bra = nOrb
      do i=nC+1,nO
        do mu=1,nS
          bra = bra + 1
 
          H(p  ,bra) = sqrt(2d0)*rho(p,i,mu)
          H(bra,p  ) = sqrt(2d0)*rho(p,i,mu)
 
        end do
      end do
     
    end do
  
  end if
  
  call wall_time(end_time)

  time = end_time - start_time

  write(*,*)
  write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for construction of U1_2h1p = ',time,' seconds'
  write(*,*)
     
  ! Second-order terms

  call wall_time(start_time)

  if(add_U2_2h1p) then

    do p=nC+1,nOrb-nR

      bra = nOrb
      do i=nC+1,nO
        do mu=1,nS
          bra = bra + 1
 
          do k=nC+1,nO
            do c=nO+1,nOrb-nR
         
              num = sqrt(2d0)*ERI(i,c,k,p)*rho(k,c,mu)
              dem = eHF(c) - eHF(k) + Om(mu)
              reg = (1d0 - exp(-2d0*flow*dem*dem))/dem
         
              H(p  ,bra) = H(p  ,bra) + num*reg
              H(bra,p  ) = H(bra,p  ) + num*reg
         
              num = sqrt(2d0)*ERI(i,k,c,p)*rho(k,c,mu)
              dem = eHF(c) - eHF(k) - Om(mu)
              reg = (1d0 - exp(-2d0*flow*dem*dem))/dem
         
              H(p  ,bra) = H(p  ,bra) + num*reg
              H(bra,p  ) = H(bra,p  ) + num*reg
         
            end do
          end do

        end do
      end do

    end do

  end if

  call wall_time(end_time)

  time = end_time - start_time

  write(*,*)
  write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for construction of U2_2h1p = ',time,' seconds'
  write(*,*)

  ! Third-order terms

  call wall_time(start_time)

  if(add_U3_2h1p) then
 
    do p=nC+1,nOrb-nR
      bra = nOrb
      do i=nC+1,nO
        do mu=1,nS
          bra = bra + 1
    
          do k=nC+1,nO
            do c=nO+1,nOrb-nR
              do nu=1,nS
    
                num = 2d0*sqrt(2d0)*rho(k,c,mu)*rho(c,i,nu)*rho(k,p,nu)
                dem1 = eHF(c) - eHF(k) - Om(mu)
                dem2 = eHF(c) - eHF(i) + Om(nu)
    
                reg = SRG_reg2(dem1,dem2,flow)
    
                H(p  ,bra) = H(p  ,bra) - num*reg
                H(bra,p  ) = H(bra,p  ) - num*reg
    
                num = 2d0*sqrt(2d0)*rho(k,c,mu)*rho(i,c,nu)*rho(p,k,nu)
                dem1 = eHF(c) - eHF(k) - Om(mu)
                dem2 = eHF(c) - eHF(i) - Om(nu)
    
                reg = SRG_reg2(dem1,dem2,flow)
    
                H(p  ,bra) = H(p  ,bra) + 0.5d0*num*reg
                H(bra,p  ) = H(bra,p  ) + 0.5d0*num*reg
    
                num = 2d0*sqrt(2d0)*rho(k,i,nu)*rho(c,k,mu)*rho(c,p,nu)
                dem1 = eHF(c) - eHF(i) + Om(mu) + Om(nu)
                dem2 = eHF(c) - eHF(k) + Om(mu)
    
                reg = SRG_reg2(dem1,dem2,flow)
    
                H(p  ,bra) = H(p  ,bra) - num*reg
                H(bra,p  ) = H(bra,p  ) - num*reg
    
              end do
            end do
          end do
    
          do a=nO+1,nOrb-nR
            do b=nO+1,nOrb-nR
              do nu=1,nS
    
                num = 2d0*sqrt(2d0)*rho(a,i,nu)*rho(b,a,mu)*rho(b,p,nu)
                dem1 = eHF(b) - eHF(i) + Om(mu) + Om(nu)
                dem2 = eHF(a) - eHF(i) + Om(nu)

                reg = SRG_reg2(dem1,dem2,flow)

                H(p  ,bra) = H(p  ,bra) + num*reg
                H(bra,p  ) = H(bra,p  ) + num*reg
 
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
  write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for construction of U3_2h1p = ',time,' seconds'
  write(*,*)
     
  !--------------!
  ! Block U_2p1h !
  !--------------!
  
  ! First-order terms

  call wall_time(start_time)

  if(add_U1_2p1h) then

    do p=nC+1,nOrb-nR

      bra = nOrb + n2h1p
      do a=nO+1,nOrb-nR
        do mu=1,nS
          bra = bra + 1
 
          H(p  ,bra) = sqrt(2d0)*rho(p,a,mu)
          H(bra,p  ) = sqrt(2d0)*rho(p,a,mu)
 
        end do
      end do

    end do

  end if
   
  call wall_time(end_time)

  time = end_time - start_time

  write(*,*)
  write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for construction of U1_2p1h = ',time,' seconds'
  write(*,*)  

  ! Second-order terms

  call wall_time(start_time)

  if(add_U2_2p1h) then

    do p=nC+1,nOrb-nR

      bra = nOrb + n2h1p
      do a=nO+1,nOrb-nR
        do mu=1,nS
          bra = bra + 1
 
          do k=nC+1,nO
            do c=nO+1,nOrb-nR

              num = sqrt(2d0)*ERI(a,k,c,p)*rho(k,c,mu)
              dem = eHF(c) - eHF(k) + Om(mu)
              reg = (1d0 - exp(-2d0*flow*dem*dem))/dem

              H(p  ,bra) = H(p  ,bra) + num*reg
              H(bra,p  ) = H(bra,p  ) + num*reg

              num = sqrt(2d0)*ERI(a,c,k,p)*rho(k,c,mu)
              dem = eHF(c) - eHF(k) - Om(mu)
              reg = (1d0 - exp(-2d0*flow*dem*dem))/dem

              H(p  ,bra) = H(p  ,bra) + num*reg
              H(bra,p  ) = H(bra,p  ) + num*reg

            end do
          end do
 
        end do
      end do

    end do

  end if

  call wall_time(end_time)

  time = end_time - start_time

  write(*,*)
  write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for construction of U2_2p1h = ',time,' seconds'
  write(*,*)
 
  ! Third-order terms

  call wall_time(start_time)

  if(add_U3_2p1h) then

    do p=nC+1,nOrb-nR
      bra = nOrb + n2h1p
      do a=nO+1,nOrb-nR
        do mu=1,nS
          bra = bra + 1
    
          do k=nC+1,nO
            do c=nO+1,nOrb-nR
              do nu=1,nS
    
                num = 2d0*sqrt(2d0)*rho(k,c,mu)*rho(a,k,nu)*rho(p,c,nu)
                dem1 = eHF(c) - eHF(k) - Om(mu)
                dem2 = eHF(a) - eHF(k) + Om(nu)

                reg = SRG_reg2(dem1,dem2,flow)
 
                H(p  ,bra) = H(p  ,bra) - num*reg
                H(bra,p  ) = H(bra,p  ) - num*reg
    
                num = 2d0*sqrt(2d0)*rho(k,c,mu)*rho(k,a,nu)*rho(c,p,nu)
                dem1 = eHF(c) - eHF(k) - Om(mu)
                dem2 = eHF(a) - eHF(k) - Om(nu)

                reg = SRG_reg2(dem1,dem2,flow)
    
                H(p  ,bra) = H(p  ,bra) + 0.5d0*num*reg
                H(bra,p  ) = H(bra,p  ) + 0.5d0*num*reg
    
                num = 2d0*sqrt(2d0)*rho(a,c,nu)*rho(c,k,mu)*rho(p,k,nu)
                dem1 = eHF(a) - eHF(k) + Om(mu) + Om(nu)
                dem2 = eHF(c) - eHF(k) + Om(mu)
    
                reg = SRG_reg2(dem1,dem2,flow)

                H(p  ,bra) = H(p  ,bra) - num*reg
                H(bra,p  ) = H(bra,p  ) - num*reg
    
              end do
            end do
          end do
    
          do i=nC+1,nO
            do j=nC+1,nO
              do nu=1,nS
    
                num = 2d0*sqrt(2d0)*rho(a,j,nu)*rho(j,i,mu)*rho(p,i,nu)
                dem1 = eHF(a) - eHF(i) + Om(mu) + Om(nu)
                dem2 = eHF(a) - eHF(j) + Om(nu)

                reg = SRG_reg2(dem1,dem2,flow)
                
                H(p  ,bra) = H(p  ,bra) + num*reg
                H(bra,p  ) = H(bra,p  ) + num*reg
    
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
  write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for construction of U3_2p1h = ',time,' seconds'
  write(*,*)
  
  !------------------!
  ! Block (K+C)_2h1p !
  !------------------!
 
  ! Zeroth-order terms

  if(add_K_2h1p) then

    bra = nOrb
    do i=nC+1,nO
      do mu=1,nS
        bra = bra + 1
 
        H(bra,bra) = eHF(i) - Om(mu) 
 
      end do
    end do

  end if

 ! First-order terms

  call wall_time(start_time)

  if(add_C1_2h1p) then

    bra = nOrb
    do i=nC+1,nO
      do mu=1,nS
        bra = bra + 1

        ket = nOrb
        do k=nC+1,nO
          do nu=1,nS
            ket = ket + 1
 
            do j=nC+1,nO

              num = rho(k,j,mu)*rho(i,j,nu)
              dem = eHF(i) - eHF(j) + Om(nu)
              reg = (1d0 - exp(-2d0*flow*dem*dem))/dem

              H(bra,ket) = H(bra,ket) + num*reg

              num = rho(k,j,mu)*rho(i,j,nu)
              dem = eHF(k) - eHF(j) + Om(mu)
              reg = (1d0 - exp(-2d0*flow*dem*dem))/dem

              H(bra,ket) = H(bra,ket) + num*reg

            end do
    
          end do
        end do

      end do
    end do

  end if

  call wall_time(end_time)

  time = end_time - start_time

  write(*,*)
  write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for construction of C1_2h1p = ',time,' seconds'
  write(*,*)

  !------------------!
  ! Block (K+C)_2p1h !
  !------------------!

  ! Zeroth-order terms

  if(add_K_2p1h) then

    bra = nOrb + n2h1p
    do a=nO+1,nOrb-nR
      do mu=1,nS
        bra = bra + 1
 
        H(bra,bra) = eHF(a) + Om(mu)

      end do
    end do

  end if

  ! First-order terms

  call wall_time(start_time)

  if(add_C1_2p1h) then

    bra = nOrb + n2h1p
    do a=nO+1,nOrb-nR
      do mu=1,nS
        bra = bra + 1
 
        ket = nOrb + n2h1p
        do c=nO+1,nOrb-nR
          do nu=1,nS
            ket = ket + 1
 

            do b=nO+1,nOrb-nR

              num = rho(b,c,mu)*rho(b,a,nu)
              dem = eHF(c) - eHF(b) - Om(mu)
              reg = (1d0 - exp(-2d0*flow*dem*dem))/dem

              H(bra,ket) = H(bra,ket) + num*reg

              num = rho(b,c,mu)*rho(b,a,nu)
              dem = eHF(a) - eHF(b) - Om(nu)
              reg = (1d0 - exp(-2d0*flow*dem*dem))/dem

              H(bra,ket) = H(bra,ket) + num*reg

            end do
 
          end do
        end do
 
      end do
    end do

  end if

  call wall_time(end_time)

  time = end_time - start_time

  write(*,*)
  write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for construction of C1_2p1h = ',time,' seconds'
  write(*,*)

!------------------!
!   Block U1_3h2p  !
!------------------!

  call wall_time(start_time)

  if(add_U1_3h2p) then

    do p=nC+1,nOrb-nR

      bra = nOrb + n2h1p + n2p1h
      do i=nC+1,nO
        do mu=1,nS
          do nu=1,nS
            bra = bra + 1  

            do a=nO+1,nOrb-nR

              num = 2d0*rho(a,i,nu)*rho(p,a,mu)
              dem1 = eHF(a) - eHF(i) + Om(nu)
              reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))/dem1
              
              H(p  ,bra) = H(p  ,bra) - num*reg1
              H(bra,p  ) = H(bra,p  ) - num*reg1

             end do
           
           end do
         end do
       end do
       
    end do
     
  end if
  
  call wall_time(end_time)

  time = end_time - start_time

  write(*,*)
  write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for construction of U1_3h2p = ',time,' seconds'
  write(*,*)

!------------------!
!   Block U1_3p2h  !
!------------------!
 
  call wall_time(start_time)

  if(add_U1_3p2h) then

    do p=nC+1,nOrb-nR

      bra = nOrb + n2h1p + n2p1h + n3h2p
      do a=nO+1,nOrb-nR
        do mu=1,nS
          do nu=1,nS
            bra = bra + 1  

            do i=nC+1,nO

              num = 2d0*rho(a,i,nu)*rho(i,p,mu)
              dem1 = eHF(a) - eHF(i) + Om(nu)
              reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))/dem1
              
              H(p  ,bra) = H(p  ,bra) + num*reg1
              H(bra,p  ) = H(bra,p  ) + num*reg1

            end do
              
          end do
        end do
      end do

    end do
     
  end if
   
  call wall_time(end_time)

  time = end_time - start_time

  write(*,*)
  write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for construction of U1_3p2h = ',time,' seconds'
  write(*,*)
  
!------------------!
!   Block K_3h2p   !
!------------------!

  call wall_time(start_time)

  if(add_K_3h2p) then

    bra = nOrb + n2h1p + n2p1h
    do i=nC+1,nO
      do mu=1,nS
        do nu=1,nS
          bra = bra + 1  
             
          H(bra,bra) = eHF(i) - Om(mu) - Om(nu)
             
        end do
      end do
    end do

  end if
  
   
  call wall_time(end_time)

  time = end_time - start_time

  write(*,*)
  write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for construction of K_3h2p = ',time,' seconds'
  write(*,*)

!------------------!
!   Block K_3p2h   !
!------------------!

  call wall_time(start_time)

  if(add_K_3h2p) then

    bra = nOrb + n2h1p + n2p1h + n3h2p
    do a=nO+1,nOrb-nR
      do mu=1,nS
        do nu=1,nS
          bra = bra + 1  
  
          H(bra,bra) = eHF(a) + Om(mu) + Om(nu)
             
        end do
      end do
    end do

  end if
   
  call wall_time(end_time)

  time = end_time - start_time

  write(*,*)
  write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for construction of K_3p2h = ',time,' seconds'
  write(*,*)
 

!--------------------!
! Block C1_2h1p-3h2p !
!--------------------!

  call wall_time(start_time)

  if(add_C1_2h1p_3h2p) then

    bra = nOrb 
    do i=nC+1,nO
      do mu=1,nS
        bra = bra + 1

        ket = nOrb + n2h1p + n2p1h
        do k=nC+1,nO
          do nu=1,nS
            do la=1,nS
              ket = ket + 1

              H(bra,ket) = sqrt(2d0)*rho(i,k,nu)*Kronecker_delta(mu,la)
              H(ket,bra) = sqrt(2d0)*rho(i,k,nu)*Kronecker_delta(mu,la)

            end do
          end do
        end do

      end do
    end do

  end if

  call wall_time(end_time)

  time = end_time - start_time

  write(*,*)
  write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for construction of C1_2h1p-3h2p = ',time,' seconds'
  write(*,*)


!--------------------!
! Block C1_2p1h-3p2h !
!--------------------!

  call wall_time(start_time)

  if(add_C1_2p1h_3p2h) then

    bra = nOrb + n2h1p
    do a=nO+1,nOrb-nR
      do mu=1,nS
        bra = bra + 1

        ket = nOrb + n2h1p + n2p1h + n3h2p
        do c=nO+1,nOrb-nR
          do nu=1,nS
            do la=1,nS
              ket = ket + 1

              H(bra,ket) = sqrt(2d0)*rho(c,a,nu)*Kronecker_delta(mu,la)
              H(ket,bra) = sqrt(2d0)*rho(c,a,nu)*Kronecker_delta(mu,la)

            end do
          end do
        end do

      end do
    end do

  end if

  call wall_time(end_time)

  time = end_time - start_time

  write(*,*)
  write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for construction of C1_2p1h-3p2h = ',time,' seconds'
  write(*,*)

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
    do p=nC+1,nOrb-nR
      Z(s) = Z(s) + H(p,s)**2
    end do
  end do          

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
        do p=nO+1,nOrb-nR
          if(abs(H(p,s)) > cutoff2)                     &
            write(*,'(1X,A7,I3,A16,1X,F15.6,1X,F15.6)') &
            '      (',p,'    )           ',H(p,s),H(p,s)**2
        end do

        bra = 0
        do i=nC+1,nO
          do mu=1,nS
            bra = bra + 1
 
            if(abs(H(nOrb+bra,s)) > cutoff2)                  &
            write(*,'(1X,A7,I3,A1,I3,A12,1X,F15.6,1X,F15.6)') &
            '      (',i,',',mu,')           ',H(nOrb+bra,s),H(nOrb+bra,s)**2
       
          end do
        end do
       
        bra = 0
        do mu=1,nS
          do b=nO+1,nOrb-nR
            bra = bra + 1

              if(abs(H(nOrb+n2h1p+bra,s)) > cutoff2)              &
                write(*,'(1X,A7,I3,A1,I3,A12,1X,F15.6,1X,F15.6)') &
                '      (',mu,',',b,')           ',H(nOrb+n2h1p+bra,s),H(nOrb+n2h1p+bra,s)**2
              
          end do
        end do

        write(*,*)'-------------------------------------------------------------'
        write(*,*)

      end if

    end do

  end if

end subroutine 
