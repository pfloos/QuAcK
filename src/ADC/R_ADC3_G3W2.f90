subroutine R_ADC3_G3W2(dotest,TDA_W,flow,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,eHF)

! ADC(3) version of G3W2

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dotest

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

  integer                       :: p,r
  integer                       :: s
  integer                       :: i,j,k,l
  integer                       :: a,b,c,d
  integer                       :: mu,nu
  integer                       :: klc,kcd,ija,iab
  double precision              :: num,dem,reg

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

  logical                       :: verbose = .false.
  double precision,parameter    :: cutoff1 = 0.1d0
  double precision,parameter    :: cutoff2 = 0.01d0
  double precision              :: eF
  double precision,parameter    :: window = 1.5d0
  double precision,external     :: SRG_reg

  double precision              :: start_timing,end_timing,timing

! Output variables

! Hello world

  write(*,*)
  write(*,*)'**************************************'
  write(*,*)'* Restricted ADC(3)-G3W2 Calculation *'
  write(*,*)'**************************************'
  write(*,*)

! Dimension of the supermatrix

  n2h1p = nO*nO*nV
  n2p1h = nV*nV*nO
  nH = nOrb + n2h1p + n2p1h

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

! Initialization

  H(:,:) = 0d0

  !-------------------------------------------------!
  !     Compute ADC-G3W2 matrix up to 2h1p/2p1h     !
  !-------------------------------------------------!
  !                                                 !
  !     | F      U_2h1p          U_2p1h           | ! 
  !     |                                         | ! 
  ! H = | U_2h1p (K+C)_2h1p-2h1p C_2p1h-2h1p      | ! 
  !     |                                         | ! 
  !     | U_2p1  C_2h1p-2p1h     (K+C)_2p1h-2p1h  | ! 
  !                                                 !
  !-------------------------------------------------!

  call wall_time(start_timing)

  !---------!
  ! Block F !
  !---------!

  do p=nC+1,nOrb-nR
    H(p,p) = eHF(p)
  end do

  !--------------!
  ! Block U_2h1p !
  !--------------!

  do p=nC+1,nOrb-nR

    ija = 0
    do i=nC+1,nO
      do mu=1,nS
        ija = ija + 1

        H(p       ,nOrb+ija) = sqrt(2d0)*rho(p,i,mu)
        H(nOrb+ija,p       ) = sqrt(2d0)*rho(p,i,mu)

        do k=nC+1,nO
          do c=nO+1,nOrb-nR

            num = sqrt(2d0)*ERI(p,c,k,i)*rho(k,c,mu)
            dem = eHF(c) - eHF(k) - Om(mu)
            reg = SRG_reg(flow,dem)

            H(p       ,nOrb+ija) = H(p       ,nOrb+ija) + num*reg
            H(nOrb+ija,p       ) = H(nOrb+ija,p       ) + num*reg

            num = sqrt(2d0)*ERI(p,k,c,i)*rho(c,k,mu)
            dem = eHF(c) - eHF(k) + Om(mu)
            reg = SRG_reg(flow,dem)

            H(p       ,nOrb+ija) = H(p       ,nOrb+ija) + num*reg
            H(nOrb+ija,p       ) = H(nOrb+ija,p       ) + num*reg

          end do
        end do

      end do
    end do

  end do

  !--------------!
  ! Block U_2p1h !
  !--------------!

  do p=nC+1,nOrb-nR

    iab = 0
    do a=nO+1,nOrb-nR
      do mu=1,nS
        iab = iab + 1

        H(p             ,nOrb+n2h1p+iab) = sqrt(2d0)*rho(p,a,mu)
        H(nOrb+n2h1p+iab,p             ) = sqrt(2d0)*rho(p,a,mu)

        do k=nC+1,nO
          do c=nO+1,nOrb-nR

            num = sqrt(2d0)*ERI(p,k,c,a)*rho(c,k,mu)
            dem = eHF(c) - eHF(k) - Om(mu)
            reg = SRG_reg(flow,dem)

            H(p             ,nOrb+n2h1p+iab) = H(p             ,nOrb+n2h1p+iab) + num*reg
            H(nOrb+n2h1p+iab,p             ) = H(nOrb+n2h1p+iab,p             ) + num*reg

            num = sqrt(2d0)*ERI(p,c,k,a)*rho(k,c,mu)
            dem = eHF(c) - eHF(k) + Om(mu)
            reg = SRG_reg(flow,dem)

            H(p             ,nOrb+n2h1p+iab) = H(p             ,nOrb+n2h1p+iab) + num*reg
            H(nOrb+n2h1p+iab,p             ) = H(nOrb+n2h1p+iab,p             ) + num*reg

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

      H(nOrb+ija,nOrb+ija) = eHF(i) - Om(mu) 

     ! First-order terms

      klc = 0
      do k=nC+1,nO
        do nu=1,nS
          klc = klc + 1

          do r=nC+1,nOrb-nR

            num = 0.5d0*rho(k,r,mu)*rho(i,r,nu)
            dem = eHF(i) - eHF(r) + Om(nu)
            reg = SRG_reg(flow,dem)

            H(nOrb+ija,nOrb+klc) = H(nOrb+ija,nOrb+klc) + num*reg

            num = 0.5d0*rho(k,r,mu)*rho(i,r,nu)
            dem = eHF(k) - eHF(r) + Om(mu)
            reg = SRG_reg(flow,dem)

            H(nOrb+ija,nOrb+klc) = H(nOrb+ija,nOrb+klc) + num*reg

          end do

        end do
      end do

    end do
  end do

  !------------------!
  ! Block (K+C)_2p1h !
  !------------------!

  iab = 0
  do a=nO+1,nOrb-nR
    do mu=1,nS
      iab = iab + 1

      H(nOrb+n2h1p+iab,nOrb+n2h1p+iab) = eHF(a) + Om(mu)

      kcd = 0
      do c=nO+1,nOrb-nR
        do nu=1,nS
          kcd = kcd + 1

          do r=nC+1,nOrb-nR

            num = 0.5d0*rho(r,c,mu)*rho(r,a,nu)
            dem = eHF(c) - eHF(r) - Om(mu)
            reg = SRG_reg(flow,dem)

            H(nOrb+n2h1p+iab,nOrb+n2h1p+kcd) = H(nOrb+n2h1p+iab,nOrb+n2h1p+kcd) + num*reg

            num = 0.5d0*rho(r,c,mu)*rho(r,a,nu)
            dem = eHF(a) - eHF(r) - Om(nu)
            reg = SRG_reg(flow,dem)

            H(nOrb+n2h1p+iab,nOrb+n2h1p+kcd) = H(nOrb+n2h1p+iab,nOrb+n2h1p+kcd) + num*reg

          end do

        end do
      end do

    end do
  end do

  !-------------------!
  ! Block C_2h1p-2p1h !
  !-------------------!

  ija = 0
  do i=nC+1,nO
    do mu=1,nS
      ija = ija + 1

      kcd = 0
      do a=nO+1,nOrb-nR
        do nu=1,nS
          kcd = kcd + 1

          ! First-order terms

            do k=nC+1,nO

              num = 1d0*rho(k,i,nu)*rho(a,k,mu)
              dem = eHF(a) - eHF(k) + Om(nu)
              reg = SRG_reg(flow,dem)

              H(nOrb+ija      ,nOrb+n2h1p+kcd) = H(nOrb+ija      ,nOrb+n2h1p+kcd) + num*reg

              H(nOrb+n2h1p+kcd,nOrb+ija      ) = H(nOrb+n2h1p+kcd,nOrb+ija      ) + num*reg

            end do

            do c=nO+1,nOrb-nR

              num = 1d0*rho(a,c,nu)*rho(c,i,mu)
              dem = eHF(i) - eHF(c) - Om(mu)
              reg = SRG_reg(flow,dem)

              H(nOrb+ija      ,nOrb+n2h1p+kcd) = H(nOrb+ija      ,nOrb+n2h1p+kcd) + num*reg

              H(nOrb+n2h1p+kcd,nOrb+ija      ) = H(nOrb+n2h1p+kcd,nOrb+ija      ) + num*reg

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
  write(*,'(1X,A45)')'| ADC(3)-G3W2 energies for all orbitals     |'
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

end subroutine 
