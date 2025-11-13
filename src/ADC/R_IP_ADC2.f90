subroutine R_IP_ADC2(dotest,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,eHF)

! Non-Dyson version of IP-ADC(2)

  implicit none
  include 'parameters.h'
  
! Input variables

  logical,intent(in)            :: dotest

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
  integer                       :: klc,kcd,ija,ijb,iab,jab

  integer                       :: n2h1p,nH
  double precision,external     :: Kronecker_delta
  double precision,allocatable  :: H(:,:)
  double precision,allocatable  :: eGF(:)
  double precision,allocatable  :: Z(:)

  logical                       :: verbose = .false.
  double precision,parameter    :: cutoff1 = 0.01d0
  double precision,parameter    :: cutoff2 = 0.01d0
  double precision              :: eF
  double precision,parameter    :: window = 2.5d0
  
  double precision,allocatable  :: Reigv(:,:) ! Right eigenvectors

  double precision              :: start_timing,end_timing,timing

! Output variables

! Hello world

  write(*,*)
  write(*,*)'************************************'
  write(*,*)'* Restricted IP-ADC(2) Calculation *'
  write(*,*)'************************************'
  write(*,*)

! Dimension of the supermatrix

  n2h1p = nO*nO*nV
  nH = 1 + n2h1p 

! Memory allocation

  allocate(H(nH,nH),eGF(nH),Z(nH),Reigv(nH,nH))

  eF = 0.5d0*(eHF(nO) + eHF(nO+1))

!-------------------------!
! Main loop over orbitals !
!-------------------------!

  do p=nO,nO

     H(:,:) = 0d0
     Reigv(:,:) = 0d0

    !---------------------------------!
    ! Compute IP-ADC2 supermatrix     !
    !---------------------------------!
    !                                 !
    !     | (K+C)_1h    (C)_1h-2h1p | ! 
    ! H = |                         | ! 
    !     | (C)_2h1p-1h (K+C)_2h1p  | ! 
    !                                 !
    !---------------------------------!

    call wall_time(start_timing)

    !----------------!
    ! Block (K+C)_1h !
    !----------------!

    H(1,1) = eHF(p) 

    do i=nC+1,nO
      do a=nO+1,nOrb-nR
        do b=nO+1,nOrb-nR

        H(1,1) = H(1,1) - (2d0*ERI(p,i,a,b) - ERI(p,i,b,a))*ERI(p,i,a,b)/(eHF(a) + eHF(b) - eHF(i) - eHF(p))

        end do
      end do
    end do

    !------------------!
    ! Block (K+C)_2h1p !
    !------------------!
 
    ija = 0
    do i=nC+1,nO
      do j=nC+1,nO
        do a=nO+1,nOrb-nR
          ija = ija + 1
               
          H(1+ija,1+ija) = eHF(i) + eHF(j) - eHF(a)
 
        end do
      end do
    end do
 
    !-------------------!
    ! Block (C)_2h1p-1h !
    !-------------------!
      
    ija = 0
    do i=nC+1,nO
      do j=nC+1,nO
        do a=nO+1,nOrb-nR
          ija = ija + 1
                   
          H(1+ija,1    ) = 2d0*ERI(p,a,i,j) - ERI(p,a,j,i)
          H(1    ,1+ija) = ERI(p,a,i,j)
        
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

    call diagonalize_general_matrix(nH,H,eGF,Reigv)
 
    call wall_time(end_timing)

    timing = end_timing - start_timing
    write(*,*)
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for diagonalization of supermatrix = ',timing,' seconds'
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
    write(*,'(1X,A37,I3,A3)')'| IP-ADC(2) energies (eV) for orbital',p,'  |'
    write(*,*)'-------------------------------------------'
    write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X)') &
              '|','#','|','e_QP','|','Z','|'
    write(*,*)'-------------------------------------------'
  
    do s=1,nH
      if(eGF(s) < eF .and. eGF(s) > eF - window) then
      ! if(Z(s) > cutoff1) then
        write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
        '|',s,'|',eGF(s)*HaToeV,'|',Z(s),'|'
      end if
    end do
  
    write(*,*)'-------------------------------------------'
    write(*,*)
 
    if(verbose) then

      do s=1,nH
       
        if(eGF(s) < eF .and. eGF(s) > eF - window) then
       
          write(*,*)'------------------------------------------------------------------------------'
          write(*,'(1X,A7,1X,I3,A6,I3,A1,1X,A7,F12.6,A13,F6.4,1X)') & 
               'Orbital',p,' and #',s,':','e_QP = ',eGF(s)*HaToeV,' eV and Z = ',Z(s)
          write(*,*)'------------------------------------------------------------------------------'
          write(*,'(1X,A20,1X,A20,1X,A15,1X,A20,1X)') &
               ' Configuration ',' Coefficient ',' Weight ',' Zeroth-order ' 
          write(*,*)'------------------------------------------------------------------------------'
          
          if(p <= nO) & 
               write(*,'(1X,A7,I3,A16,1X,F15.6,1X,F15.6,1X,F12.6)') &
               '      (',p,')               ',Reigv(1,s),Reigv(1,s)**2,-eHF(p)*HaToeV
    
          ija = 0
          do i=nC+1,nO
            do j=nC+1,nO
              do a=nO+1,nOrb-nR
                ija = ija + 1
  
                if(abs(Reigv(1+ija,s)) > cutoff2)               &
                     write(*,'(1X,A3,I3,A1,I3,A6,I3,A7,1X,F15.6,1X,F15.6,1X,F12.6)') &
                     '  (',i,',',j,') -> (',a,')      ',Reigv(1+ija,s),Reigv(1+ija,s)**2, & 
                                                        (eHF(i) + eHF(j) - eHF(a))*HaToeV
           
              end do
            end do
          end do
           
          write(*,*)'------------------------------------------------------------------------------'
          write(*,*)

        end if ! If state s should be print

      end do ! Loop on s
       
    end if ! If verbose

  end do ! Loop on the orbital in the e block
  
end subroutine 
