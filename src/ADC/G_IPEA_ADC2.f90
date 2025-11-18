subroutine G_IPEA_ADC2_single_state(dotest,nBas,nBas2,nC,nO,nV,nR,nS,ENuc,EGHF,ERI,eHF)

! Dyson version of IP/EA-ADC(2)

  implicit none
  include 'parameters.h'
  
! Input variables

  logical,intent(in)            :: dotest

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nBas2
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: EGHF
  double precision,intent(in)   :: ERI(nBas2,nBas2,nBas2,nBas2)
  double precision,intent(in)   :: eHF(nBas2)

! Local variables

  integer                       :: p
  integer                       :: s
  integer                       :: i,j,k,l
  integer                       :: a,b,c,d
  integer                       :: jb,kc,ia,ja
  integer                       :: klc,kcd,ija,ijb,iab,jab

  integer                       :: n2h1p,n2p1h,nH
  double precision,external     :: Kronecker_delta
  double precision,allocatable  :: H(:,:)
  double precision,allocatable  :: eGF(:)
  double precision,allocatable  :: Z(:)

  logical                       :: verbose = .false.
  double precision,parameter    :: cutoff1 = 0.05d0
  double precision,parameter    :: cutoff2 = 0.05d0
  double precision              :: eF
  double precision,parameter    :: window = 0.5d0

  double precision              :: start_timing,end_timing,timing

! Output variables

! Hello world

  write(*,*)
  write(*,*)'*****************************************************'
  write(*,*)'* Generalized single-state IP/EA-ADC(2) Calculation *'
  write(*,*)'*****************************************************'
  write(*,*)

! Dimension of the supermatrix

  n2h1p = nO*nO*nV
  n2p1h = nV*nV*nO
  nH = 1 + n2h1p + n2p1h

! Memory allocation

  allocate(H(nH,nH),eGF(nH),Z(nH))

  eF = 0.5d0*(eHF(nO) + eHF(nO+1))

!-------------------------!
! Main loop over orbitals !
!-------------------------!

  do p=nO,nO

     H(:,:) = 0d0

    !--------------------------------------!
    ! Compute IP/EA-ADC(2) supermatrix     !
    !--------------------------------------!
    !                                      !
    !     |   e    U_2h1p U_2p1h |         ! 
    !     |                      |         ! 
    ! H = | U_2h1p K_2h1p 0      |         ! 
    !     |                      |         ! 
    !     | U_2p1h 0      K_2p1h |         ! 
    !                                      !
    !--------------------------------------!

    call wall_time(start_timing)

    !---------!
    ! Block e !
    !---------!
      
    H(1,1) = eHF(p)
    
    !--------------!
    ! Block U_2h1p !
    !--------------!

    ija = 0
    do i=nC+1,nO
      do j=nC+1,nO
        do a=nO+1,nBas2-nR
          ija = ija + 1
             
          H(1    ,1+ija) = (ERI(p,a,i,j) - ERI(p,a,j,i))/sqrt(2d0)
          H(1+ija,1    ) = (ERI(p,a,i,j) - ERI(p,a,j,i))/sqrt(2d0)
             
        end do
      end do
    end do

    !--------------!
    ! Block U_2p1h !
    !--------------!     
 
    iab = 0
    do i=nC+1,nO
      do a=nO+1,nBas2-nR
        do b=nO+1,nBas2-nR
          iab = iab + 1   
 
          H(1          ,1+n2h1p+iab) = (ERI(p,i,a,b) - ERI(p,i,b,a))/sqrt(2d0)
          H(1+n2h1p+iab,1          ) = (ERI(p,i,a,b) - ERI(p,i,b,a))/sqrt(2d0)
               
          end do
        end do
      end do
 
    !--------------!
    ! Block K_2h1p !
    !--------- ----!
 
    ija = 0
    do i=nC+1,nO
      do j=nC+1,nO
        do a=nO+1,nBas2-nR
          ija = ija + 1
               
          H(1+ija,1+ija) = eHF(i) + eHF(j) - eHF(a)

 
        end do
      end do
    end do
 
    !--------------!
    ! Block K_2p1h !
    !--------------!
      
    iab = 0
    do i=nC+1,nO
      do a=nO+1,nBas2-nR
        do b=nO+1,nBas2-nR
          iab = iab + 1
               
          H(1+n2h1p+iab,1+n2h1p+iab) = eHF(a) + eHF(b) - eHF(i)
        
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

    call diagonalize_matrix(nH,H,eGF)
 
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
    write(*,'(1X,A38,I3,A2)')'| IPEA-ADC(2) energies (eV) for orbital',p,' |'
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
               '      (',p,')               ',H(1,s),H(1,s)**2,-eHF(p)*HaToeV
          if(p > nO) & 
               write(*,'(1X,A16,I3,A7,1X,F15.6,1X,F15.6)') &
               '               (',p,')      ',H(1,s),H(1,s)**2,-eHF(p)*HaToeV
    
          ija = 0
          do i=nC+1,nO
            do j=nC+1,nO
              do a=nO+1,nBas2-nR
                ija = ija + 1
  
                if(abs(H(1+ija,s)) > cutoff2)               &
                     write(*,'(1X,A3,I3,A1,I3,A6,I3,A7,1X,F15.6,1X,F15.6,1X,F12.6)') &
                     '  (',i,',',j,') -> (',a,')      ',H(1+ija,s),H(1+ija,s)**2, & 
                                                        (eHF(i) + eHF(j) - eHF(a))*HaToeV
           
              end do
            end do
          end do
           
          iab = 0
          do i=nC+1,nO
            do a=nO+1,nBas2-nR
              do b=nO+1,nBas2-nR
                iab = iab + 1
 
                if(abs(H(1+n2h1p+iab,s)) > cutoff2)           &
                     write(*,'(1X,A7,I3,A6,I3,A1,I3,A3,1X,F15.6,1X,F15.6,1X,F12.6)') &
                     '      (',i,') -> (',a,',',b,')  ',H(1+n2h1p+iab,s),H(1+n2h1p+iab,s)**2, & 
                                                        (eHF(a) + eHF(b) - eHF(i))*HaToeV
                  
              end do
            end do
          end do

          write(*,*)'------------------------------------------------------------------------------'
          write(*,*)

        end if ! If state s should be print

      end do ! Loop on s
       
    end if ! If verbose

  end do ! Loop on the orbital in the e block
  
end subroutine G_IPEA_ADC2_single_state

subroutine G_IPEA_ADC2(dotest,nBas,nBas2,nC,nO,nV,nR,nS,ENuc,EGHF,ERI,eHF)

! Dyson version of IP/EA-ADC(2)

  implicit none
  include 'parameters.h'
  
! Input variables

  logical,intent(in)            :: dotest

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nBas2
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: EGHF
  double precision,intent(in)   :: ERI(nBas2,nBas2,nBas2,nBas2)
  double precision,intent(in)   :: eHF(nBas2)

! Local variables

  integer                       :: p
  integer                       :: s
  integer                       :: i,j,k,l
  integer                       :: a,b,c,d
  integer                       :: jb,kc,ia,ja
  integer                       :: klc,kcd,ija,ijb,iab,jab

  integer                       :: n2h1p,n2p1h,nH
  double precision,external     :: Kronecker_delta
  double precision,allocatable  :: H(:,:)
  double precision,allocatable  :: eGF(:)
  double precision,allocatable  :: Z(:)

  logical                       :: verbose = .false.
  double precision,parameter    :: cutoff1 = 0.05d0
  double precision,parameter    :: cutoff2 = 0.05d0
  double precision              :: eF
  double precision,parameter    :: window = 1.5d0

  double precision              :: start_timing,end_timing,timing

! Output variables

! Hello world

  write(*,*)
  write(*,*)'****************************************'
  write(*,*)'* Generalized IP/EA-ADC(2) Calculation *'
  write(*,*)'****************************************'
  write(*,*)

! Dimension of the supermatrix

  n2h1p = nO*nO*nV
  n2p1h = nV*nV*nO
  nH = nBas2 + n2h1p + n2p1h

! Memory allocation

  allocate(H(nH,nH),eGF(nH),Z(nH))

  eF = 0.5d0*(eHF(nO) + eHF(nO+1))

  H(:,:) = 0d0

  call vecout(nO,eHF)

  !--------------------------------------!
  ! Compute IP/EA-ADC(2) supermatrix     !
  !--------------------------------------!
  !                                      !
  !     |   F    U_2h1p U_2p1h |         ! 
  !     |                      |         ! 
  ! H = | U_2h1p K_2h1p 0      |         ! 
  !     |                      |         ! 
  !     | U_2p1h 0      K_2p1h |         ! 
  !                                      !
  !--------------------------------------!

  call wall_time(start_timing)

  !---------!
  ! Block F !
  !---------!

  do p=1,nBas2
     H(p,p) = eHF(p)
  end do   

  !--------------!
  ! Block U_2h1p !
  !--------------!

  do p=1,nBas2
    ija = 0
    do i=nC+1,nO
      do j=nC+1,nO
        do a=nO+1,nBas2-nR
          ija = ija + 1
             
          H(p    ,nBas2+ija) = (ERI(p,a,i,j) - ERI(p,a,j,i))/sqrt(2d0)
          H(nBas2+ija,p    ) = (ERI(p,a,i,j) - ERI(p,a,j,i))/sqrt(2d0)
             
        end do
      end do
    end do
  end do   

  !--------------!
  ! Block U_2p1h !
  !--------------!     

  do p=1,nBas2
    iab = 0
    do i=nC+1,nO
      do a=nO+1,nBas2-nR
        do b=nO+1,nBas2-nR
          iab = iab + 1   
 
          H(p          ,nBas2+n2h1p+iab) = (ERI(p,i,a,b) - ERI(p,i,b,a))/sqrt(2d0)
          H(nBas2+n2h1p+iab,p          ) = (ERI(p,i,a,b) - ERI(p,i,b,a))/sqrt(2d0)
               
          end do
        end do
      end do
   end do

   !--------------!
   ! Block K_2h1p !
   !--------- ----!


   ija = 0
   do i=nC+1,nO
     do j=nC+1,nO
       do a=nO+1,nBas2-nR
         ija = ija + 1
            
         H(nBas2+ija,nBas2+ija) = eHF(i) + eHF(j) - eHF(a)

       end do
     end do
   end do
 
   !--------------!
   ! Block K_2p1h !
   !--------------!
     
   iab = 0
   do i=nC+1,nO
     do a=nO+1,nBas2-nR
       do b=nO+1,nBas2-nR
         iab = iab + 1
              
         H(nBas2+n2h1p+iab,nBas2+n2h1p+iab) = eHF(a) + eHF(b) - eHF(i)
       
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

   call diagonalize_matrix(nH,H,eGF)
 
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
      do p=nC+1,nBas2-nR
         Z(s) = Z(s) + H(p,s)**2
      end do
   end do

   !--------------!
   ! Dump results !
   !--------------!
   
   write(*,*)'-------------------------------------------'
   write(*,'(1X,A38,I3,A2)')'| IPEA-ADC(2) energies (eV) |'
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
 
          write(*,*)'-------------------------------------------------------------'
          write(*,'(1X,A12,1X,I3,A1,1X,A7,F12.6,A13,F6.4,1X)') & 
           'Eigenvalue #',s,':','e_QP = ',eGF(s)*HaToeV,' eV and Z = ',Z(s)
          write(*,*)'-------------------------------------------------------------'
          write(*,'(1X,A20,1X,A20,1X,A15,1X)') &
                    ' Configuration ',' Coefficient ',' Weight ' 
          write(*,*)'-------------------------------------------------------------'
         
          do p=nC+1,nO 
            if(abs(H(p,s)) > cutoff2)                   &
            write(*,'(1X,A7,I3,A16,1X,F15.6,1X,F15.6)') &
            '      (',p,')               ',H(p,s),H(p,s)**2
          end do
          do p=nO+1,nBas2-nR
            if(abs(H(p,s)) > cutoff2)                 &
          write(*,'(1X,A16,I3,A7,1X,F15.6,1X,F15.6)') &
          '               (',p,')      ',H(p,s),H(p,s)**2
          end do

          ija = 0
          do i=nC+1,nO
            do j=nC+1,nO
              do a=nO+1,nBas2-nR
                ija = ija + 1
 
                if(abs(H(nBas2+ija,s)) > cutoff2)               &
                write(*,'(1X,A3,I3,A1,I3,A6,I3,A7,1X,F15.6,1X,F15.6)') &
                '  (',i,',',j,') -> (',a,')      ',H(nBas2+ija,s),H(nBas2+ija,s)**2
         
              end do
            end do
          end do
         
          iab = 0
          do i=nC+1,nO
            do a=nO+1,nBas2-nR
              do b=nO+1,nBas2-nR
                iab = iab + 1

                if(abs(H(nBas2+n2h1p+iab,s)) > cutoff2)           &
                  write(*,'(1X,A7,I3,A6,I3,A1,I3,A3,1X,F15.6,1X,F15.6)') &
                  '      (',i,') -> (',a,',',b,')  ',H(nBas2+n2h1p+iab,s),H(nBas2+n2h1p+iab,s)**2
                
              end do
            end do
          end do

          write(*,*)'-------------------------------------------------------------'
          write(*,*)

        end if

      end do
       
    end if ! If verbose
  
end subroutine G_IPEA_ADC2
