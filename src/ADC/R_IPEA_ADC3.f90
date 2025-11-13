subroutine R_IPEA_ADC3(dotest,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,e)

! Dyson version of IP/EA-ADC(3) 

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
  double precision,intent(in)   :: e(nOrb)

! Local variables

  integer                       :: p
  integer                       :: s
  integer                       :: i,j,k,l
  integer                       :: a,b,c,d
  integer                       :: jb,kc,ia,ja
  integer                       :: akl,jab,bij,icd,ija,iab

  integer                       :: n2h1p,n2p1h,nH
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
  write(*,*)'***************************************'
  write(*,*)'* Restricted IP/EA-ADC(3) Calculation *'
  write(*,*)'***************************************'
  write(*,*)

! Dimension of the supermatrix

  n2h1p = nO*nO*nV
  n2p1h = nV*nV*nO
  nH = 1 + n2h1p + n2p1h

! Memory allocation

  allocate(H(nH,nH),eGF(nH),Z(nH),Reigv(nH,nH))

  eF = 0.5d0*(e(nO) + e(nO+1))

!-------------------------!
! Main loop over orbitals !
!-------------------------!

  do p=nO,nO

     H(:,:) = 0d0
     Reigv(:,:) = 0d0

    !--------------------------------------!
    ! Compute IP/EA-ADC(3) supermatrix     !
    !--------------------------------------!
    !                                      !
    !     |   F      U_2h1p     U_2p1h   | ! 
    !     |                              | ! 
    ! H = | U_2h1p (K+C)_2h1p   0        | ! 
    !     |                              | ! 
    !     | U_2p1h   0        (K+C)_2p1h | ! 
    !                                      !
    !--------------------------------------!

    call wall_time(start_timing)

    !---------!
    ! Block F !
    !---------!
      
    H(1,1) = e(p)

    !--------------!
    ! Block U_2h1p !
    !--------------!

    akl = 0
    do a=nO+1,nOrb-nR
      do k=nC+1,nO
        do l=nC+1,nO
          akl = akl + 1
             
          H(1    ,1+akl) = ERI(k,l,p,a)  - ERI(k,l,a,p)

          H(1+akl,1    ) = ERI(k,l,p,a)  - ERI(k,l,a,p)

          do b=nO+1,nOrb-nR
            do c=nO+1,nOrb-nR
      
              H(1    ,1+akl) = H(1    ,1+akl) & 
                             + 0.5d0*(ERI(k,l,b,c) - ERI(k,l,c,b))/(e(k) + e(l) - e(b) - e(c))*(ERI(b,c,p,a) - ERI(b,c,a,p))

              H(1+akl,1    ) = H(1+akl,1    ) &
                             + 0.5d0*(ERI(k,l,b,c) - ERI(k,l,c,b))/(e(k) + e(l) - e(b) - e(c))*(ERI(b,c,p,a) - ERI(b,c,a,p))

            end do
          end do

          do i=nC+1,nO
            do b=nO+1,nOrb-nR
      
              H(1    ,1+akl) = H(1    ,1+akl) & 
                             - 0.5d0*(ERI(k,i,b,a) - ERI(k,i,a,b))/(e(k) + e(i) - e(b) - e(a))*(ERI(b,l,p,i) - ERI(b,l,i,p)) &
                             + 0.5d0*(ERI(l,i,b,a) - ERI(l,i,a,b))/(e(l) + e(i) - e(b) - e(a))*(ERI(b,k,p,i) - ERI(b,k,i,p))

              H(1+akl,1    ) = H(1+akl,1    ) &
                             - 0.5d0*(ERI(k,i,b,a) - ERI(k,i,a,b))/(e(k) + e(i) - e(b) - e(a))*(ERI(b,l,p,i) - ERI(b,l,i,p)) &
                             + 0.5d0*(ERI(l,i,b,a) - ERI(l,i,a,b))/(e(l) + e(i) - e(b) - e(a))*(ERI(b,k,p,i) - ERI(b,k,i,p))

            end do
          end do
             
        end do
      end do
    end do

    !--------------!
    ! Block U_2p1h !
    !--------------!     
 
    jab = 0
    do j=nC+1,nO
      do a=nO+1,nOrb-nR
        do b=nO+1,nOrb-nR
          jab = jab + 1   
 
          H(1          ,1+n2h1p+jab) = ERI(a,b,p,j) - ERI(a,b,j,p)

          H(1+n2h1p+jab,1          ) = ERI(a,b,p,j) - ERI(a,b,j,p)

          do k=nC+1,nO
            do l=nC+1,nO
      
              H(1          ,1+n2h1p+jab) = H(1          ,1+n2h1p+jab) &
                                         - 0.5d0*(ERI(a,b,k,l) - ERI(a,b,l,k))/(e(a) + e(b) - e(k) - e(l))*(ERI(k,l,p,j) - ERI(k,l,j,p))
              H(1+n2h1p+jab,1          ) = H(1+n2h1p+jab,1          ) &
                                         - 0.5d0*(ERI(a,b,k,l) - ERI(a,b,l,k))/(e(a) + e(b) - e(k) - e(l))*(ERI(k,l,p,j) - ERI(k,l,j,p))

            end do
          end do

          do k=nC+1,nO
            do c=nO+1,nOrb-nR
      
              H(1          ,1+n2h1p+jab) = H(1          ,1+n2h1p+jab) &
                                         + 0.5d0*(ERI(a,c,k,j) - ERI(a,c,j,k))/(e(a) + e(c) - e(k) - e(j))*(ERI(k,b,p,c) - ERI(k,b,c,p)) &
                                         - 0.5d0*(ERI(b,c,k,j) - ERI(b,c,j,k))/(e(b) + e(c) - e(k) - e(j))*(ERI(k,a,p,c) - ERI(k,a,c,p))  
              H(1+n2h1p+jab,1          ) = H(1+n2h1p+jab,1          ) &
                                         + 0.5d0*(ERI(a,c,k,j) - ERI(a,c,j,k))/(e(a) + e(c) - e(k) - e(j))*(ERI(k,b,p,c) - ERI(k,b,c,p)) &
                                         - 0.5d0*(ERI(b,c,k,j) - ERI(b,c,j,k))/(e(b) + e(c) - e(k) - e(j))*(ERI(k,a,p,c) - ERI(k,a,c,p))  

            end do
          end do
               
        end do
      end do
    end do
 
    !------------------!
    ! Block (K+C)_2h1p !
    !------------------!
 
    akl = 0
    do a=nO+1,nOrb-nR
      do k=nC+1,nO
        do l=nC+1,nO
          akl = akl + 1
               
          H(1+akl,1+akl) = e(k) + e(l) - e(a)

          bij = 0
          do b=nO+1,nOrb-nR
            do i=nC+1,nO
              do j=nC+1,nO
                bij = bij + 1
        
                H(1+akl,1+bij) = H(1+akl,1+bij) &
                               - Kronecker_delta(a,b)*(ERI(k,l,i,j) - ERI(k,l,j,i)) &
                               + Kronecker_delta(k,i)*(ERI(b,l,a,j) - ERI(b,l,j,a)) &
                               + Kronecker_delta(l,j)*(ERI(b,k,a,i) - ERI(b,k,i,a)) &
                               - Kronecker_delta(k,j)*(ERI(b,l,a,i) - ERI(b,l,i,a)) &
                               - Kronecker_delta(l,i)*(ERI(b,k,a,j) - ERI(b,k,j,a))  

              end do
            end do
          end do

        end do
      end do
    end do
 
    !------------------!
    ! Block (K+C)_2p1h !
    !------------------!
      
    jab = 0
    do j=nC+1,nO
      do a=nO+1,nOrb-nR
        do b=nO+1,nOrb-nR
          jab = jab + 1
               
          H(1+n2h1p+jab,1+n2h1p+jab) = e(a) + e(b) - e(i)

          icd = 0
          do i=nC+1,nO
            do c=nO+1,nOrb-nR
              do d=nO+1,nOrb-nR
                icd = icd + 1

                H(1+n2h1p+jab,1+n2h1p+icd) = H(1+n2h1p+jab,1+n2h1p+icd) &
                                           + Kronecker_delta(j,i)*(ERI(a,b,c,d) - ERI(a,b,d,c)) &
                                           - Kronecker_delta(a,c)*(ERI(i,b,j,d) - ERI(i,b,d,j)) &
                                           - Kronecker_delta(b,d)*(ERI(i,a,j,c) - ERI(i,a,c,j)) &
                                           + Kronecker_delta(a,d)*(ERI(i,b,j,c) - ERI(i,b,c,j)) &
                                           + Kronecker_delta(b,c)*(ERI(i,a,j,d) - ERI(i,a,d,j))  
        
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
    write(*,'(1X,A38,I3,A2)')'| IPEA-ADC(3) energies (eV) for orbital',p,' |'
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
               '      (',p,')               ',Reigv(1,s),Reigv(1,s)**2,-e(p)*HaToeV
          if(p > nO) & 
               write(*,'(1X,A16,I3,A7,1X,F15.6,1X,F15.6)') &
               '               (',p,')      ',Reigv(1,s),Reigv(1,s)**2,-e(p)*HaToeV
    
          ija = 0
          do i=nC+1,nO
            do j=nC+1,nO
              do a=nO+1,nOrb-nR
                ija = ija + 1
  
                if(abs(Reigv(1+ija,s)) > cutoff2)               &
                     write(*,'(1X,A3,I3,A1,I3,A6,I3,A7,1X,F15.6,1X,F15.6,1X,F12.6)') &
                     '  (',i,',',j,') -> (',a,')      ',Reigv(1+ija,s),Reigv(1+ija,s)**2, & 
                                                        (e(i) + e(j) - e(a))*HaToeV
           
              end do
            end do
          end do
           
          iab = 0
          do i=nC+1,nO
            do a=nO+1,nOrb-nR
              do b=nO+1,nOrb-nR
                iab = iab + 1
 
                if(abs(Reigv(1+n2h1p+iab,s)) > cutoff2)           &
                     write(*,'(1X,A7,I3,A6,I3,A1,I3,A3,1X,F15.6,1X,F15.6,1X,F12.6)') &
                     '      (',i,') -> (',a,',',b,')  ',Reigv(1+n2h1p+iab,s),Reigv(1+n2h1p+iab,s)**2, & 
                                                        (e(a) + e(b) - e(i))*HaToeV
                  
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
