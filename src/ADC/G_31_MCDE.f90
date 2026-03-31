subroutine G_31_MCDE(dotest,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,EGHF,ERI,eHF)

! (3,1) multichannel Dyson Equation

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
  double precision,intent(in)   :: EGHF
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: eHF(nOrb)

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
  double precision,parameter    :: cutoff1 = 0.1d0
  double precision,parameter    :: cutoff2 = 0.01d0
  double precision              :: eF
  double precision,parameter    :: window = 1.5d0

  double precision              :: start_timing,end_timing,timing

! Output variables

! Hello world

  write(*,*)
  write(*,*)'**************************************'
  write(*,*)'* Generalized (3,1)-MCDE Calculation *'
  write(*,*)'**************************************'
  write(*,*)

! Dimension of the supermatrix

  n2h1p = nO*(nO-1)*nV/2
  n2p1h = nV*(nV-1)*nO/2
  nH = nOrb + n2h1p + n2p1h

! Memory allocation

  allocate(H(nH,nH),eGF(nH),Z(nH))

  eF = 0.5d0*(eHF(nO) + eHF(nO+1))

  H(:,:) = 0d0

  !--------------------------------------!
  ! Compute IP/EA-ADC(3) supermatrix     !
  !--------------------------------------!
  !                                      !
  !     |   F    U_2h1p   U_2p1h   |     ! 
  !     |                          |     ! 
  ! H = | U_2h1p K+C_2h1p 0        |     ! 
  !     |                          |     ! 
  !     | U_2p1h 0        K+C_2p1h |     ! 
  !                                      !
  !--------------------------------------!

  call wall_time(start_timing)

  !---------!
  ! Block F !
  !---------!

  do p=1,nOrb
     H(p,p) = eHF(p)
  end do   

  !--------------!
  ! Block U_2h1p !
  !--------------!

  do p=1,nOrb
     
     ija = nOrb
     do i=nC+1,nO
        do j=i+1,nO
           do a=nO+1,nOrb-nR
              ija = ija + 1

              ! First-order contribution          

              H(p,ija) = ERI(i,j,p,a) - ERI(i,j,a,p)

              ! Symmetrize

              H(ija,p) = H(p,ija)
              
           end do
        end do
     end do
   
  end do ! p

  !--------------!
  ! Block U_2p1h !
  !--------------!     

  do p=1,nOrb
     iab = nOrb + n2h1p
     do i=nC+1,nO
        do a=nO+1,nOrb-nR
           do b=a+1,nOrb-nR
              iab = iab + 1   
              
              ! First-order contribution

              H(p,iab) = ERI(a,b,p,i) - ERI(a,b,i,p)

              ! Symmetrize

              H(iab,p) = H(p,iab)

           end do
        end do
     end do
     
  end do ! p

   !--------------!
   ! Block K_2h1p !
   !--------------!

   ija = nOrb
   do i=nC+1,nO
     do j=i+1,nO
       do a=nO+1,nOrb-nR
         ija = ija + 1
            
         H(ija,ija) = eHF(i) + eHF(j) - eHF(a)

       end do
     end do
   end do
 
   !--------------!
   ! Block K_2p1h !
   !--------------!

   iab = nOrb + n2h1p
   do i=nC+1,nO
     do a=nO+1,nOrb-nR
       do b=a+1,nOrb-nR
         iab = iab + 1
              
         H(iab,iab) = eHF(a) + eHF(b) - eHF(i)
       
       end do
     end do
   end do

   !--------------!
   ! Block C_2h1p !
   !--------------!

   ija = nOrb
   do i=nC+1,nO
      do j=i+1,nO
         do a=nO+1,nOrb-nR
            ija = ija + 1
          
            klc = nOrb
            do k=nC+1,nO
               do l=k+1,nO
                  do c=nO+1,nOrb-nR
                     klc = klc + 1
          
                     H(ija,klc) = H(ija,klc) &
                                - kronecker_delta(a,c) * (ERI(i,j,k,l) - ERI(i,j,l,k)) &
                                + kronecker_delta(i,k) * (ERI(c,j,a,l) - ERI(c,j,l,a)) &
                                + kronecker_delta(j,l) * (ERI(c,i,a,k) - ERI(c,i,k,a)) &
                                - kronecker_delta(i,l) * (ERI(c,j,a,k) - ERI(c,j,k,a)) &
                                - kronecker_delta(j,k) * (ERI(c,i,a,l) - ERI(c,i,l,a))

                  end do
               end do
            end do
          
         end do
      end do
   end do
 
   !--------------!
   ! Block C_2p1h !
   !--------------!

   iab = nOrb + n2h1p
   do i=nC+1,nO
      do a=nO+1,nOrb-nR
         do b=a+1,nOrb-nR
            iab = iab + 1
            
              kcd = nOrb + n2h1p
              do k=nC+1,nO
                 do c=nO+1,nOrb-nR
                    do d=c+1,nOrb-nR
                       kcd = kcd + 1
         
                       H(iab,kcd) = H(iab,kcd) &
                                  + kronecker_delta(i,k) * (ERI(a,b,c,d) - ERI(a,b,d,c)) &
                                  - kronecker_delta(a,c) * (ERI(k,b,i,d) - ERI(k,b,d,i)) &
                                  - kronecker_delta(b,d) * (ERI(k,a,i,c) - ERI(k,a,c,i)) &
                                  + kronecker_delta(a,d) * (ERI(k,b,i,c) - ERI(k,b,c,i)) &
                                  + kronecker_delta(b,c) * (ERI(k,a,i,d) - ERI(k,a,d,i))

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
      do p=nC+1,nOrb-nR
         Z(s) = Z(s) + H(p,s)**2
      end do
   end do

   !--------------!
   ! Dump results !
   !--------------!
   
   write(*,*)'-------------------------------------------'
   write(*,'(1X,A43)')'| (3,1)-MCDE energies for all orbitals    |'
   write(*,*)'-------------------------------------------'
   write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X)') &
        '|','#','|','e_QP (eV)','|','Z','|'
   write(*,*)'-------------------------------------------'
   
   do s=1,nH
!     if(eGF(s) < eF .and. eGF(s) > eF - window) then
      if(Z(s) > cutoff1) then
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
          do p=nO+1,nOrb-nR
            if(abs(H(p,s)) > cutoff2)                 &
          write(*,'(1X,A16,I3,A7,1X,F15.6,1X,F15.6)') &
          '               (',p,')      ',H(p,s),H(p,s)**2
          end do

          ija = 0
          do i=nC+1,nO
            do j=i+1,nO
              do a=nO+1,nOrb-nR
                ija = ija + 1
 
                if(abs(H(nOrb+ija,s)) > cutoff2)               &
                write(*,'(1X,A3,I3,A1,I3,A6,I3,A7,1X,F15.6,1X,F15.6)') &
                '  (',i,',',j,') -> (',a,')      ',H(nOrb+ija,s),H(nOrb+ija,s)**2
         
              end do
            end do
          end do
         
          iab = 0
          do i=nC+1,nO
            do a=nO+1,nOrb-nR
              do b=a+1,nOrb-nR
                iab = iab + 1

                if(abs(H(nOrb+n2h1p+iab,s)) > cutoff2)           &
                  write(*,'(1X,A7,I3,A6,I3,A1,I3,A3,1X,F15.6,1X,F15.6)') &
                  '      (',i,') -> (',a,',',b,')  ',H(nOrb+n2h1p+iab,s),H(nOrb+n2h1p+iab,s)**2
                
              end do
            end do
          end do

          write(*,*)'-------------------------------------------------------------'
          write(*,*)

        end if

      end do
       
    end if ! If verbose
  
end subroutine 
