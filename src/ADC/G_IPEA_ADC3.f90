subroutine G_IPEA_ADC3(dotest,nBas,nBas2,nC,nO,nV,nR,nS,ENuc,EGHF,ERI,eHF)

! Dyson version of ADC(3)

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
  write(*,*)'****************************************'
  write(*,*)'* Generalized IP/EA-ADC(3) Calculation *'
  write(*,*)'****************************************'
  write(*,*)

! Dimension of the supermatrix

! Note that ADC(3) is implemented using i<j and a<b restriction while ADC(2) is not.
  n2h1p = nO*(nO-1)*nV/2
  n2p1h = nV*(nV-1)*nO/2
  nH = nBas2 + n2h1p + n2p1h

! Memory allocation

  allocate(H(nH,nH),eGF(nH),Z(nH))

  eF = 0.5d0*(eHF(nO) + eHF(nO+1))

  H(:,:) = 0d0

  call vecout(nO,eHF)

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

  do p=1,nBas2
     H(p,p) = eHF(p)
  end do   

  !--------------!
  ! Block U_2h1p !
  !--------------!
  do p=1,nBas2
     
     ija = 0
     do i=nC+1,nO
        do j=i+1,nO
           do a=nO+1,nBas2-nR
              ija = ija + 1

              ! First-order contribution          
              H(p    ,nBas2+ija) = (ERI(p,a,i,j) - ERI(p,a,j,i))
              ! Second-order contribution        
              do c=nO+1,nBas2-nR
                 do d=nO+1,nBas2-nR
                    H(p    ,nBas2+ija) = H(p    ,nBas2+ija) - 0.5d0 * (ERI(p,a,c,d) - ERI(p,a,d,c)) * (ERI(c,d,i,j) - ERI(c,d,j,i)) / (eHF(c) + eHF(d) - eHF(i) - eHF(j))
                 end do
              end do
              
              do k=nC+1,nO
                 do c=nO+1,nBas2-nR
                    H(p    ,nBas2+ija) = H(p    ,nBas2+ija) + (ERI(p,k,c,j) - ERI(p,k,j,c)) * (ERI(c,a,i,k) - ERI(c,a,k,i)) / (eHF(c) + eHF(a) - eHF(i) - eHF(k))
                    H(p    ,nBas2+ija) = H(p    ,nBas2+ija) - (ERI(p,k,c,i) - ERI(p,k,i,c)) * (ERI(c,a,j,k) - ERI(c,a,k,j)) / (eHF(c) + eHF(a) - eHF(j) - eHF(k))
                 end do
              end do
              
              ! Symmetrize
              H(nBas2+ija,p    ) = H(p    ,nBas2+ija)
              
           end do
        end do
     end do
   
  end do ! p

  !--------------!
  ! Block U_2p1h !
  !--------------!     

  do p=1,nBas2
     iab = 0
     do i=nC+1,nO
        do a=nO+1,nBas2-nR
           do b=a+1,nBas2-nR
              iab = iab + 1   
              
              ! First-order contribution
              H(p          ,nBas2+n2h1p+iab) = (ERI(p,i,a,b) - ERI(p,i,b,a))
              ! Second-order contribution
              do k=nC+1,nO
                 do l=nC+1,nO
                    H(p    ,nBas2+n2h1p+iab) = H(p    ,nBas2+n2h1p+iab) + 0.5d0 * (ERI(p,i,k,l) - ERI(p,i,l,k)) * (ERI(k,l,a,b) - ERI(k,l,b,a)) / (eHF(k) + eHF(l) - eHF(a) - eHF(b))
                 end do
              end do
              
              do k=nC+1,nO
                 do c=nO+1,nBas2-nR
                    H(p    ,nBas2+n2h1p+iab) = H(p    ,nBas2+n2h1p+iab) - (ERI(p,c,k,b) - ERI(p,c,b,k)) * (ERI(k,i,a,c) - ERI(k,i,c,a)) / (eHF(a) + eHF(c) - eHF(i) - eHF(k))
                    H(p    ,nBas2+n2h1p+iab) = H(p    ,nBas2+n2h1p+iab) + (ERI(p,c,k,a) - ERI(p,c,a,k)) * (ERI(k,i,b,c) - ERI(k,i,c,b)) / (eHF(b) + eHF(c) - eHF(i) - eHF(k))
                 end do
              end do
              ! Symmetrize
              H(nBas2+n2h1p+iab,p          ) = H(p          ,nBas2+n2h1p+iab)
           end do
        end do
     end do
     
  end do ! p

   !--------------!
   ! Block K_2h1p !
   !--------- ----!
   ija = 0
   do i=nC+1,nO
     do j=i+1,nO
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
       do b=a+1,nBas2-nR
         iab = iab + 1
              
         H(nBas2+n2h1p+iab,nBas2+n2h1p+iab) = eHF(a) + eHF(b) - eHF(i)
       
       end do
     end do
   end do

   !--------------!
   ! Block C_2h1p !
   !--------- ----!
   ! ija = 0
   ! do i=nC+1,nO
   !    do j=i+1,nO
   !       do a=nO+1,nBas2-nR
   !          ija = ija + 1
            
   !          klc = 0
   !          do k=nC+1,nO
   !             do l=k+1,nO
   !                do c=nO+1,nBas2-nR
   !                   klc = klc + 1
            
   !                   H(nBas2+ija,nBas2+klc) = H(nBas2+ija,nBas2+klc) &
   !                                          - kronecker_delta(a,c) * (ERI(i,j,k,l) - ERI(i,j,l,k)) &
   !                                          + kronecker_delta(i,k) * (ERI(c,j,a,l) - ERI(c,j,l,a)) &
   !                                          + kronecker_delta(j,l) * (ERI(c,i,a,k) - ERI(c,i,k,a)) &
   !                                          - kronecker_delta(i,l) * (ERI(c,j,a,k) - ERI(c,j,k,a)) &
   !                                          - kronecker_delta(j,k) * (ERI(c,i,a,l) - ERI(c,i,l,a))

   !                end do
   !             end do
   !          end do
            
   !       end do
   !    end do
   ! end do
 
   !--------------!
   ! Block C_2p1h !
   !--------------!
   iab = 0
   do i=nC+1,nO
      do a=nO+1,nBas2-nR
         do b=a+1,nBas2-nR
            iab = iab + 1
            
            ! kcd = 0
            ! do k=nC+1,nO
            !    do c=nO+1,nBas2-nR
            !       do d=c+1,nBas2-nR
            !          kcd = kcd + 1
         
            !          H(nBas2+n2h1p+iab,nBas2+n2h1p+kcd) = H(nBas2+n2h1p+iab,nBas2+n2h1p+kcd) &
            !                                             + kronecker_delta(i,k) * (ERI(a,b,c,d) - ERI(a,b,d,c)) &
            !                                             - kronecker_delta(a,c) * (ERI(k,b,i,d) - ERI(k,b,d,i)) &
            !                                             - kronecker_delta(b,d) * (ERI(k,a,i,c) - ERI(k,a,c,i)) &
            !                                             + kronecker_delta(a,d) * (ERI(k,b,i,c) - ERI(k,b,c,i)) &
            !                                             + kronecker_delta(b,c) * (ERI(k,a,i,d) - ERI(k,a,d,i))

            !       end do
            !    end do
            ! end do
            
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
   write(*,'(1X,A38,I3,A2)')'| IPEA-ADC(3) energies (eV) |'
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
  
end subroutine G_IPEA_ADC3
