subroutine eomRG0W0(dotest,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,eHF)

! EOM version of G0W0 

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

  logical                       :: print_W = .false.
  logical                       :: dRPA
  integer                       :: isp_W
  double precision              :: EcRPA
  integer                       :: n2h1p,n2p1h,nH
  double precision,external     :: Kronecker_delta
  double precision,allocatable  :: H(:,:)
  double precision,allocatable  :: cGW(:,:)
  double precision,allocatable  :: eGW(:)
  double precision,allocatable  :: Z(:)
  integer,allocatable           :: order(:)

  logical                       :: verbose = .false.
  double precision,parameter    :: cutoff1 = 0.01d0
  double precision,parameter    :: cutoff2 = 0.01d0
  double precision              :: eF
  double precision,parameter    :: window = 2.5d0

  double precision              :: start_timing,end_timing,timing

! Output variables

! Hello world

  write(*,*)
  write(*,*)'***********************************'
  write(*,*)'* Restricted EOM-G0W0 Calculation *'
  write(*,*)'***********************************'
  write(*,*)

! Dimension of the supermatrix

  n2h1p = nO*nO*nV
  n2p1h = nV*nV*nO
  nH = 1 + n2h1p + n2p1h

! Memory allocation

  allocate(H(nH,nH),eGW(nH),cGW(nH,nH),Z(nH),order(nH))
  
! Initialization

  dRPA = .true.
  EcRPA = 0d0

  eF = 0.5d0*(eHF(nO+1) + eHF(nO))

!-------------------------!
! Main loop over orbitals !
!-------------------------!

  do p=nO,nO+1

    H(:,:) = 0d0
 
!-----------------------------------------!
!  Compute BSE supermatrix                !
!-----------------------------------------!
!                                         !
!     |  A    V2h1p V2p1h   0     0   |   ! 
!     |                               |   ! 
!     | V2h1p A2h2p   0   B2h1p   0   |   ! 
!     |                               |   ! 
! H = | V2p1h   0   A2p2h   0   B2p1h |   ! 
!     |                               |   ! 
!     |  0      0     0     0     0   |   ! 
!     |                               |   ! 
!     |  0      0     0     0     0   |   ! 
!                                         !
!-----------------------------------------!
 
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
        do a=nO+1,nOrb-nR
          ija = ija + 1
             
          H(1    ,1+ija) = sqrt(2d0)*ERI(p,a,i,j)
          H(1+ija,1    ) = sqrt(2d0)*ERI(p,a,i,j)
!         H(1+n2h1p+n2p1h+ija,1    ) = sqrt(2d0)*ERI(p,a,i,j)
!         H(1+ija,1+n2h1p+n2p1h    ) = sqrt(2d0)*ERI(p,a,i,j)
 
        end do
      end do
    end do

    !-------------!
    ! Block V2p1h !
    !-------------!     
 
    iab = 0
    do i=nC+1,nO
      do a=nO+1,nOrb-nR
        do b=nO+1,nOrb-nR
          iab = iab + 1   
 
          H(1          ,1+n2h1p+iab) = sqrt(2d0)*ERI(p,i,b,a)
          H(1+n2h1p+iab,1          ) = sqrt(2d0)*ERI(p,i,b,a)
!         H(1          ,1+2*n2h1p+n2p1h+iab) = sqrt(2d0)*ERI(p,i,b,a)
!         H(1+2*n2h1p+n2p1h+iab,1          ) = sqrt(2d0)*ERI(p,i,b,a)
             
        end do
      end do
    end do
 
    !-------------!
    ! Block A2h1p !
    !-------------!
 
    ija = 0
    do i=nC+1,nO
      do j=nC+1,nO
        do a=nO+1,nOrb-nR
          ija = ija + 1
             
          klc = 0
          do k=nC+1,nO
            do l=nC+1,nO
              do c=nO+1,nOrb-nR
                klc = klc + 1
                   
                H(1+ija,1+klc) & 
                     = ((eHF(i) + eHF(j) - eHF(a))*Kronecker_delta(j,l)*Kronecker_delta(a,c) & 
                     - 2d0*ERI(j,c,a,l) - 2d0*ERI(j,l,a,c))*Kronecker_delta(i,k)
                   
!               H(1+n2h1p+n2p1h+ija,1+n2h1p+n2p1h+klc) & 
!                    = ((eHF(i) + eHF(j) - eHF(a))*Kronecker_delta(j,l)*Kronecker_delta(a,c) & 
!                    - 2d0*ERI(j,c,a,l))*Kronecker_delta(i,k)
                   
              end do
            end do
          end do
 
        end do
      end do
    end do
 
    !-------------!
    ! Block A2p1h !
    !-------------!
    
    iab = 0
    do i=nC+1,nO
      do a=nO+1,nOrb-nR
        do b=nO+1,nOrb-nR
          iab = iab + 1
             
          kcd = 0
          do k=nC+1,nO
            do c=nO+1,nOrb-nR
              do d=nO+1,nOrb-nR
                kcd = kcd + 1
                   
                H(1+n2h1p+iab,1+n2h1p+kcd) &
                     = ((eHF(a) + eHF(b) - eHF(i))*Kronecker_delta(i,k)*Kronecker_delta(a,c) & 
                     + 2d0*ERI(a,k,i,c) + 2d0*ERI(a,c,i,k))*Kronecker_delta(b,d)
!               H(1+2*n2h1p+n2p1h+iab,1+2*n2h1p+n2p1h+kcd) &
!                    = ((eHF(a) + eHF(b) - eHF(i))*Kronecker_delta(i,k)*Kronecker_delta(a,c) & 
!                    + 2d0*ERI(a,k,i,c))*Kronecker_delta(b,d)
                   
              end do
            end do
          end do
      
        end do
      end do
    end do
     
    !-------------!
    ! Block B2h1p !
    !-------------!
 
!   ija = 0
!   do i=nC+1,nO
!     do j=nC+1,nO
!       do a=nO+1,nOrb-nR
!         ija = ija + 1

!         kcd = 0
!         do k=nC+1,nO
!           do c=nO+1,nOrb-nR
!             do d=nO+1,nOrb-nR
!               kcd = kcd + 1
!                  
!               H(1+ija,1+n2h1p+kcd) = - 2d0*ERI(j,k,a,c)
!                  
!             end do
!           end do
!         end do
!
!       end do
!     end do
!   end do
 
    !-------------!
    ! Block B2p1h !
    !-------------!
    
!   iab = 0
!   do i=nC+1,nO
!     do a=nO+1,nOrb-nR
!       do b=nO+1,nOrb-nR
!         iab = iab + 1

!         klc = 0
!         do k=nC+1,nO
!           do l=nC+1,nO
!             do c=nO+1,nOrb-nR
!               klc = klc + 1            

!               H(1+n2h1p+iab,1+klc) = - 2d0*ERI(a,c,i,l)
!                  
!             end do
!           end do
!         end do
!     
!       end do
!     end do
!   end do
     
    !-------------------------!
    ! Diagonalize supermatrix !
    !-------------------------!

    call wall_time(start_timing)

    call diagonalize_general_matrix(nH,H,eGW,cGW)

    do s=1,nH
      order(s) = s
    end do

    call quick_sort(eGW,order,nH)
    call set_order(cGW,order,nH,nH)

    call wall_time(end_timing)
    timing = end_timing - start_timing

    write(*,*)
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for construction of supermatrix = ',timing,' seconds'
    write(*,*)

    !-----------------!
    ! Compute weights !
    !-----------------!

    do s=1,nH
      Z(s) = cGW(1,s)**2
    end do
 
    write(*,*)'-------------------------------------------'
    write(*,'(1X,A32,I3,A8)')'| G0W0 energies (eV) for orbital',p,'      |'
    write(*,*)'-------------------------------------------'
    write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X)') &
              '|','#','|','e_QP','|','Z','|'
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
 
  end do ! Loop on the orbital in the e block

end subroutine 
