subroutine R_IPEA_ADC3(dotest,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,EGHF,ERI,eHF)

! Dyson version of ADC(3)

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

  integer                       :: nI_2h1p,nII_2h1p,nIII_2h1p
  integer                       :: nI_2p1h,nII_2p1h,nIII_2p1h
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

  double precision,allocatable  :: Reigv(:,:) 

  double precision              :: start_timing,end_timing,timing

! Output variables

! Hello world

  write(*,*)
  write(*,*)'***************************************'
  write(*,*)'* Restricted IP/EA-ADC(3) Calculation *'
  write(*,*)'***************************************'
  write(*,*)

! Dimension of the supermatrix

  nI_2h1p   = nO*nV
  nII_2h1p  = nO*(nO-1)/2*nV
  nIII_2h1p = nO*(nO-1)/2*nV

  nI_2p1h   = nO*nV
  nII_2p1h  = nV*(nV-1)/2*nO
  nIII_2p1h = nV*(nV-1)/2*nO

  n2h1p = nI_2h1p + nII_2h1p + nIII_2h1p
  n2p1h = nI_2p1h + nII_2p1h + nIII_2p1h
  nH = nOrb + n2h1p + n2p1h

! Memory allocation

  allocate(H(nH,nH),eGF(nH),Z(nH),Reigv(nH,nH))

  eF = 0.5d0*(eHF(nO) + eHF(nO+1))

  H(:,:) = 0d0
  Reigv(:,:) = 0d0

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

  !----------------------!
  ! Block U_2h1p: Part I !
  !----------------------!

  do p=1,nOrb
     
    ija = nOrb
    do i=nC+1,nO
      do a=nO+1,nOrb-nR
        ija = ija + 1

        ! First-order contribution          

        H(p,ija) = ERI(i,i,p,a)

        ! Second-order contribution        

        do c=nO+1,nOrb-nR
          do d=nO+1,nOrb-nR

            H(p,ija) = H(p,ija) & 
                     - 0.5d0*ERI(i,i,c,d)*(ERI(c,d,p,a) + ERI(c,d,a,p))/(eHF(c) + eHF(d) - eHF(i) - eHF(i))

          end do
        end do
        
        do k=nC+1,nO
          do c=nO+1,nOrb-nR

            H(p,ija) = H(p,ija) & 
                     - (ERI(i,k,c,a) + ERI(i,k,a,c))*(ERI(c,i,p,k) + ERI(c,i,k,p))/(eHF(i) + eHF(k) - eHF(a) - eHF(c))

            H(p,ija) = H(p,ija) & 
                     + 3d0*ERI(i,k,a,c)*ERI(c,i,k,p)/(eHF(i) + eHF(k) - eHF(a) - eHF(c))

          end do
        end do
        
        ! Symmetrize

        H(ija,p) = H(p,ija)
              
      end do
    end do
   
  end do ! p

  !-----------------------!
  ! Block U_2h1p: Part II !
  !-----------------------!

  do p=1,nOrb
     
    ija = nOrb + nI_2h1p
    do i=nC+1,nO
      do j=i+1,nO
        do a=nO+1,nOrb-nR
          ija = ija + 1

          ! First-order contribution          

          H(p,ija) = sqrt(1.5d0)*(ERI(i,j,p,a) - ERI(i,j,a,p))

          ! Second-order contribution        

          do c=nO+1,nOrb-nR
            do d=nO+1,nOrb-nR

              H(p,ija) = H(p,ija) & 
                       - sqrt(1.5d0)*ERI(i,j,c,d)*(ERI(c,d,p,a) - ERI(c,d,a,p))/(eHF(c) + eHF(d) - eHF(i) - eHF(j))

            end do
          end do
         
          do k=nC+1,nO
            do c=nO+1,nOrb-nR

              H(p,ija) = H(p,ija) & 
                       - sqrt(1.5d0)*(ERI(i,k,c,a) - ERI(i,k,a,c))*(ERI(c,j,p,k) - ERI(c,j,k,p))/(eHF(i) + eHF(k) - eHF(a) - eHF(c))

              H(p,ija) = H(p,ija) & 
                       - sqrt(1.5d0)*ERI(i,k,a,c)*ERI(c,j,k,p)/(eHF(i) + eHF(k) - eHF(a) - eHF(c))

              H(p,ija) = H(p,ija) & 
                       + sqrt(1.5d0)*(ERI(j,k,c,a) - ERI(j,k,a,c))*(ERI(c,i,p,k) - ERI(c,i,k,p))/(eHF(j) + eHF(k) - eHF(a) - eHF(c))

              H(p,ija) = H(p,ija) & 
                       + sqrt(1.5d0)*ERI(j,k,a,c)*ERI(c,i,k,p)/(eHF(j) + eHF(k) - eHF(a) - eHF(c))

            end do
          end do
         
          ! Symmetrize

          H(ija,p) = H(p,ija)
              
        end do
      end do
    end do
   
  end do ! p

  !------------------------!
  ! Block U_2h1p: Part III !
  !------------------------!

  do p=1,nOrb
     
    ija = nOrb + nI_2h1p + nII_2h1p
    do i=nC+1,nO
      do j=i+1,nO
        do a=nO+1,nOrb-nR
          ija = ija + 1

          ! First-order contribution          

          H(p,ija) = sqrt(0.5d0)*(ERI(i,j,p,a) + ERI(i,j,a,p))

          ! Second-order contribution        

          do c=nO+1,nOrb-nR
            do d=nO+1,nOrb-nR

              H(p,ija) = H(p,ija) & 
                       - sqrt(0.5d0)*ERI(i,j,c,d)*(ERI(c,d,p,a) + ERI(c,d,a,p))/(eHF(c) + eHF(d) - eHF(i) - eHF(j))

            end do
          end do
         
          do k=nC+1,nO
            do c=nO+1,nOrb-nR

              H(p,ija) = H(p,ija) & 
                       - sqrt(0.5d0)*(ERI(j,k,c,a) + ERI(j,k,a,c))*(ERI(c,i,p,k) + ERI(c,i,k,p))/(eHF(j) + eHF(k) - eHF(a) - eHF(c))

              H(p,ija) = H(p,ija) & 
                       + 3d0*sqrt(0.5d0)*ERI(j,k,a,c)*ERI(c,i,k,p)/(eHF(j) + eHF(k) - eHF(a) - eHF(c))

              H(p,ija) = H(p,ija) & 
                       - sqrt(0.5d0)*(ERI(i,k,c,a) + ERI(i,k,a,c))*(ERI(c,j,p,k) + ERI(c,j,k,p))/(eHF(i) + eHF(k) - eHF(a) - eHF(c))

              H(p,ija) = H(p,ija) & 
                       + 3d0*sqrt(0.5d0)*ERI(i,k,a,c)*ERI(c,j,k,p)/(eHF(i) + eHF(k) - eHF(a) - eHF(c))

            end do
          end do
         
          ! Symmetrize

          H(ija,p) = H(p,ija)
              
        end do
      end do
    end do
   
  end do ! p

  !----------------------!
  ! Block U_2p1h: Part I !
  !----------------------!

  do p=1,nOrb
     
    iab = nOrb + n2h1p
    do i=nC+1,nO
      do a=nO+1,nOrb-nR
        iab = iab + 1

        ! First-order contribution          

        H(p,iab) = ERI(a,a,p,i)

        ! Second-order contribution        

        do k=nC+1,nO
          do l=nC+1,nO

            H(p,iab) = H(p,iab) & 
                     + 0.5d0*ERI(a,a,k,l)*(ERI(k,l,p,i) + ERI(k,l,i,p))/(eHF(k) + eHF(l) - eHF(a) - eHF(a))

          end do
        end do
        
        do k=nC+1,nO
          do c=nO+1,nOrb-nR

            H(p,iab) = H(p,iab) &
                     + (ERI(a,c,k,i) + ERI(a,c,i,k))*(ERI(k,a,p,c) + ERI(k,a,c,p))/(eHF(a) + eHF(c) - eHF(i) - eHF(k))

            H(p,iab) = H(p,iab) & 
                     - 3d0*ERI(a,c,i,k)*ERI(k,a,c,p)/(eHF(a) + eHF(c) - eHF(i) - eHF(k))

          end do
        end do
        
        ! Symmetrize

        H(iab,p) = H(p,iab)
              
      end do
    end do
   
  end do ! p

  !-----------------------!
  ! Block U_2p1h: Part II !
  !-----------------------!

  do p=1,nOrb
     
    iab = nOrb + n2h1p + nI_2p1h
    do i=nC+1,nO
      do a=nO+1,nOrb-nR
        do b=a+1,nOrb-nR
          iab = iab + 1

          ! First-order contribution          

          H(p,iab) = sqrt(1.5d0)*(ERI(a,b,p,i) - ERI(a,b,i,p))

          ! Second-order contribution        

          do k=nC+1,nO
            do l=nC+1,nO

              H(p,iab) = H(p,iab) & 
                       + sqrt(1.5d0)*ERI(a,b,k,l)*(ERI(k,l,p,i) - ERI(k,l,i,p))/(eHF(k) + eHF(l) - eHF(a) - eHF(b))

            end do
          end do
         
          do k=nC+1,nO
            do c=nO+1,nOrb-nR

              H(p,iab) = H(p,iab) & 
                       + sqrt(1.5d0)*(ERI(a,c,k,i) - ERI(a,c,i,k))*(ERI(k,b,p,c) - ERI(k,b,c,p))/(eHF(a) + eHF(c) - eHF(i) - eHF(k))

              H(p,iab) = H(p,iab) & 
                       + sqrt(1.5d0)*ERI(a,c,i,k)*ERI(k,b,c,p)/(eHF(a) + eHF(c) - eHF(i) - eHF(k))

              H(p,iab) = H(p,iab) & 
                       - sqrt(1.5d0)*(ERI(b,c,k,i) - ERI(b,c,i,k))*(ERI(k,a,p,c) - ERI(k,a,c,p))/(eHF(b) + eHF(c) - eHF(i) - eHF(k))

              H(p,iab) = H(p,iab) &
                       - sqrt(1.5d0)*ERI(b,c,i,k)*ERI(k,a,c,p)/(eHF(b) + eHF(c) - eHF(i) - eHF(k))

            end do
          end do
         
          ! Symmetrize

          H(iab,p) = H(p,iab)
              
        end do
      end do
    end do
   
  end do ! p

  !------------------------!
  ! Block U_2p1h: Part III !
  !------------------------!

  do p=1,nOrb
     
    iab = nOrb + n2h1p + nI_2p1h + nII_2p1h
    do i=nC+1,nO
      do a=nO+1,nOrb-nR
        do b=a+1,nOrb-nR
          iab = iab + 1

          ! First-order contribution          

          H(p,iab) = sqrt(0.5d0)*(ERI(a,b,p,i) + ERI(a,b,i,p))

          ! Second-order contribution        

          do k=nC+1,nO
            do l=nC+1,nO

              H(p,iab) = H(p,iab) &
                       + sqrt(0.5d0)*ERI(a,b,k,l)*(ERI(k,l,p,i) + ERI(k,l,i,p))/(eHF(k) + eHF(l) - eHF(a) - eHF(b))

            end do
          end do
         
          do k=nC+1,nO
            do c=nO+1,nOrb-nR

              H(p,iab) = H(p,iab) & 
                       + sqrt(0.5d0)*(ERI(b,c,k,i) + ERI(b,c,i,k))*(ERI(k,a,p,c) + ERI(k,a,c,p))/(eHF(b) + eHF(c) - eHF(i) - eHF(k))

              H(p,iab) = H(p,iab) & 
                       - 3d0*sqrt(0.5d0)*ERI(b,c,i,k)*ERI(k,a,c,p)/(eHF(b) + eHF(c) - eHF(i) - eHF(k))

              H(p,iab) = H(p,iab) & 
                       + sqrt(0.5d0)*(ERI(a,c,k,i) + ERI(a,c,i,k))*(ERI(k,b,p,c) + ERI(k,b,c,p))/(eHF(a) + eHF(c) - eHF(i) - eHF(k))

              H(p,iab) = H(p,iab) & 
                       - 3d0*sqrt(0.5d0)*ERI(a,c,i,k)*ERI(k,b,c,p)/(eHF(a) + eHF(c) - eHF(i) - eHF(k))

            end do
          end do
         
          ! Symmetrize

          H(iab,p) = H(p,iab)
              
        end do
      end do
    end do
   
  end do ! p

 !--------------!
 ! Block K_2h1p !
 !--------- ----!

  ija = nOrb
  do i=nC+1,nO
    do a=nO+1,nOrb-nR
      ija = ija + 1
           
      H(ija,ija) = eHF(i) + eHF(i) - eHF(a)

    end do
  end do
 
  ija = nOrb + nI_2h1p
  do i=nC+1,nO
    do j=i+1,nO
      do a=nO+1,nOrb-nR
        ija = ija + 1
           
        H(ija,ija) = eHF(i) + eHF(j) - eHF(a)

      end do
    end do
  end do
 
  ija = nOrb + nI_2h1p + nII_2h1p
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
      iab = iab + 1
             
      H(iab,iab) = eHF(a) + eHF(a) - eHF(i)
      
    end do
  end do

  iab = nOrb + n2h1p + nI_2p1h
  do i=nC+1,nO
    do a=nO+1,nOrb-nR
      do b=a+1,nOrb-nR
        iab = iab + 1

        H(iab,iab) = eHF(a) + eHF(b) - eHF(i)
      
      end do
    end do
  end do

  iab = nOrb + n2h1p + nI_2p1h + nII_2p1h
  do i=nC+1,nO
    do a=nO+1,nOrb-nR
      do b=a+1,nOrb-nR
        iab = iab + 1
             
        H(iab,iab) = eHF(a) + eHF(b) - eHF(i)
      
      end do
    end do
  end do

  !------------------------!
  ! Block C_2h1p: Part I-I !
  !------------------------!

  ija = nOrb
  do i=nC+1,nO
    do a=nO+1,nOrb-nR
      ija = ija + 1

      klc = nOrb
      do k=nC+1,nO
        do c=nO+1,nOrb-nR
          klc = klc + 1
         
          H(ija,klc) = H(ija,klc) &
                     - kronecker_delta(a,c)*ERI(i,i,k,k) &
                     - kronecker_delta(i,k)*(ERI(c,i,k,a) - 2d0*ERI(c,i,a,k))

        end do
      end do
    end do
  end do
 
  !--------------------------!
  ! Block C_2h1p: Part II-II !
  !--------------------------!

  ija = nOrb + nI_2h1p
  do i=nC+1,nO
    do j=i+1,nO
      do a=nO+1,nOrb-nR
        ija = ija + 1
         
        klc = nOrb + nI_2h1p
        do k=nC+1,nO
          do l=k+1,nO
            do c=nO+1,nOrb-nR
              klc = klc + 1
         
              H(ija,klc) = H(ija,klc) &
                         - kronecker_delta(a,c)*(ERI(i,j,k,l) - ERI(i,j,l,k)) &
                         + kronecker_delta(j,l)*ERI(c,i,a,k) &
                         - kronecker_delta(j,k)*ERI(c,i,a,l) &
                         - 1.5d0*kronecker_delta(i,k)*ERI(c,j,l,a) &
                         + 1.5d0*kronecker_delta(i,l)*ERI(c,j,k,a) &
                         - kronecker_delta(i,l)*ERI(c,j,a,k) &
                         + kronecker_delta(i,k)*ERI(c,j,a,l) &
                         - 1.5d0*kronecker_delta(j,l)*ERI(c,i,k,a) &
                         + 1.5d0*kronecker_delta(j,k)*ERI(c,i,l,a)  

            end do
          end do
        end do
         
      end do
    end do
  end do


  !----------------------------!
  ! Block C_2h1p: Part III-III !
  !----------------------------!

  ija = nOrb + nI_2h1p + nII_2h1p
  do i=nC+1,nO
    do j=i+1,nO
      do a=nO+1,nOrb-nR
        ija = ija + 1
         
        klc = nOrb + nI_2h1p + nII_2h1p
        do k=nC+1,nO
          do l=k+1,nO
            do c=nO+1,nOrb-nR
              klc = klc + 1
         
              H(ija,klc) = H(ija,klc) &
                         - kronecker_delta(a,c)*ERI(i,j,k,l) &
                         + kronecker_delta(i,k)*ERI(c,j,a,l) &
                         + kronecker_delta(j,l)*ERI(c,i,a,k) &
                         - 0.5d0*kronecker_delta(i,k)*ERI(c,j,l,a) &
                         - 0.5d0*kronecker_delta(j,l)*ERI(c,i,k,a) &
                         - kronecker_delta(a,c)*ERI(i,j,l,k) &
                         + kronecker_delta(i,l)*ERI(c,j,a,k) &
                         + kronecker_delta(j,k)*ERI(c,i,a,l) &
                         - 0.5d0*kronecker_delta(i,l)*ERI(c,j,k,a) &
                         - 0.5d0*kronecker_delta(j,k)*ERI(c,i,l,a)  

            end do
          end do
        end do
         
      end do
    end do
  end do

  !--------------------------------!
  ! Block C_2h1p: Part I-II & II-I !
  !--------------------------------!

  ija = nOrb
  do i=nC+1,nO
    do a=nO+1,nOrb-nR
      ija = ija + 1
         
      klc = nOrb + nI_2h1p
      do k=nC+1,nO
        do l=k+1,nO
          do c=nO+1,nOrb-nR
            klc = klc + 1
       
            H(ija,klc) = H(ija,klc) &
                       - sqrt(1.5d0)*kronecker_delta(i,k)*ERI(c,i,l,a) &
                       + sqrt(1.5d0)*kronecker_delta(i,l)*ERI(c,i,k,a)  

          end do
        end do
      end do

      H(klc,ija) = H(ija,klc)
         
    end do
  end do

  !----------------------------------!
  ! Block C_2h1p: Part I-III & III-I !
  !----------------------------------!

  ija = nOrb
  do i=nC+1,nO
    do a=nO+1,nOrb-nR
      ija = ija + 1
         
      klc = nOrb + nI_2h1p + nII_2h1p
      do k=nC+1,nO
        do l=k+1,nO
          do c=nO+1,nOrb-nR
            klc = klc + 1
       
            H(ija,klc) = H(ija,klc) &
                       - sqrt(2d0)*kronecker_delta(a,c)*ERI(i,i,k,l) &
                       + sqrt(2d0)*kronecker_delta(i,k)*ERI(c,i,a,l) &
                       - 0.5d0*sqrt(2d0)*kronecker_delta(i,k)*ERI(c,i,l,a) &
                       + sqrt(2d0)*kronecker_delta(i,l)*ERI(c,i,a,k) &
                       - 0.5d0*sqrt(2d0)*kronecker_delta(i,l)*ERI(c,i,k,a)  

          end do
        end do
      end do

      H(klc,ija) = H(ija,klc)
         
    end do
  end do

  !---------------------------!
  ! Block C_2h1p: Part II-III !
  !---------------------------!

  ija = nOrb + nI_2h1p
  do i=nC+1,nO
    do j=i+1,nO
      do a=nO+1,nOrb-nR
        ija = ija + 1
         
        klc = nOrb + nI_2h1p + nII_2h1p
        do k=nC+1,nO
          do l=k+1,nO
            do c=nO+1,nOrb-nR
              klc = klc + 1
         
              H(ija,klc) = H(ija,klc) &
                         - 0.5d0*sqrt(3d0)*kronecker_delta(i,k)*ERI(c,j,l,a) &
                         + 0.5d0*sqrt(3d0)*kronecker_delta(j,l)*ERI(c,i,k,a) &
                         - 0.5d0*sqrt(3d0)*kronecker_delta(i,l)*ERI(c,j,k,a) &
                         + 0.5d0*sqrt(3d0)*kronecker_delta(j,k)*ERI(c,i,l,a)

            end do
          end do
        end do
         
      end do
    end do
  end do

  !---------------------------!
  ! Block C_2h1p: Part III-II !
  !---------------------------!

  ija = nOrb + nI_2h1p + nII_2h1p
  do i=nC+1,nO
    do j=i+1,nO
      do a=nO+1,nOrb-nR
        ija = ija + 1
         
        klc = nOrb + nI_2h1p
        do k=nC+1,nO
          do l=k+1,nO
            do c=nO+1,nOrb-nR
              klc = klc + 1
         
              H(ija,klc) = H(ija,klc) &
                         - 0.5d0*sqrt(3d0)*kronecker_delta(i,k)*ERI(c,j,l,a) &
                         + 0.5d0*sqrt(3d0)*kronecker_delta(j,l)*ERI(c,i,k,a) &
                         + 0.5d0*sqrt(3d0)*kronecker_delta(i,l)*ERI(c,j,k,a) &
                         - 0.5d0*sqrt(3d0)*kronecker_delta(j,k)*ERI(c,i,l,a)

            end do
          end do
        end do
         
      end do
    end do
  end do
 
  !------------------------!
  ! Block C_2p1h: Part I-I !
  !------------------------! 

  iab = nOrb + n2h1p
  do i=nC+1,nO
    do a=nO+1,nOrb-nR
      iab = iab + 1
           
      kcd = nOrb + n2h1p
      do k=nC+1,nO
        do c=nO+1,nOrb-nR
          kcd = kcd + 1
      
          H(iab,kcd) = H(iab,kcd) &
                     + kronecker_delta(i,k)*ERI(a,a,c,c) &
                     + kronecker_delta(a,c)*(ERI(k,a,c,i) - 2d0*ERI(k,a,i,c))

        end do
      end do
           
    end do
  end do
   
  !--------------------------!
  ! Block C_2p1h: Part II-II !
  !--------------------------! 

  iab = nOrb + n2h1p + nI_2p1h
  do i=nC+1,nO
    do a=nO+1,nOrb-nR
      do b=a+1,nOrb-nR
        iab = iab + 1
           
        kcd = nOrb + n2h1p + nI_2p1h
        do k=nC+1,nO
          do c=nO+1,nOrb-nR
            do d=c+1,nOrb-nR
              kcd = kcd + 1
        
              H(iab,kcd) = H(iab,kcd) &
                         + kronecker_delta(i,k)*(ERI(a,b,c,d) - ERI(a,b,d,c)) &
                         - kronecker_delta(b,d)*ERI(k,a,i,c) &
                         + kronecker_delta(b,c)*ERI(k,a,i,d) &
                         + 1.5d0*kronecker_delta(a,c)*ERI(k,b,d,i) &
                         - 1.5d0*kronecker_delta(a,d)*ERI(k,b,c,i) &
                         + kronecker_delta(a,d)*ERI(k,b,i,c) &
                         - kronecker_delta(a,c)*ERI(k,b,i,d) &
                         - 1.5d0*kronecker_delta(b,c)*ERI(k,a,d,i) &
                         + 1.5d0*kronecker_delta(b,d)*ERI(k,a,c,i)  

            end do
          end do
        end do
           
      end do
    end do
  end do
   
  !----------------------------!
  ! Block C_2p1h: Part III-III !
  !----------------------------! 

  iab = nOrb + n2h1p + nI_2p1h + nII_2p1h
  do i=nC+1,nO
    do a=nO+1,nOrb-nR
      do b=a+1,nOrb-nR
        iab = iab + 1
           
        kcd = nOrb + n2h1p + nI_2p1h + nII_2p1h
        do k=nC+1,nO
          do c=nO+1,nOrb-nR
            do d=c+1,nOrb-nR
              kcd = kcd + 1
        
              H(iab,kcd) = H(iab,kcd) &
                         + kronecker_delta(i,k)*ERI(a,b,c,d) &
                         - kronecker_delta(a,c)*ERI(k,b,i,d) &
                         - kronecker_delta(b,d)*ERI(k,a,i,c) &
                         + 0.5d0*kronecker_delta(a,c)*ERI(k,b,d,i) &
                         + 0.5d0*kronecker_delta(b,d)*ERI(k,a,c,i) &
                         + kronecker_delta(i,k)*ERI(a,b,d,c) &
                         - kronecker_delta(a,d)*ERI(k,b,i,c) &
                         - kronecker_delta(b,c)*ERI(k,a,i,d) &
                         + 0.5d0*kronecker_delta(a,d)*ERI(k,b,c,i) &
                         + 0.5d0*kronecker_delta(b,c)*ERI(k,a,d,i)  

            end do
          end do
        end do
           
      end do
    end do
  end do
   
  !-------------------------!
  ! Block C_2p1h: Part I-II !
  !-------------------------! 

  iab = nOrb + n2h1p 
  do i=nC+1,nO
    do a=nO+1,nOrb-nR
      iab = iab + 1
           
      kcd = nOrb + n2h1p + nI_2p1h 
      do k=nC+1,nO
        do c=nO+1,nOrb-nR
          do d=c+1,nOrb-nR
            kcd = kcd + 1
      
            H(iab,kcd) = H(iab,kcd) &
                       + sqrt(1.5d0)*kronecker_delta(a,c)*ERI(k,a,d,i) &
                       - sqrt(1.5d0)*kronecker_delta(a,d)*ERI(k,a,c,i) 

          end do
        end do

        H(kcd,iab) = H(iab,kcd)
           
      end do
    end do
  end do
   
  !--------------------------!
  ! Block C_2p1h: Part I-III !
  !--------------------------! 

  iab = nOrb + n2h1p 
  do i=nC+1,nO
    do a=nO+1,nOrb-nR
      iab = iab + 1
           
      kcd = nOrb + n2h1p + nI_2p1h + nII_2p1h
      do k=nC+1,nO
        do c=nO+1,nOrb-nR
          do d=c+1,nOrb-nR
            kcd = kcd + 1
      
            H(iab,kcd) = H(iab,kcd) &
                       + sqrt(2d0)*kronecker_delta(i,k)*ERI(a,a,c,d) &
                       - sqrt(2d0)*kronecker_delta(a,c)*ERI(k,a,i,d) &
                       + 0.5d0*sqrt(2d0)*kronecker_delta(a,c)*ERI(k,a,d,i) &
                       - sqrt(2d0)*kronecker_delta(a,d)*ERI(k,a,i,c) &
                       + 0.5d0*sqrt(2d0)*kronecker_delta(a,d)*ERI(k,a,c,i)  

          end do
        end do

        H(kcd,iab) = H(iab,kcd)
           
      end do
    end do
  end do

  !---------------------------!
  ! Block C_2p1h: Part II-III !
  !---------------------------! 

  iab = nOrb + n2h1p + nI_2p1h 
  do i=nC+1,nO
    do a=nO+1,nOrb-nR
      do b=a+1,nOrb-nR
        iab = iab + 1
           
        kcd = nOrb + n2h1p + nI_2p1h + nII_2p1h
        do k=nC+1,nO
          do c=nO+1,nOrb-nR
            do d=c+1,nOrb-nR
              kcd = kcd + 1
        
              H(iab,kcd) = H(iab,kcd) &
                         + 0.5d0*sqrt(3d0)*kronecker_delta(a,c)*ERI(k,b,d,i) &
                         - 0.5d0*sqrt(3d0)*kronecker_delta(b,d)*ERI(k,a,c,i) &
                         + 0.5d0*sqrt(3d0)*kronecker_delta(a,d)*ERI(k,b,c,i) &
                         - 0.5d0*sqrt(3d0)*kronecker_delta(b,c)*ERI(k,a,d,i)  

            end do
          end do
        end do
           
      end do
    end do
  end do

  !---------------------------!
  ! Block C_2p1h: Part III-II !
  !---------------------------! 

  iab = nOrb + n2h1p + nI_2p1h + nII_2p1h
  do i=nC+1,nO
    do a=nO+1,nOrb-nR
      do b=a+1,nOrb-nR
        iab = iab + 1
           
        kcd = nOrb + n2h1p + nI_2p1h 
        do k=nC+1,nO
          do c=nO+1,nOrb-nR
            do d=c+1,nOrb-nR
              kcd = kcd + 1
        
              H(iab,kcd) = H(iab,kcd) &
                         + 0.5d0*sqrt(3d0)*kronecker_delta(a,c)*ERI(k,b,d,i) &
                         - 0.5d0*sqrt(3d0)*kronecker_delta(b,d)*ERI(k,a,c,i) &
                         - 0.5d0*sqrt(3d0)*kronecker_delta(a,d)*ERI(k,b,c,i) &
                         + 0.5d0*sqrt(3d0)*kronecker_delta(b,c)*ERI(k,a,d,i)  

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

  Z(:) = 0d0
  do s=1,nH
    do p=nC+1,nOrb-nR
      Z(s) = Z(s) + Reigv(p,s)**2
    end do
  end do

  !--------------!
  ! Dump results !
  !--------------!

  write(*,*)'-------------------------------------------'
  write(*,'(1X,A43)')'| IPEA-ADC(3) energies for all orbitals   |'
  write(*,*)'-------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X)') &
            '|','#','|','e_QP (eV)','|','Z','|'
  write(*,*)'-------------------------------------------'

  do s=1,nH
!   if(eGF(s) < eF .and. eGF(s) > eF - window) then
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
        if(p > nO) & 
             write(*,'(1X,A16,I3,A7,1X,F15.6,1X,F15.6)') &
             '               (',p,')      ',Reigv(1,s),Reigv(1,s)**2,-eHF(p)*HaToeV
  
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
         
        iab = 0
        do i=nC+1,nO
          do a=nO+1,nOrb-nR
            do b=nO+1,nOrb-nR
              iab = iab + 1

              if(abs(Reigv(1+n2h1p+iab,s)) > cutoff2)           &
                   write(*,'(1X,A7,I3,A6,I3,A1,I3,A3,1X,F15.6,1X,F15.6,1X,F12.6)') &
                   '      (',i,') -> (',a,',',b,')  ',Reigv(1+n2h1p+iab,s),Reigv(1+n2h1p+iab,s)**2, & 
                                                      (eHF(a) + eHF(b) - eHF(i))*HaToeV
                
            end do
          end do
        end do

        write(*,*)'------------------------------------------------------------------------------'
        write(*,*)

      end if ! If state s should be print

    end do ! Loop on s
     
  end if ! If verbose
  
end subroutine 
