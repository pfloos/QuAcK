subroutine R_G3W2_self_energy_diag_alt(eta,nBas,nOrb,nC,nO,nV,nR,nS,eHF,Om,rho,ERI,EcGM,Sig,Z)

! Alternative form of the G3W2 self-energy

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: eHF(nOrb)
  double precision,intent(in)   :: Om(nS)
  double precision,intent(in)   :: rho(nOrb,nOrb,nS)
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)

! Local variables

  integer                       :: p,r,s
  integer                       :: i,j,k,l
  integer                       :: a,b,c,d
  integer                       :: jb,kc,ia,ja
  integer                       :: n2h1p,n2p1h
  integer                       :: mu,nu
  integer                       :: klc,kcd,ija,ijb,iab,jab
  double precision              :: w
  double precision              :: num,num1,num2
  double precision              :: dem,dem1,dem2,dem3
  double precision              :: reg,reg1,reg2,reg3

  double precision,allocatable  :: C1_2h1p_2h1p(:,:)
  double precision,allocatable  :: C1_2p1h_2p1h(:,:)
  double precision,allocatable  :: C1_2h1p_2p1h(:,:)
  double precision,allocatable  :: K_2h1p_2h1p(:,:)
  double precision,allocatable  :: K_2p1h_2p1h(:,:)
  double precision,allocatable  :: dK_2h1p_2h1p(:,:)
  double precision,allocatable  :: dK_2p1h_2p1h(:,:)
  double precision,allocatable  :: U1_2h1p(:)
  double precision,allocatable  :: U1_2p1h(:)
  double precision,allocatable  :: U2_2h1p(:)
  double precision,allocatable  :: U2_2p1h(:)
  double precision,allocatable  :: KxU1_2h1p(:)
  double precision,allocatable  :: KxU1_2p1h(:)
  double precision,allocatable  :: dKxU1_2h1p(:)
  double precision,allocatable  :: dKxU1_2p1h(:)
  double precision,allocatable  :: U1xK_2h1p(:)
  double precision,allocatable  :: U1xK_2p1h(:)
  double precision,allocatable  :: U1xdK_2h1p(:)
  double precision,allocatable  :: U1xdK_2p1h(:)


  double precision              :: flow = 1d6

! Output variables

  double precision,intent(out)  :: Sig(nOrb)
  double precision,intent(out)  :: Z(nOrb)
  double precision,intent(out)  :: EcGM

! Dimension of the 2h1p and 2p1h subspaces

  n2h1p = nO*nO*nV
  n2p1h = nV*nV*nO

! Initialization

  Sig(:) = 0d0
  Z(:)   = 0d0

! Memory allocation

  allocate(C1_2h1p_2h1p(n2h1p,n2h1p))
  allocate(C1_2p1h_2p1h(n2p1h,n2p1h))
  allocate(C1_2h1p_2p1h(n2h1p,n2p1h))
  
  allocate(K_2h1p_2h1p(n2h1p,n2h1p))
  allocate(K_2p1h_2p1h(n2p1h,n2p1h))
  
  allocate(dK_2h1p_2h1p(n2h1p,n2h1p))
  allocate(dK_2p1h_2p1h(n2p1h,n2p1h))
  
  allocate(U1_2h1p(n2h1p))
  allocate(U1_2p1h(n2p1h))
  allocate(U2_2h1p(n2h1p))
  allocate(U2_2p1h(n2p1h))
  
!--------------------!
! Block C1_2h1p-2h1p !
!--------------------!
  
  C1_2h1p_2h1p(:,:) = 0d0

  ija = 0
  do i=nC+1,nO
    do mu=1,nS
      ija = ija + 1
  
      ! First-order terms
 
      klc = 0
      do k=nC+1,nO
        do nu=1,nS
          klc = klc + 1
     
          do r=nC+1,nOrb-nR

            num = rho(k,r,mu)*rho(i,r,nu)
            dem = eHF(i) - eHF(r) + Om(nu)
            reg = (1d0 - exp(-2d0*flow*dem*dem))/dem
           
            C1_2h1p_2h1p(ija,klc) = C1_2h1p_2h1p(ija,klc) + num*reg
           
            num = rho(k,r,mu)*rho(i,r,nu)
            dem = eHF(k) - eHF(r) + Om(mu)
            reg = (1d0 - exp(-2d0*flow*dem*dem))/dem
           
            C1_2h1p_2h1p(ija,klc) = C1_2h1p_2h1p(ija,klc) + num*reg

          end do
  
        end do
      end do
  
    end do
  end do
 
!--------------------!
! Block C1_2p1h-2p1h !
!--------------------!

  C1_2p1h_2p1h(:,:) = 0d0

  iab = 0
  do a=nO+1,nOrb-nR
    do mu=1,nS
      iab = iab + 1
 
      ! First-order terms
 
      kcd = 0
      do c=nO+1,nOrb-nR
        do nu=1,nS
          kcd = kcd + 1
     
          do r=nC+1,nOrb-nR

            num = rho(r,c,mu)*rho(r,a,nu)
            dem = eHF(c) - eHF(r) - Om(mu)
            reg = (1d0 - exp(-2d0*flow*dem*dem))/dem
           
            C1_2p1h_2p1h(iab,kcd) = C1_2p1h_2p1h(iab,kcd) + num*reg
           
            num = rho(r,c,mu)*rho(r,a,nu)
            dem = eHF(a) - eHF(r) - Om(nu)
            reg = (1d0 - exp(-2d0*flow*dem*dem))/dem
           
            C1_2p1h_2p1h(iab,kcd) = C1_2p1h_2p1h(iab,kcd) + num*reg

          end do
  
        end do
      end do
  
    end do
  end do
  
!--------------------!
! Block C1_2h1p-2p1h !
!--------------------!

  C1_2h1p_2p1h(:,:) = 0d0

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

            num = 2d0*rho(k,i,mu)*rho(a,k,nu)
            dem = eHF(a) - eHF(k) + Om(nu)
            reg = (1d0 - exp(-2d0*flow*dem*dem))/dem
           
            C1_2h1p_2p1h(ija,kcd) = C1_2h1p_2p1h(ija,kcd) + num*reg

          end do
         
          do c=nO+1,nOrb-nR

            num = 2d0*rho(c,i,mu)*rho(a,c,nu)
            dem = eHF(i) - eHF(c) - Om(mu)
            reg = (1d0 - exp(-2d0*flow*dem*dem))/dem
          
            C1_2h1p_2p1h(ija,kcd) = C1_2h1p_2p1h(ija,kcd) + num*reg
         
          end do
  
        end do
      end do
  
    end do
  end do
  
!-------------------------!
! Main loop over orbitals !
!-------------------------!

  do p=nC+1,nOrb-nR

    w = eHF(p)

    ! Downfolding the 3h2p configurations

    do i=nC+1,nO
      do mu=1,nS
      do nu=1,nS
          do r=nC+1,nOrb-nR
          do s=nC+1,nOrb-nR

            num1 = 2d0*rho(r,i,mu)*rho(p,r,nu)
            num2 = 2d0*rho(s,i,mu)*rho(p,s,nu)
            dem1 = eHF(r) - eHF(i) + Om(mu)
            dem2 = w - eHF(i) + Om(nu) + Om(mu)
            dem3 = eHF(s) - eHF(i) + Om(mu)

            reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))/dem1
            reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))/dem2
            reg3 = (1d0 - exp(-2d0*flow*dem3*dem3))/dem3

            Sig(p) = Sig(p) + num1*num2*reg1*reg2*reg3
            Z(p)   = Z(p)   - num1*num2*reg1*reg2*reg3/dem2

         end do
         end do
       end do
       end do
     end do

    ! Downfolding the 3p2h configurations

    do a=nO+1,nOrb-nR
      do mu=1,nS
      do nu=1,nS
          do r=nC+1,nOrb-nR
          do s=nC+1,nOrb-nR

            num1 = 2d0*rho(a,r,mu)*rho(r,p,nu)
            num2 = 2d0*rho(a,s,mu)*rho(s,p,nu)
            dem1 = eHF(r) - eHF(a) - Om(mu)
            dem2 = w - eHF(a) - Om(nu) - Om(mu)
            dem3 = eHF(s) - eHF(a) - Om(mu)

            reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))/dem1
            reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))/dem2
            reg3 = (1d0 - exp(-2d0*flow*dem3*dem3))/dem3

            Sig(p) = Sig(p) + num1*num2*reg1*reg2*reg3
            Z(p)   = Z(p)   - num1*num2*reg1*reg2*reg3/dem2

         end do
         end do
       end do
       end do
     end do

    !---------------!
    ! Blocks U_2h1p !
    !---------------!

     U2_2h1p(:) = 0d0
     U1_2h1p(:) = 0d0

    ija = 0
    do i=nC+1,nO
      do mu=1,nS
        ija = ija + 1

        ! First-order terms

        U1_2h1p(ija) = sqrt(2d0)*rho(p,i,mu)

        ! Second-order terms

        do k=nC+1,nO
          do c=nO+1,nOrb-nR

          num = sqrt(2d0)*rho(k,c,mu)*ERI(i,k,c,p)
          dem = eHF(c) - eHF(k) - Om(mu)
          reg = (1d0 - exp(-2d0*flow*dem*dem))/dem

          U2_2h1p(ija) = U2_2h1p(ija) + num*reg

          num = sqrt(2d0)*rho(c,k,mu)*ERI(i,c,k,p)
          dem = eHF(c) - eHF(k) + Om(mu)
          reg = (1d0 - exp(-2d0*flow*dem*dem))/dem

          U2_2h1p(ija) = U2_2h1p(ija) + num*reg

          end do
        end do

        ! Third-order terms

        do k=nC+1,nO
          do c=nO+1,nOrb-nR
            do nu=1,nS

              num = 2d0*sqrt(2d0)*rho(c,k,mu)*rho(i,k,nu)*rho(p,c,nu)
              dem1 = eHF(c) - eHF(k) + Om(mu)
              dem2 = eHF(k) - eHF(i) - Om(nu)

              reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))/dem1
              reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))/dem2

              U2_2h1p(ija) = U2_2h1p(ija) + num*reg1*reg2

              num = 2d0*sqrt(2d0)*rho(k,c,mu)*rho(c,i,nu)*rho(k,p,nu)
              dem1 = eHF(k) - eHF(c) + Om(mu)
              dem2 = eHF(i) - eHF(c) - Om(nu)

              reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))/dem1
              reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))/dem2

              U2_2h1p(ija) = U2_2h1p(ija) - num*reg1*reg2

              num = 2d0*sqrt(2d0)*rho(k,c,mu)*rho(i,c,nu)*rho(p,k,nu)
              dem1 = eHF(k) - eHF(c) + Om(mu)
              dem2 = eHF(c) - eHF(i) - Om(nu)

              reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))/dem1
              reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))/dem2

              U2_2h1p(ija) = U2_2h1p(ija) - 0.5d0*num*reg1*reg2

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

              U2_2h1p(ija) = U2_2h1p(ija) + 0.5d0*num*reg1*reg2

            end do
          end do
        end do

      end do
    end do

    !---------------!
    ! Blocks U_2p1h !
    !---------------!

    U2_2p1h(:) = 0d0
    U1_2p1h(:) = 0d0

    iab = 0
    do a=nO+1,nOrb-nR
      do mu=1,nS
        iab = iab + 1

        ! First-order terms

        U1_2p1h(iab) = sqrt(2d0)*rho(a,p,mu)

        ! Second-order terms

        do k=nC+1,nO
          do c=nO+1,nOrb-nR

          num = sqrt(2d0)*rho(k,c,mu)*ERI(a,c,k,p)
          dem = eHF(c) - eHF(k) - Om(mu)
          reg = (1d0 - exp(-2d0*flow*dem*dem))/dem

          U2_2p1h(iab) = U2_2p1h(iab) + num*reg

          num = sqrt(2d0)*rho(c,k,mu)*ERI(a,k,c,p)
          dem = eHF(c) - eHF(k) + Om(mu)
          reg = (1d0 - exp(-2d0*flow*dem*dem))/dem

          U2_2p1h(iab) = U2_2p1h(iab) + num*reg

          end do
        end do

        ! Third-order terms

        do k=nC+1,nO
          do c=nO+1,nOrb-nR
            do nu=1,nS

              num = 2d0*sqrt(2d0)*rho(c,k,mu)*rho(c,a,nu)*rho(k,p,nu)
              dem1 = eHF(c) - eHF(k) + Om(mu)
              dem2 = eHF(a) - eHF(c) - Om(nu)

              reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))/dem1
              reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))/dem2

              U2_2p1h(iab) = U2_2p1h(iab) + num*reg1*reg2

              num = 2d0*sqrt(2d0)*rho(k,c,mu)*rho(a,k,nu)*rho(p,c,nu)
              dem1 = eHF(k) - eHF(c) + Om(mu)
              dem2 = eHF(k) - eHF(a) - Om(nu)

              reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))/dem1
              reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))/dem2

              U2_2p1h(iab) = U2_2p1h(iab) - num*reg1*reg2

              num = 2d0*sqrt(2d0)*rho(k,c,mu)*rho(k,a,nu)*rho(c,p,nu)
              dem1 = eHF(k) - eHF(c) + Om(mu)
              dem2 = eHF(a) - eHF(k) - Om(nu)

              reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))/dem1
              reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))/dem2

              U2_2p1h(iab) = U2_2p1h(iab) - 0.5d0*num*reg1*reg2

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

              U2_2p1h(iab) = U2_2p1h(iab) + 0.5d0*num*reg1*reg2

            end do
          end do
        end do

      end do
    end do

    !-------------------!
    ! Block K_2h1p-2h1p !
    !-------------------!

    K_2h1p_2h1p(:,:) = 0d0

    ija = 0
    do i=nC+1,nO
      do mu=1,nS
        ija = ija + 1

        ! Zeroth-order terms
   
        K_2h1p_2h1p(ija,ija) = 1d0/(w - eHF(i) + Om(mu))

      end do
    end do

    !-------------------!
    ! Block K_2p1h-2p1h !
    !-------------------!

    K_2p1h_2p1h(:,:) = 0d0

    iab = 0
    do a=nO+1,nOrb-nR
      do mu=1,nS
        iab = iab + 1

        ! Zeroth-order terms

        K_2p1h_2p1h(iab,iab) = 1d0/(w - eHF(a) - Om(mu))

      end do
    end do

    allocate(U1xK_2h1p(n2h1p),U1xK_2p1h(n2p1h),KxU1_2h1p(n2h1p),KxU1_2p1h(n2p1h))

    U1xK_2h1p = matmul(U1_2h1p,K_2h1p_2h1p)
    U1xK_2p1h = matmul(U1_2p1h,K_2p1h_2p1h)
    KxU1_2h1p = matmul(K_2h1p_2h1p,U1_2h1p)
    KxU1_2p1h = matmul(K_2p1h_2p1h,U1_2p1h)

    Sig(p) = Sig(p) &
 
           + dot_product(U1xK_2h1p,U1_2h1p) &
           + dot_product(U1xK_2p1h,U1_2p1h) &
 
           + dot_product(U1xK_2h1p,U2_2h1p) &
           + dot_product(U1xK_2p1h,U2_2p1h) &
 
           + dot_product(U2_2h1p,KxU1_2h1p) &
           + dot_product(U2_2p1h,KxU1_2p1h) &
 
           + dot_product(U1xK_2p1h,matmul(C1_2p1h_2p1h,KxU1_2p1h)) &
           + dot_product(U1xK_2h1p,matmul(C1_2h1p_2h1p,KxU1_2h1p)) &
            
           + dot_product(U1xK_2h1p,matmul(C1_2h1p_2p1h,KxU1_2p1h)) &
           + dot_product(U1xK_2p1h,matmul(transpose(C1_2h1p_2p1h),KxU1_2h1p))

    !-------------------!
    ! Block K_2h1p-2h1p !
    !-------------------!

    dK_2h1p_2h1p(:,:) = 0d0

    ija = 0
    do i=nC+1,nO
      do mu=1,nS
        ija = ija + 1

        ! Zeroth-order terms
   
        dK_2h1p_2h1p(ija,ija) = - 1d0/(w - eHF(i) + Om(mu))**2

      end do
    end do

    !-------------------!
    ! Block K_2p1h-2p1h !
    !-------------------!

    dK_2p1h_2p1h(:,:) = 0d0

    iab = 0
    do a=nO+1,nOrb-nR
      do mu=1,nS
        iab = iab + 1

        ! Zeroth-order terms

        dK_2p1h_2p1h(iab,iab) = - 1d0/(w - eHF(a) - Om(mu))**2

      end do
    end do

    allocate(U1xdK_2h1p(n2h1p),U1xdK_2p1h(n2p1h),dKxU1_2h1p(n2h1p),dKxU1_2p1h(n2p1h))

    U1xdK_2h1p = matmul(U1_2h1p,dK_2h1p_2h1p)
    U1xdK_2p1h = matmul(U1_2p1h,dK_2p1h_2p1h)
    dKxU1_2h1p = matmul(dK_2h1p_2h1p,U1_2h1p)
    dKxU1_2p1h = matmul(dK_2p1h_2p1h,U1_2p1h)

    Z(p) = Z(p) &
 
         + dot_product(U1xdK_2h1p,U1_2h1p) &
         + dot_product(U1xdK_2p1h,U1_2p1h) &
 
         + dot_product(U1xdK_2h1p,U2_2h1p) &
         + dot_product(U1xdK_2p1h,U2_2p1h) &
 
         + dot_product(U2_2h1p,dKxU1_2h1p) &
         + dot_product(U2_2p1h,dKxU1_2p1h) &
         
         + dot_product(U1xdK_2p1h,matmul(C1_2p1h_2p1h,KxU1_2p1h)) &
         + dot_product(U1xdK_2h1p,matmul(C1_2h1p_2h1p,KxU1_2h1p)) &
          
         + dot_product(U1xdK_2h1p,matmul(C1_2h1p_2p1h,KxU1_2p1h)) &
         + dot_product(U1xdK_2p1h,matmul(transpose(C1_2h1p_2p1h),KxU1_2h1p)) &

         + dot_product(U1xK_2p1h,matmul(C1_2p1h_2p1h,dKxU1_2p1h)) &
         + dot_product(U1xK_2h1p,matmul(C1_2h1p_2h1p,dKxU1_2h1p)) &
          
         + dot_product(U1xK_2h1p,matmul(C1_2h1p_2p1h,dKxU1_2p1h)) &
         + dot_product(U1xK_2p1h,matmul(transpose(C1_2h1p_2p1h),dKxU1_2h1p))

    deallocate(U1xK_2h1p,U1xK_2p1h,KxU1_2h1p,KxU1_2p1h)
    deallocate(U1xdK_2h1p,U1xdK_2p1h,dKxU1_2h1p,dKxU1_2p1h)
 
  end do

  print*,'Alternative form of the self-energy'
  call vecout(nOrb,Sig)

  Z(:) = 1d0/(1d0 - Z(:))
  call vecout(nOrb,Z)

end subroutine 
