subroutine R_psdG3W2_self_energy_diag(eta,flow,nBas,nOrb,nC,nO,nV,nR,nS,eHF,Om,rho,ERI,EcGM,SigC,Z)

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
  double precision,intent(in)   :: flow

! Local variables

  integer                       :: p,r
  integer                       :: i,j,k,l
  integer                       :: a,b,c,d
  integer                       :: jb,kc,ia,ja
  integer                       :: inu,jmu,anu,bmu
  integer                       :: n2h1p,n2p1h
  integer                       :: mu,nu,s,t
  integer                       :: klc,kcd,ija,ijb,iab,jab
  double precision              :: w
  double precision              :: num,num1,num2,num3
  double precision              :: dem,dem1,dem2,dem3
  double precision              :: reg,reg1,reg2,reg3
  double precision              :: tmpA,tmpB,tmpBZ

  double precision,allocatable  :: C1_2h1p(:,:)
  double precision,allocatable  :: C1_2p1h(:,:)
  double precision,allocatable  :: K_2h1p(:)
  double precision,allocatable  :: K_2p1h(:)
  double precision,allocatable  :: dK_2h1p(:)
  double precision,allocatable  :: dK_2p1h(:)
  double precision,allocatable  :: U1_2h1p(:)
  double precision,allocatable  :: U1_2p1h(:)
  double precision,allocatable  :: U2_2h1p(:)
  double precision,allocatable  :: U2_2p1h(:)
  double precision,allocatable  :: U3_2h1p(:)
  double precision,allocatable  :: U3_2p1h(:)



! Output variables

  double precision,intent(out)  :: SigC(nOrb)
  double precision,intent(out)  :: Z(nOrb)
  double precision,intent(out)  :: EcGM

! Dimension of the 2h1p and 2p1h subspaces

  n2h1p = nO*nO*nV
  n2p1h = nV*nV*nO

! Initialization

  SigC(:) = 0d0
  Z(:)   = 0d0
  EcGM = 0d0

! Memory allocation

  allocate(C1_2h1p(n2h1p,n2h1p))
  allocate(C1_2p1h(n2p1h,n2p1h))
  
  allocate(K_2h1p(n2h1p))
  allocate(K_2p1h(n2p1h))
  
  allocate(dK_2h1p(n2h1p))
  allocate(dK_2p1h(n2p1h))
  
  allocate(U1_2h1p(n2h1p))
  allocate(U1_2p1h(n2p1h))
  allocate(U2_2h1p(n2h1p))
  allocate(U2_2p1h(n2p1h))
  allocate(U3_2h1p(n2h1p))
  allocate(U3_2p1h(n2p1h))
  
!---------------!
! Block C1_2h1p !
!---------------!
  
  C1_2h1p(:,:) = 0d0

  inu = 0
  do i=nC+1,nO
    do nu=1,nS
      inu = inu + 1
  
      ! First-order terms
 
      jmu = 0
      do j=nC+1,nO
        do mu=1,nS
          jmu = jmu + 1
   
          do k=nC+1,nO

            num = rho(i,k,mu)*rho(j,k,nu)
            dem = eHF(i) - eHF(k) + Om(mu)
            reg = (1d0 - exp(-2d0*flow*dem*dem))
         
            C1_2h1p(inu,jmu) = C1_2h1p(inu,jmu) + num*reg/dem
         
            num = rho(i,k,mu)*rho(j,k,nu)
            dem = eHF(j) - eHF(k) + Om(nu)
            reg = (1d0 - exp(-2d0*flow*dem*dem))
         
            C1_2h1p(inu,jmu) = C1_2h1p(inu,jmu) + num*reg/dem

          end do
  
        end do
      end do
  
    end do
  end do
 
!---------------!
! Block C1_2p1h !
!---------------!

  C1_2p1h(:,:) = 0d0

  anu = 0
  do a=nO+1,nOrb-nR
    do nu=1,nS
      anu = anu + 1
 
      ! First-order terms
 
      bmu = 0
      do b=nO+1,nOrb-nR
        do mu=1,nS
          bmu = bmu + 1
   
          do c=nO+1,nOrb-nR

            num = rho(c,a,mu)*rho(c,b,nu)
            dem = eHF(a) - eHF(c) - Om(mu)
            reg = (1d0 - exp(-2d0*flow*dem*dem))
         
            C1_2p1h(anu,bmu) = C1_2p1h(anu,bmu) + num*reg/dem
         
            num = rho(c,a,mu)*rho(c,b,nu)
            dem = eHF(b) - eHF(c) - Om(nu)
            reg = (1d0 - exp(-2d0*flow*dem*dem))
         
            C1_2p1h(anu,bmu) = C1_2p1h(anu,bmu) + num*reg/dem

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
       do nu=1,nS
          do mu=1,nS
                   
             do j=nC+1,nO
                do k=nC+1,nO

                   num1 = 2d0*rho(k,i,nu)*rho(p,k,mu)
                   dem1 = w - eHF(k) + Om(mu)
                   reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))
                   
                   dem2 = w - eHF(i) + Om(nu) + Om(mu)
                   reg2 = (1d0 - exp(-2d0*flow*dem2*dem2)) 
                
                   num3 = 2d0*rho(j,i,mu)*rho(p,j,nu)
                   dem3 = w - eHF(j) + Om(nu)
                   reg3 = (1d0 - exp(-2d0*flow*dem3*dem3))
                   
                   SigC(p) = SigC(p) + num1*num3*reg1*reg2*reg3/dem1/dem2/dem3
                   
                   Z(p)    = Z(p)    - num1*num3*reg1*reg2*reg3/dem1/dem2/dem3/dem1
                   Z(p)    = Z(p)    - num1*num3*reg1*reg2*reg3/dem1/dem2/dem3/dem2
                   Z(p)    = Z(p)    - num1*num3*reg1*reg2*reg3/dem1/dem2/dem3/dem3  
                   
                end do
             end do

             do k=nC+1,nO
                do a=nO+1,nOrb-nR

                   num1 = 2d0*rho(k,i,nu)*rho(p,k,mu)
                   dem1 = w - eHF(k) + Om(mu)
                   reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))
                   
                   dem2 = w - eHF(i) + Om(nu) + Om(mu)
                   reg2 = (1d0 - exp(-2d0*flow*dem2*dem2)) 
                
                   num3 = 2d0*rho(a,i,mu)*rho(p,a,nu)
                   dem3 = eHF(i) - eHF(a) - Om(mu)
                   reg3 = (1d0 - exp(-2d0*flow*dem3*dem3))
                   
                   SigC(p) = SigC(p) + 2d0*num1*num3*reg1*reg2*reg3/dem1/dem2/dem3
                   
                   Z(p)    = Z(p)    - 2d0*num1*num3*reg1*reg2*reg3/dem1/dem2/dem3/dem1
                   Z(p)    = Z(p)    - 2d0*num1*num3*reg1*reg2*reg3/dem1/dem2/dem3/dem2
                   
                end do
             end do

             do b=nO+1,nOrb-nR
                do a=nO+1,nOrb-nR

                   num1 = 2d0*rho(b,i,nu)*rho(p,b,mu)
                   dem1 = eHF(i) - eHF(b) - Om(nu)
                   reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))
                   
                   dem2 = w - eHF(i) + Om(nu) + Om(mu)
                   reg2 = (1d0 - exp(-2d0*flow*dem2*dem2)) 
                
                   num3 = 2d0*rho(a,i,mu)*rho(p,a,nu)
                   dem3 = eHF(i) - eHF(a) - Om(mu)
                   reg3 = (1d0 - exp(-2d0*flow*dem3*dem3))
                   
                   SigC(p) = SigC(p) + num1*num3*reg1*reg2*reg3/dem1/dem2/dem3

                   Z(p)    = Z(p)    - num1*num3*reg1*reg2*reg3/dem1/dem2/dem3/dem2
                   
                end do
             end do

                
          end do
       end do
    end do

    ! Downfolding the 3p2h configurations

    do a=nO+1,nOrb-nR
       do nu=1,nS
          do mu=1,nS
             
             do b=nO+1,nOrb-nR
                do c=nO+1,nOrb-nR
                   
                   num1 = 2d0*rho(a,c,nu)*rho(c,p,mu)
                   dem1 = w - eHF(c) - Om(mu)
                   reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))
                   
                   dem2 = w - eHF(a) - Om(nu) - Om(mu)
                   reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))  
                   
                   num3 = 2d0*rho(a,b,mu)*rho(b,p,nu)
                   dem3 = w - eHF(b) - Om(nu)
                   reg3 = (1d0 - exp(-2d0*flow*dem3*dem3))


                   SigC(p) = SigC(p) + num1*num3*reg1*reg2*reg3/dem1/dem2/dem3
                   
                   Z(p)    = Z(p)    - num1*num3*reg1*reg2*reg3/dem1/dem2/dem3/dem1
                   Z(p)    = Z(p)    - num1*num3*reg1*reg2*reg3/dem1/dem2/dem3/dem2
                   Z(p)    = Z(p)    - num1*num3*reg1*reg2*reg3/dem1/dem2/dem3/dem3
                   
                end do
             end do

             do c=nO+1,nOrb-nR
                do i=nC+1,nO
                   
                   num1 = 2d0*rho(a,c,nu)*rho(c,p,mu)
                   dem1 = w - eHF(c) - Om(mu)
                   reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))
                   
                   dem2 = w - eHF(a) - Om(nu) - Om(mu)
                   reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))  
                   
                   num3 = 2d0*rho(a,i,mu)*rho(i,p,nu)
                   dem3 = eHF(a) - eHF(i) + Om(mu)
                   reg3 = (1d0 - exp(-2d0*flow*dem3*dem3))


                   SigC(p) = SigC(p) + 2d0*num1*num3*reg1*reg2*reg3/dem1/dem2/dem3
                   
                   Z(p)    = Z(p)    - 2d0*num1*num3*reg1*reg2*reg3/dem1/dem2/dem3/dem1
                   Z(p)    = Z(p)    - 2d0*num1*num3*reg1*reg2*reg3/dem1/dem2/dem3/dem2
                   
                end do
             end do

             do j=nC+1,nO
                do i=nC+1,nO
                   
                   num1 = 2d0*rho(a,j,nu)*rho(j,p,mu)
                   dem1 = eHF(a) - eHF(j) + Om(nu)
                   reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))
                   
                   dem2 = w - eHF(a) - Om(nu) - Om(mu)
                   reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))  
                   
                   num3 = 2d0*rho(a,i,mu)*rho(i,p,nu)
                   dem3 = eHF(a) - eHF(i) + Om(mu)
                   reg3 = (1d0 - exp(-2d0*flow*dem3*dem3))


                   SigC(p) = SigC(p) + num1*num3*reg1*reg2*reg3/dem1/dem2/dem3
                   
                   Z(p)    = Z(p)    - num1*num3*reg1*reg2*reg3/dem1/dem2/dem3/dem2
                   
                end do
             end do
             
             
          end do
       end do
    end do
     
     !---------------!
     ! Blocks U_2h1p !
     !---------------!
     
     U1_2h1p(:) = 0d0
     U2_2h1p(:) = 0d0
     U3_2h1p(:) = 0d0
     
     inu = 0
     do i=nC+1,nO
        do nu=1,nS
           inu = inu + 1
           
           ! First-order terms
           
           U1_2h1p(inu) = sqrt(2d0)*rho(p,i,nu)

           ! Second-order terms

           do k=nC+1,nO
              do c=nO+1,nOrb-nR
                 
                 num = sqrt(2d0)*rho(k,c,nu)*ERI(i,k,c,p)
                 dem = eHF(c) - eHF(k) - Om(nu)
                 reg = (1d0 - exp(-2d0*flow*dem*dem))
                 
                 U2_2h1p(inu) = U2_2h1p(inu) + num*reg/dem
          
                 num = sqrt(2d0)*rho(c,k,nu)*ERI(i,c,k,p)
                 dem = eHF(c) - eHF(k) + Om(nu)
                 reg = (1d0 - exp(-2d0*flow*dem*dem))

                 U2_2h1p(inu) = U2_2h1p(inu) + num*reg/dem
                 
              end do
           end do

           ! Third-order terms

           do k=nC+1,nO
              do c=nO+1,nOrb-nR
                 do mu=1,nS

                    num = 2d0*sqrt(2d0)*rho(k,c,nu)*rho(c,i,mu)*rho(k,p,mu)
                    dem1 = eHF(c) - eHF(k) - Om(nu)
                    dem2 = eHF(c) - eHF(i) + Om(mu)

                    reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))
                    reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))

                    U3_2h1p(inu) = U3_2h1p(inu) - num*reg1*reg2/dem1/dem2

                    num = 2d0*sqrt(2d0)*rho(k,c,nu)*rho(i,c,mu)*rho(p,k,mu)
                    dem1 = eHF(c) - eHF(k) - Om(nu)
                    dem2 = eHF(c) - eHF(i) - Om(mu)

                    reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))
                    reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))

                    U3_2h1p(inu) = U3_2h1p(inu) + 0.5d0*num*reg1*reg2/dem1/dem2

                    num = 2d0*sqrt(2d0)*rho(k,i,mu)*rho(c,k,nu)*rho(c,p,mu)
                    dem1 = eHF(c) - eHF(i) + Om(mu) + Om(nu)
                    dem2 = eHF(c) - eHF(k) + Om(nu)

                    reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))
                    reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))
                  
                    U3_2h1p(inu) = U3_2h1p(inu) - num*reg1*reg2/dem1/dem2

                 end do
              end do
           end do

           do a=nO+1,nOrb-nR
            do b=nO+1,nOrb-nR
              do mu=1,nS
    
                num = 2d0*sqrt(2d0)*rho(a,i,mu)*rho(b,a,nu)*rho(b,p,mu)
                dem1 = eHF(b) - eHF(i) + Om(mu) + Om(nu)
                dem2 = eHF(a) - eHF(i) + Om(mu)

                reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))
                reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))
                  
                U3_2h1p(inu) = U3_2h1p(inu) + num*reg1*reg2/dem1/dem2

 
              end do
            end do
          end do

        end do
     end do

     !---------------!
     ! Blocks U_2p1h !
     !---------------!

     U1_2p1h(:) = 0d0
     U2_2p1h(:) = 0d0
     U3_2p1h(:) = 0d0

     anu = 0
     do a=nO+1,nOrb-nR
        do nu=1,nS
           anu = anu + 1

           ! First-order terms

           U1_2p1h(anu) = sqrt(2d0)*rho(a,p,nu)

           ! Second-order terms

           do k=nC+1,nO
              do c=nO+1,nOrb-nR

                 num = sqrt(2d0)*rho(k,c,nu)*ERI(a,c,k,p)
                 dem = eHF(c) - eHF(k) - Om(nu)
                 reg = (1d0 - exp(-2d0*flow*dem*dem))

                 U2_2p1h(anu) = U2_2p1h(anu) + num*reg/dem

                 num = sqrt(2d0)*rho(c,k,nu)*ERI(a,k,c,p)
                 dem = eHF(c) - eHF(k) + Om(nu)
                 reg = (1d0 - exp(-2d0*flow*dem*dem))

                 U2_2p1h(anu) = U2_2p1h(anu) + num*reg/dem

              end do
           end do

           ! Third-order terms
           
           do k=nC+1,nO
              do c=nO+1,nOrb-nR
                 do mu=1,nS

                    num = 2d0*sqrt(2d0)*rho(k,c,nu)*rho(a,k,mu)*rho(p,c,mu)
                    dem1 = eHF(c) - eHF(k) - Om(nu)
                    dem2 = eHF(a) - eHF(k) + Om(mu)

                    reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))
                    reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))

                    U3_2p1h(anu) = U3_2p1h(anu) - num*reg1*reg2/dem1/dem2
                  
                    num = 2d0*sqrt(2d0)*rho(k,c,nu)*rho(k,a,mu)*rho(c,p,mu)
                    dem1 = eHF(c) - eHF(k) - Om(nu)
                    dem2 = eHF(a) - eHF(k) - Om(mu)

                    reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))
                    reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))

                    U3_2p1h(anu) = U3_2p1h(anu) + 0.5d0*num*reg1*reg2/dem1/dem2

                    num = 2d0*sqrt(2d0)*rho(a,c,mu)*rho(c,k,nu)*rho(p,k,mu)
                    dem1 = eHF(a) - eHF(k) + Om(mu) + Om(nu)
                    dem2 = eHF(c) - eHF(k) + Om(nu)

                    reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))
                    reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))

                    U3_2p1h(anu) = U3_2p1h(anu) - num*reg1*reg2/dem1/dem2

                 end do
              end do
           end do

           do i=nC+1,nO
              do j=nC+1,nO
                 do mu=1,nS
                    
                    num = 2d0*sqrt(2d0)*rho(a,j,mu)*rho(j,i,nu)*rho(p,i,mu)
                    dem1 = eHF(a) - eHF(i) + Om(mu) + Om(nu)
                    dem2 = eHF(a) - eHF(j) + Om(mu)

                    reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))
                    reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))
            
                    U3_2p1h(anu) = U3_2p1h(anu) + num*reg1*reg2/dem1/dem2

                 end do
              end do
           end do

        end do
    end do

    !--------------!
    ! Block K_2h1p !
    !--------------!

    K_2h1p(:)  = 0d0
    dK_2h1p(:) = 0d0

    inu = 0
    do i=nC+1,nO
       do nu=1,nS
          inu = inu + 1
          
          ! Zeroth-order terms
          dem = w - eHF(i) + Om(nu)
          reg = (1d0 - exp(-2d0*flow*dem*dem))
          K_2h1p(inu)  =   reg/dem
          dK_2h1p(inu) = - reg/dem/dem

      end do
    end do

    !--------------!
    ! Block K_2p1h !
    !--------------!

    K_2p1h(:)  = 0d0
    dK_2p1h(:) = 0d0

    anu = 0
    do a=nO+1,nOrb-nR
       do nu=1,nS
          anu = anu + 1
          
          ! Zeroth-order terms
          dem = w - eHF(a) - Om(nu)
          reg = (1d0 - exp(-2d0*flow*dem*dem))
          K_2p1h(anu)  =   reg/dem
          dK_2p1h(anu) = - reg/dem/dem

       end do
    end do

    !----------------------------!
    !    Building self-energy    !
    !----------------------------!
    inu = 0
    do i=nC+1,nO
       do nu=1,nS
          inu = inu + 1
             
             SigC(p) = SigC(p) + U1_2h1p(inu) *  K_2h1p(inu) * U1_2h1p(inu)
             Z(p)    = Z(p)    + U1_2h1p(inu) * dK_2h1p(inu) * U1_2h1p(inu)
             
             SigC(p) = SigC(p) + U2_2h1p(inu) *  K_2h1p(inu) * U1_2h1p(inu)
             Z(p)    = Z(p)    + U2_2h1p(inu) * dK_2h1p(inu) * U1_2h1p(inu)
             
             SigC(p) = SigC(p) + U1_2h1p(inu) *  K_2h1p(inu) * U2_2h1p(inu)
             Z(p)    = Z(p)    + U1_2h1p(inu) * dK_2h1p(inu) * U2_2h1p(inu)
             
             SigC(p) = SigC(p) + U3_2h1p(inu) *  K_2h1p(inu) * U1_2h1p(inu)
             Z(p)    = Z(p)    + U3_2h1p(inu) * dK_2h1p(inu) * U1_2h1p(inu)
             
             SigC(p) = SigC(p) + U1_2h1p(inu) *  K_2h1p(inu) * U3_2h1p(inu)
             Z(p)    = Z(p)    + U1_2h1p(inu) * dK_2h1p(inu) * U3_2h1p(inu)

             jmu = 0
             do j=nC+1,nO
                do mu=1,nS
                   jmu = jmu + 1
             
                   SigC(p) = SigC(p) + U1_2h1p(inu) *  K_2h1p(inu) * C1_2h1p(inu,jmu) *  K_2h1p(jmu) * U1_2h1p(jmu)
                   Z(p)    = Z(p)    + U1_2h1p(inu) * dK_2h1p(inu) * C1_2h1p(inu,jmu) *  K_2h1p(jmu) * U1_2h1p(jmu)
                   Z(p)    = Z(p)    + U1_2h1p(inu) *  K_2h1p(inu) * C1_2h1p(inu,jmu) * dK_2h1p(jmu) * U1_2h1p(jmu)

                end do
             end do
              
       end do
    end do

    anu = 0
    do a=nO+1,nOrb-nR
       do nu=1,nS
          anu = anu + 1
           
          SigC(p) = SigC(p) + U1_2p1h(anu) *  K_2p1h(anu) * U1_2p1h(anu)
          Z(p)    = Z(p)    + U1_2p1h(anu) * dK_2p1h(anu) * U1_2p1h(anu)

          SigC(p) = SigC(p) + U2_2p1h(anu) *  K_2p1h(anu) * U1_2p1h(anu)
          Z(p)    = Z(p)    + U2_2p1h(anu) * dK_2p1h(anu) * U1_2p1h(anu)
          
          SigC(p) = SigC(p) + U1_2p1h(anu) *  K_2p1h(anu) * U2_2p1h(anu)
          Z(p)    = Z(p)    + U1_2p1h(anu) * dK_2p1h(anu) * U2_2p1h(anu)
          
          SigC(p) = SigC(p) + U3_2p1h(anu) *  K_2p1h(anu) * U1_2p1h(anu)
          Z(p)    = Z(p)    + U3_2p1h(anu) * dK_2p1h(anu) * U1_2p1h(anu)
          
          SigC(p) = SigC(p) + U1_2p1h(anu) *  K_2p1h(anu) * U3_2p1h(anu)
          Z(p)    = Z(p)    + U1_2p1h(anu) * dK_2p1h(anu) * U3_2p1h(anu)

          bmu = 0
          do b=nO+1,nOrb-nR
             do mu=1,nS
                bmu = bmu + 1
                
                SigC(p) = SigC(p) + U1_2p1h(anu) *  K_2p1h(anu) * C1_2p1h(anu,bmu) *  K_2p1h(bmu) * U1_2p1h(bmu)
                Z(p)    = Z(p)    + U1_2p1h(anu) * dK_2p1h(anu) * C1_2p1h(anu,bmu) *  K_2p1h(bmu) * U1_2p1h(bmu)
                Z(p)    = Z(p)    + U1_2p1h(anu) *  K_2p1h(anu) * C1_2p1h(anu,bmu) * dK_2p1h(bmu) * U1_2p1h(bmu)
                
             end do
          end do
             
       end do
    end do

 end do
 

 print*,'Alternative form of the self-energy'
 
 Z(:) = 1d0/(1d0 - Z(:))

end subroutine 
