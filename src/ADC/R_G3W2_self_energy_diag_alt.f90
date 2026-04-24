subroutine R_G3W2_self_energy_diag_alt(eta,nBas,nOrb,nC,nO,nV,nR,nS,eHF,Om,rho,ERI,EcGM,SigC,Z)

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
  integer                       :: inu,jmu,anu,bmu
  integer                       :: n2h1p,n2p1h
  integer                       :: mu,nu
  integer                       :: klc,kcd,ija,ijb,iab,jab
  double precision              :: w
  double precision              :: num,num1,num2
  double precision              :: dem,dem1,dem2,dem3
  double precision              :: reg,reg1,reg2,reg3

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


  double precision              :: flow = 1d5

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
            reg = (1d0 - exp(-2d0*flow*dem*dem))/dem
         
            C1_2h1p(inu,jmu) = C1_2h1p(inu,jmu) + num*reg
         
            num = rho(i,k,mu)*rho(j,k,nu)
            dem = eHF(j) - eHF(k) + Om(nu)
            reg = (1d0 - exp(-2d0*flow*dem*dem))/dem
         
            C1_2h1p(inu,jmu) = C1_2h1p(inu,jmu) + num*reg

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
            reg = (1d0 - exp(-2d0*flow*dem*dem))/dem
         
            C1_2p1h(anu,bmu) = C1_2p1h(anu,bmu) + num*reg
         
            num = rho(c,a,mu)*rho(c,b,nu)
            dem = eHF(b) - eHF(c) - Om(nu)
            reg = (1d0 - exp(-2d0*flow*dem*dem))/dem
         
            C1_2p1h(anu,bmu) = C1_2p1h(anu,bmu) + num*reg

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
             do c=nO+1,nOrb-nR
                do d=nO+1,nOrb-nR
                   
                   num1 = 2d0*rho(c,i,mu)*rho(p,c,nu)
                   num2 = 2d0*rho(d,i,mu)*rho(p,d,nu)
                   dem1 = eHF(c) - eHF(i) + Om(mu)
                   dem2 = w - eHF(i) + Om(nu) + Om(mu)
                   dem3 = eHF(d) - eHF(i) + Om(mu)
                   
                   reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))/dem1
                   reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))/dem2
                   reg3 = (1d0 - exp(-2d0*flow*dem3*dem3))/dem3
                   
                   SigC(p) = SigC(p) + num1*num2*reg1*reg2*reg3
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
             do i=nC+1,nO
                do j=nC+1,nO

                   num1 = 2d0*rho(a,i,mu)*rho(i,p,nu)
                   num2 = 2d0*rho(a,j,mu)*rho(j,p,nu)
                   dem1 = eHF(i) - eHF(a) - Om(mu)
                   dem2 = w - eHF(a) - Om(nu) - Om(mu)
                   dem3 = eHF(j) - eHF(a) - Om(mu)

                   reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))/dem1
                   reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))/dem2
                   reg3 = (1d0 - exp(-2d0*flow*dem3*dem3))/dem3
                   
                   SigC(p) = SigC(p) + num1*num2*reg1*reg2*reg3
                   Z(p)   = Z(p)   - num1*num2*reg1*reg2*reg3/dem2

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
                 reg = (1d0 - exp(-2d0*flow*dem*dem))/dem
                 
                 U2_2h1p(inu) = U2_2h1p(inu) + num*reg
          
                 num = sqrt(2d0)*rho(c,k,nu)*ERI(i,c,k,p)
                 dem = eHF(c) - eHF(k) + Om(nu)
                 reg = (1d0 - exp(-2d0*flow*dem*dem))/dem

                 U2_2h1p(inu) = U2_2h1p(inu) + num*reg
                 
              end do
           end do

           ! Third-order terms

           do k=nC+1,nO
              do c=nO+1,nOrb-nR
                 do mu=1,nS

                    num = 2d0*sqrt(2d0)*rho(k,c,nu)*rho(c,i,mu)*rho(k,p,mu)
                    dem1 = eHF(c) - eHF(k) - Om(nu)
                    dem2 = eHF(c) - eHF(i) + Om(mu)

                    reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))/dem1
                    reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))/dem2

                    U3_2h1p(inu) = U3_2h1p(inu) - num*reg1*reg2

                    num = 2d0*sqrt(2d0)*rho(k,c,nu)*rho(i,c,mu)*rho(p,k,mu)
                    dem1 = eHF(c) - eHF(k) - Om(nu)
                    dem2 = eHF(c) - eHF(i) - Om(mu)

                    reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))/dem1
                    reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))/dem2

                    U3_2h1p(inu) = U3_2h1p(inu) + 0.5d0*num*reg1*reg2

                    num = 2d0*sqrt(2d0)*rho(k,i,mu)*rho(c,k,nu)*rho(c,p,mu)
                    dem1 = eHF(c) - eHF(i) + Om(mu) + Om(nu)
                    dem2 = eHF(c) - eHF(k) + Om(nu)

                    reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))/dem1
                    reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))/dem2
                  
                    U3_2h1p(inu) = U3_2h1p(inu) - num*reg1*reg2

                 end do
              end do
           end do

           do a=nO+1,nOrb-nR
            do b=nO+1,nOrb-nR
              do mu=1,nS
    
                num = 2d0*sqrt(2d0)*rho(a,i,mu)*rho(b,a,nu)*rho(b,p,mu)
                dem1 = eHF(b) - eHF(i) + Om(mu) + Om(nu)
                dem2 = eHF(a) - eHF(i) + Om(mu)

                reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))/dem1
                reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))/dem2
                  
                U3_2h1p(inu) = U3_2h1p(inu) + num*reg1*reg2

 
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
                 reg = (1d0 - exp(-2d0*flow*dem*dem))/dem

                 U2_2p1h(anu) = U2_2p1h(anu) + num*reg

                 num = sqrt(2d0)*rho(c,k,nu)*ERI(a,k,c,p)
                 dem = eHF(c) - eHF(k) + Om(nu)
                 reg = (1d0 - exp(-2d0*flow*dem*dem))/dem

                 U2_2p1h(anu) = U2_2p1h(anu) + num*reg

              end do
           end do

           ! Third-order terms
           
           do k=nC+1,nO
              do c=nO+1,nOrb-nR
                 do mu=1,nS

                    num = 2d0*sqrt(2d0)*rho(k,c,nu)*rho(a,k,mu)*rho(p,c,mu)
                    dem1 = eHF(c) - eHF(k) - Om(nu)
                    dem2 = eHF(a) - eHF(k) + Om(mu)

                    reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))/dem1
                    reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))/dem2

                    U3_2p1h(anu) = U3_2p1h(anu) - num*reg1*reg2
                  
                    num = 2d0*sqrt(2d0)*rho(k,c,nu)*rho(k,a,mu)*rho(c,p,mu)
                    dem1 = eHF(c) - eHF(k) - Om(nu)
                    dem2 = eHF(a) - eHF(k) - Om(mu)

                    reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))/dem1
                    reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))/dem2

                    U3_2p1h(anu) = U3_2p1h(anu) + 0.5d0*num*reg1*reg2

                    num = 2d0*sqrt(2d0)*rho(a,c,mu)*rho(c,k,nu)*rho(p,k,mu)
                    dem1 = eHF(a) - eHF(k) + Om(mu) + Om(nu)
                    dem2 = eHF(c) - eHF(k) + Om(nu)

                    reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))/dem1
                    reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))/dem2

                    U3_2p1h(anu) = U3_2p1h(anu) - num*reg1*reg2

                 end do
              end do
           end do

           do i=nC+1,nO
              do j=nC+1,nO
                 do mu=1,nS
                    
                    num = 2d0*sqrt(2d0)*rho(a,j,mu)*rho(j,i,nu)*rho(p,i,mu)
                    dem1 = eHF(a) - eHF(i) + Om(mu) + Om(nu)
                    dem2 = eHF(a) - eHF(j) + Om(mu)

                    reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))/dem1
                    reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))/dem2
            
                    U3_2p1h(anu) = U3_2p1h(anu) + num*reg1*reg2

                 end do
              end do
           end do

        end do
    end do

    !--------------!
    ! Block K_2h1p !
    !--------------!

    K_2h1p(:) = 0d0
    dK_2h1p(:) = 0d0

    ija = 0
    do i=nC+1,nO
       do mu=1,nS
          ija = ija + 1
          
          ! Zeroth-order terms
          K_2h1p(ija)  =   1d0/(w - eHF(i) + Om(mu))
          dK_2h1p(ija) = - 1d0/(w - eHF(i) + Om(mu))**2

      end do
    end do

    !--------------!
    ! Block K_2p1h !
    !--------------!

    K_2p1h(:)  = 0d0
    dK_2p1h(:) = 0d0

    iab = 0
    do a=nO+1,nOrb-nR
       do mu=1,nS
          iab = iab + 1
          
          ! Zeroth-order terms
          K_2p1h(iab)  =   1d0/(w - eHF(a) - Om(mu))
          dK_2p1h(iab) = - 1d0/(w - eHF(a) - Om(mu))**2

       end do
    end do

!----------------------------!
!    Building self-energy    !
!----------------------------!
    ija = 0
    do i=nC+1,nO
       do j=nC+1,nO
          do a=nO+1,nBas-nR
             ija = ija + 1
             
             SigC(p) = SigC(p) + U1_2h1p(ija) *  K_2h1p(ija) * U1_2h1p(ija)
             Z(p)    = Z(p)    + U1_2h1p(ija) * dK_2h1p(ija) * U1_2h1p(ija)
             
             SigC(p) = SigC(p) + U2_2h1p(ija) *  K_2h1p(ija) * U1_2h1p(ija)
             Z(p)    = Z(p)    + U2_2h1p(ija) * dK_2h1p(ija) * U1_2h1p(ija)
             
             SigC(p) = SigC(p) + U1_2h1p(ija) *  K_2h1p(ija) * U2_2h1p(ija)
             Z(p)    = Z(p)    + U1_2h1p(ija) * dK_2h1p(ija) * U2_2h1p(ija)
             
             SigC(p) = SigC(p) + U3_2h1p(ija) *  K_2h1p(ija) * U1_2h1p(ija)
             Z(p)    = Z(p)    + U3_2h1p(ija) * dK_2h1p(ija) * U1_2h1p(ija)
             
             SigC(p) = SigC(p) + U1_2h1p(ija) *  K_2h1p(ija) * U3_2h1p(ija)
             Z(p)    = Z(p)    + U1_2h1p(ija) * dK_2h1p(ija) * U3_2h1p(ija)

             klc = 0
             do k=nC+1,nO
                do l=nC+1,nO
                   do c=nO+1,nBas-nR
                      klc = klc + 1
             
                      SigC(p) = SigC(p) + U1_2h1p(ija) *  K_2h1p(ija) * C1_2h1p(ija,klc) *  K_2h1p(klc) * U1_2h1p(klc)
                      Z(p)    = Z(p)    + U1_2h1p(ija) * dK_2h1p(ija) * C1_2h1p(ija,klc) *  K_2h1p(klc) * U1_2h1p(klc)
                      Z(p)    = Z(p)    + U1_2h1p(ija) *  K_2h1p(ija) * C1_2h1p(ija,klc) * dK_2h1p(klc) * U1_2h1p(klc)

                   end do
                end do
             end do
              
          end do
       end do
    end do
    iab = 0
    do i=nC+1,nO
       do a=nO+1,nBas-nR
          do b=nO+1,nBas-nR
             iab = iab + 1
             
             SigC(p) = SigC(p) + U1_2p1h(iab) *  K_2p1h(iab) * U1_2p1h(iab)
             Z(p)    = Z(p)    + U1_2p1h(iab) * dK_2p1h(iab) * U1_2p1h(iab)

             SigC(p) = SigC(p) + U2_2p1h(iab) *  K_2p1h(iab) * U1_2p1h(iab)
             Z(p)    = Z(p)    + U2_2p1h(iab) * dK_2p1h(iab) * U1_2p1h(iab)

             SigC(p) = SigC(p) + U1_2p1h(iab) *  K_2p1h(iab) * U2_2p1h(iab)
             Z(p)    = Z(p)    + U1_2p1h(iab) * dK_2p1h(iab) * U2_2p1h(iab)

             SigC(p) = SigC(p) + U3_2p1h(iab) *  K_2p1h(iab) * U1_2p1h(iab)
             Z(p)    = Z(p)    + U3_2p1h(iab) * dK_2p1h(iab) * U1_2p1h(iab)

             SigC(p) = SigC(p) + U1_2p1h(iab) *  K_2p1h(iab) * U3_2p1h(iab)
             Z(p)    = Z(p)    + U1_2p1h(iab) * dK_2p1h(iab) * U3_2p1h(iab)

             kcd = 0
             do k=nC+1,nO
                do c=nO+1,nBas-nR
                   do d=nO+1,nBas-nR
                      kcd = kcd + 1
                      
                      SigC(p) = SigC(p) + U1_2p1h(iab) *  K_2p1h(iab) * C1_2p1h(iab,kcd) *  K_2p1h(kcd) * U1_2p1h(kcd)
                      Z(p)    = Z(p)    + U1_2p1h(iab) * dK_2p1h(iab) * C1_2p1h(iab,kcd) *  K_2p1h(kcd) * U1_2p1h(kcd)
                      Z(p)    = Z(p)    + U1_2p1h(iab) *  K_2p1h(iab) * C1_2p1h(iab,kcd) * dK_2p1h(kcd) * U1_2p1h(kcd)
                      
                   end do
                end do
             end do
             
          end do
       end do
    end do

 end do
 

 print*,'Alternative form of the self-energy'
 
 Z(:) = 1d0/(1d0 - Z(:))

end subroutine 
