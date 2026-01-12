subroutine R_linDM_2SOSEX(nOrb,nC,nO,nV,nR,nS,e,Om,rho,ERI,eta,linDM)
  
! Compute the linearized GW density matrix

implicit none
include 'parameters.h'


! Input variables

  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: e(nOrb)
  double precision,intent(in)   :: Om(nS)
  double precision,intent(in)   :: rho(nOrb,nOrb,nS)
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)
  integer,intent(in)            :: eta

! Local variables

  integer                       :: p,q,r
  integer                       :: n
  integer                       :: i,j,k,l
  integer                       :: a,b,c,d
  integer                       :: n2h1p,n2p1h
  integer                       :: mu,nu
  integer                       :: klc,kcd,ija,ijb,iab,jab
  integer                       :: anu,inu,bnu,jnu
  double precision              :: num
  double precision              :: dem,dem1,dem2
  double precision              :: reg,reg1,reg2
  double precision              :: s
  double precision,allocatable  :: U1_2h1p(:,:)
  double precision,allocatable  :: U1_2p1h(:,:)
  double precision,allocatable  :: U2_2h1p(:,:)
  double precision,allocatable  :: U2_2p1h(:,:)
  
! Output variables

  double precision,intent(inout)  :: linDM(nOrb,nOrb)

! Dimension of the 2h1p and 2p1h subspaces

  n2h1p = nO*nO*nV
  n2p1h = nV*nV*nO

! Initialization
  
  linDM(:,:) = 0d0
  s = 500d0
  
! Memory allocation  
  allocate(U1_2h1p(n2h1p,nOrb))
  allocate(U1_2p1h(n2p1h,nOrb))
  allocate(U2_2h1p(n2h1p,nOrb))
  allocate(U2_2p1h(n2p1h,nOrb))

! Computing and storing intermediates (i.e. ADC block) for computation of linDM
  
  !---------------!
  ! Blocks U_2h1p !
  !---------------!
  
  U2_2h1p(:,:) = 0d0
  U1_2h1p(:,:) = 0d0
  
  do p=nC+1,nOrb-nR
     
    inu = 0
    do i=nC+1,nO
      do nu=1,nS
        inu = inu + 1

        ! First-order terms

        U1_2h1p(inu,p) = sqrt(2d0)*rho(p,i,nu)

        ! Second-order terms

        do k=nC+1,nO
          do c=nO+1,nOrb-nR

          num = sqrt(2d0)*rho(k,c,nu)*ERI(i,k,c,p)
          dem = e(c) - e(k) - Om(nu)
          reg = (1d0 - exp(-2d0*s*dem*dem))/dem

          U2_2h1p(inu,p) = U2_2h1p(inu,p) + num*reg

          num = sqrt(2d0)*rho(c,k,nu)*ERI(i,c,k,p)
          dem = e(c) - e(k) + Om(nu)
          reg = (1d0 - exp(-2d0*s*dem*dem))/dem

          U2_2h1p(inu,p) = U2_2h1p(inu,p) + num*reg

          end do
        end do

       end do ! nu
     end do ! i
    
   end do ! p


   !---------------!
   ! Blocks U_2p1h !
   !---------------!

    U2_2p1h(:,:) = 0d0
    U1_2p1h(:,:) = 0d0
  
    do p=nC+1,nOrb-nR
       
      anu = 0
      do a=nO+1,nOrb-nR
        do nu=1,nS
          anu = anu + 1

          ! First-order terms

          U1_2p1h(anu,p) = sqrt(2d0)*rho(a,p,nu)
          
          ! Second-order terms

          do k=nC+1,nO
             do c=nO+1,nOrb-nR
                
                num = sqrt(2d0)*rho(k,c,nu)*ERI(a,c,k,p)
                dem = e(c) - e(k) - Om(nu)
                reg = (1d0 - exp(-2d0*s*dem*dem))/dem
                
                U2_2p1h(anu,p) = U2_2p1h(anu,p) + num*reg
                
                num = sqrt(2d0)*rho(c,k,nu)*ERI(a,k,c,p)
                dem = e(c) - e(k) + Om(nu)
                reg = (1d0 - exp(-2d0*s*dem*dem))/dem
                
                U2_2p1h(anu,p) = U2_2p1h(anu,p) + num*reg
                
             end do
          end do

        end do ! nu
      end do ! a
     
    end do ! p

  
! OccOcc block of the density matrix
  
  do i=nC+1,nO
     do j=nC+1,nO

        ! Single pole terms
        anu = 0
        do a=nO+1,nOrb-nR
           do nu=1,nS
              anu = anu + 1
              
              num = - 2d0*U1_2p1h(anu,i)*U1_2p1h(anu,j)
              num = num - 2d0*U1_2p1h(anu,i)*U2_2p1h(anu,j)
              num = num - 2d0*U2_2p1h(anu,i)*U1_2p1h(anu,j)
              dem1 = e(i) - e(a) - Om(nu)
              dem2 = e(j) - e(a) - Om(nu)
 
              reg1 = (1d0 - exp(-2d0*s*dem1*dem1))/dem1
              reg2 = (1d0 - exp(-2d0*s*dem2*dem2))/dem2

              linDM(i,j) = linDM(i,j) + num*reg1*reg2
              
           end do ! nu
        end do ! a
        
     end do ! j 
  end do ! i

  ! VirVir block of the density matrix
  
  do a=nO+1,nOrb-nR
     do b=nO+1,nOrb-nR
        
        ! Single pole terms
        inu=0
        do i=nC+1,nO
           do nu=1,nS
              inu = inu + 1
              
              num = 2d0*U1_2h1p(inu,a)*U1_2h1p(inu,b)
              num = num + 2d0*U1_2h1p(inu,a)*U2_2h1p(inu,b)
              num = num + 2d0*U2_2h1p(inu,a)*U1_2h1p(inu,b)
              dem1 = e(i) - e(a) - Om(nu)
              dem2 = e(i) - e(b) - Om(nu)

              reg1 = (1d0 - exp(-2d0*s*dem1*dem1))/dem1
              reg2 = (1d0 - exp(-2d0*s*dem2*dem2))/dem2

              linDM(a,b) = linDM(a,b) + num*reg1*reg2
              
           end do ! nu
        end do ! i
        
     end do ! b
  end do ! a

  ! OccVir block of the density matrix
  
  do i=nC+1,nO
     do a=nO+1,nOrb-nR

        ! Single pole terms
        jnu = 0
        do j=nC+1,nO
           do nu=1,nS
              jnu = jnu + 1
              
              num = - 2d0 * U1_2h1p(jnu,i) * U1_2h1p(jnu,a)
              num = num - 2d0 * U1_2h1p(jnu,i) * U2_2h1p(jnu,a)
              num = num - 2d0 * U2_2h1p(jnu,i) * U1_2h1p(jnu,a)
              dem1 = e(i) - e(a)
              dem2 = e(j) - e(a) - Om(nu)

              reg1 = (1d0 - exp(-2d0*s*dem1*dem1))/dem1
              reg2 = (1d0 - exp(-2d0*s*dem2*dem2))/dem2

              linDM(i,a) = linDM(i,a) + num*reg1*reg2
              linDM(a,i) = linDM(a,i) + num*reg1*reg2
              
           end do
        end do

        ! Single pole terms
        bnu=0
        do b=nO+1,nOrb-nR
           do nu=1,nS
              bnu = bnu + 1
              
              num = 2d0 * U1_2p1h(bnu,i) * U1_2p1h(bnu,a)
              num = num + 2d0 * U1_2p1h(bnu,i) * U2_2p1h(bnu,a)
              num = num + 2d0 * U2_2p1h(bnu,i) * U1_2p1h(bnu,a)
              dem1 = e(i) - e(a)
              dem2 = e(i) - e(b) - Om(nu)

              reg1 = (1d0 - exp(-2d0*s*dem1*dem1))/dem1
              reg2 = (1d0 - exp(-2d0*s*dem2*dem2))/dem2

              linDM(i,a) = linDM(i,a) + num*reg1*reg2
              linDM(a,i) = linDM(a,i) + num*reg1*reg2
              
           end do
        end do
        
     end do ! a
  end do ! i

  deallocate(U1_2h1p,U2_2h1p,U1_2p1h,U2_2p1h)
  
end subroutine R_linDM_2SOSEX
