subroutine R_linDM_ADC3(nOrb,nC,nO,nV,nR,nS,e,Om,rho,ERI,eta,linDM)
  
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
  integer                       :: amu,imu,bmu,jmu,cmu,kmu
  double precision              :: num
  double precision              :: dem,dem1,dem2,dem3
  double precision              :: reg,reg1,reg2,reg3
  double precision              :: s
  double precision,allocatable  :: U1_2h1p(:,:)
  double precision,allocatable  :: U1_2p1h(:,:)
  double precision,allocatable  :: U2_2h1p(:,:)
  double precision,allocatable  :: U2_2p1h(:,:)
  double precision,allocatable  :: C1_2h1p_2h1p(:,:)
  double precision,allocatable  :: C1_2p1h_2p1h(:,:)
  
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

  allocate(C1_2h1p_2h1p(n2h1p,n2h1p))
  allocate(C1_2p1h_2p1h(n2p1h,n2p1h))

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

!--------------------!
! Block C1_2h1p-2h1p !
!--------------------!
  
  C1_2h1p_2h1p(:,:) = 0d0

  inu = 0
  do i=nC+1,nO
    do nu=1,nS
      inu = inu + 1
  
      ! First-order terms
 
      jmu = 0
      do j=nC+1,nO
        do mu=1,nS
          jmu = jmu + 1
     
          do r=nC+1,nOrb-nR

            num = rho(j,r,nu)*rho(i,r,mu)
            dem = e(i) - e(r) + Om(mu)
            reg = (1d0 - exp(-2d0*s*dem*dem))/dem
           
            C1_2h1p_2h1p(inu,jmu) = C1_2h1p_2h1p(inu,jmu) + num*reg
           
            num = rho(j,r,nu)*rho(i,r,mu)
            dem = e(j) - e(r) + Om(nu)
            reg = (1d0 - exp(-2d0*s*dem*dem))/dem
           
            C1_2h1p_2h1p(inu,jmu) = C1_2h1p_2h1p(inu,jmu) + num*reg

          end do ! r
  
        end do ! mu
      end do ! j
  
    end do ! nu
  end do ! i
 
!--------------------!
! Block C1_2p1h-2p1h !
!--------------------!

  C1_2p1h_2p1h(:,:) = 0d0

  anu = 0
  do a=nO+1,nOrb-nR
    do nu=1,nS
      anu = anu + 1
 
      ! First-order terms
 
      bmu = 0
      do b=nO+1,nOrb-nR
        do mu=1,nS
          bmu = bmu + 1
     
          do r=nC+1,nOrb-nR

            num = rho(r,b,nu)*rho(r,a,mu)
            dem = e(b) - e(r) - Om(nu)
            reg = (1d0 - exp(-2d0*s*dem*dem))/dem
           
            C1_2p1h_2p1h(anu,bmu) = C1_2p1h_2p1h(anu,bmu) + num*reg
           
            num = rho(r,b,nu)*rho(r,a,mu)
            dem = e(a) - e(r) - Om(mu)
            reg = (1d0 - exp(-2d0*s*dem*dem))/dem
           
            C1_2p1h_2p1h(anu,bmu) = C1_2p1h_2p1h(anu,bmu) + num*reg

          end do ! r
  
        end do ! mu
      end do ! b
  
    end do ! nu
  end do ! a
    
! Computing the linDM

  !----------------!
  !  OccOcc block  !
  !----------------!
  
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

        ! Double pole terms
        anu = 0
        do a=nO+1,nOrb-nR
          do nu=1,nS
            anu = anu + 1
            
            bmu = 0
            do b=nO+1,nOrb-nR
              do mu=1,nS
                bmu = bmu + 1
            
                num = 2d0*U1_2p1h(anu,i)*C1_2p1h_2p1h(anu,bmu)*U1_2p1h(bmu,j)
                dem1 = e(i) - e(a) - Om(nu)
                dem3 = e(i) - e(b) - Om(mu)
                dem2 = e(j) - e(a) - Om(nu)
                
                reg1 = (1d0 - exp(-2d0*s*dem1*dem1))/dem1
                reg2 = (1d0 - exp(-2d0*s*dem2*dem2))/dem2
                reg3 = (1d0 - exp(-2d0*s*dem3*dem3))/dem3
                
                linDM(i,j) = linDM(i,j) + num*reg1*reg2*reg3

                num = 2d0*U1_2p1h(anu,i)*C1_2p1h_2p1h(anu,bmu)*U1_2p1h(bmu,j)
                dem1 = e(i) - e(b) - Om(mu)
                dem3 = e(j) - e(a) - Om(nu)
                dem2 = e(j) - e(b) - Om(mu)
                
                reg1 = (1d0 - exp(-2d0*s*dem1*dem1))/dem1
                reg2 = (1d0 - exp(-2d0*s*dem2*dem2))/dem2
                reg3 = (1d0 - exp(-2d0*s*dem3*dem3))/dem3
                
                linDM(i,j) = linDM(i,j) + num*reg1*reg2*reg3

              end do ! mu
           end do ! b           
            
          end do ! nu
        end do ! a

        
     end do ! j 
  end do ! i

  !----------------!
  !  VirVir block  !
  !----------------!
  
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

        ! Double pole terms
        inu=0
        do i=nC+1,nO
           do nu=1,nS
              inu = inu + 1

              jmu=0
              do j=nC+1,nO
                 do mu=1,nS
                    jmu = jmu + 1
              
                    num = 2d0*U1_2h1p(inu,a)*C1_2h1p_2h1p(inu,jmu)*U1_2h1p(jmu,b)
                    dem1 = e(i) - e(a) - Om(nu)
                    dem2 = e(j) - e(a) - Om(mu)
                    dem3 = e(i) - e(b) - Om(nu)
                    
                    reg1 = (1d0 - exp(-2d0*s*dem1*dem1))/dem1
                    reg2 = (1d0 - exp(-2d0*s*dem2*dem2))/dem2
                    reg3 = (1d0 - exp(-2d0*s*dem3*dem3))/dem3

                    linDM(a,b) = linDM(a,b) + num*reg1*reg2*reg3

                    num = 2d0*U1_2h1p(inu,a)*C1_2h1p_2h1p(inu,jmu)*U1_2h1p(jmu,b)
                    dem1 = e(j) - e(a) - Om(mu)
                    dem2 = e(i) - e(b) - Om(nu)
                    dem3 = e(j) - e(b) - Om(mu)
                    
                    reg1 = (1d0 - exp(-2d0*s*dem1*dem1))/dem1
                    reg2 = (1d0 - exp(-2d0*s*dem2*dem2))/dem2
                    reg3 = (1d0 - exp(-2d0*s*dem3*dem3))/dem3

                    linDM(a,b) = linDM(a,b) + num*reg1*reg2*reg3
                    
                 end do ! mu
              end do ! j
              
           end do ! nu
        end do ! i
        
     end do ! b
  end do ! a

  !----------------!
  !  OccVir block  !
  !----------------!
  
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

        ! Double pole terms
        bnu = 0
        do b=nO+1,nOrb-nR
          do nu=1,nS
            bnu = bnu + 1
            
            cmu = 0
            do c=nO+1,nOrb-nR
              do mu=1,nS
                cmu = cmu + 1
            
                num = 2d0*U1_2p1h(bnu,i)*C1_2p1h_2p1h(bnu,cmu)*U1_2p1h(cmu,a)
                dem1 = e(i) - e(a)
                dem3 = e(i) - e(b) - Om(nu)
                dem2 = e(i) - e(c) - Om(mu)
                
                reg1 = (1d0 - exp(-2d0*s*dem1*dem1))/dem1
                reg2 = (1d0 - exp(-2d0*s*dem2*dem2))/dem2
                reg3 = (1d0 - exp(-2d0*s*dem3*dem3))/dem3
                
                linDM(i,a) = linDM(i,a) + num*reg1*reg2*reg3

                num = 2d0*U1_2p1h(bnu,a)*C1_2p1h_2p1h(bnu,cmu)*U1_2p1h(cmu,i)
                dem1 = e(i) - e(a)
                dem3 = e(i) - e(b) - Om(nu)
                dem2 = e(i) - e(c) - Om(mu)
                
                reg1 = (1d0 - exp(-2d0*s*dem1*dem1))/dem1
                reg2 = (1d0 - exp(-2d0*s*dem2*dem2))/dem2
                reg3 = (1d0 - exp(-2d0*s*dem3*dem3))/dem3
                
                linDM(a,i) = linDM(a,i) + num*reg1*reg2*reg3

              end do ! mu
           end do ! c           
            
          end do ! nu
        end do ! b
        
        jnu=0
        do j=nC+1,nO
           do nu=1,nS
              jnu = jnu + 1

              kmu=0
              do k=nC+1,nO
                 do mu=1,nS
                    kmu = kmu + 1
              
                    num = - 2d0*U1_2h1p(jnu,i)*C1_2h1p_2h1p(jnu,kmu)*U1_2h1p(kmu,a)
                    dem1 = e(a) - e(i)
                    dem2 = e(a) - e(j) + Om(nu)
                    dem3 = e(a) - e(k) + Om(mu)
                    
                    reg1 = (1d0 - exp(-2d0*s*dem1*dem1))/dem1
                    reg2 = (1d0 - exp(-2d0*s*dem2*dem2))/dem2
                    reg3 = (1d0 - exp(-2d0*s*dem3*dem3))/dem3

                    linDM(i,a) = linDM(i,a) + num*reg1*reg2*reg3

                    num = - 2d0*U1_2h1p(jnu,a)*C1_2h1p_2h1p(jnu,kmu)*U1_2h1p(kmu,i)
                    dem1 = e(a) - e(i)
                    dem2 = e(a) - e(j) + Om(nu)
                    dem3 = e(a) - e(k) + Om(mu)
                    
                    reg1 = (1d0 - exp(-2d0*s*dem1*dem1))/dem1
                    reg2 = (1d0 - exp(-2d0*s*dem2*dem2))/dem2
                    reg3 = (1d0 - exp(-2d0*s*dem3*dem3))/dem3

                    linDM(a,i) = linDM(a,i) + num*reg1*reg2*reg3
                    
                 end do ! mu
              end do ! k
              
           end do ! nu
        end do ! j
        
     end do ! a
  end do ! i

  deallocate(U1_2h1p,U2_2h1p,U1_2p1h,U2_2p1h)
  
end subroutine R_linDM_ADC3
