subroutine R_psdG3W2_self_energy_iieta(p,w,eta,nOrb,nC,nO,nV,nR,nS,eHF,Om,rho,ERI,Re_SigC,Im_SigC,Z)

! Compute diagonal of the correlation part of the self-energy and the renormalization factor
! for the G3W2 approximation

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: p
  double precision,intent(in)   :: w

  double precision,intent(in)   :: eta
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
  
  integer                       :: i,j,k,l
  integer                       :: a,b,c,d
  integer                       :: jb,kc,ia,ja
  integer                       :: inu,jmu,anu,bmu
  integer                       :: imu,knu,amu,cnu
  integer                       :: kla,cla
  integer                       :: n2h1p,n2p1h
  integer                       :: mu,nu,la,s,t
  integer                       :: klc,kcd,ija,ijb,iab,jab
  complex *16                   :: num,num1,num2,num3
  complex *16                   :: dem,dem1,dem2,dem3
  logical                       :: do_psd = .true.

  complex *16,allocatable  :: C1_2h1p(:,:)
  complex *16,allocatable  :: C1_2p1h(:,:)
  complex *16,allocatable  :: K_2h1p(:)
  complex *16,allocatable  :: K_2p1h(:)
  complex *16,allocatable  :: dK_2h1p(:)
  complex *16,allocatable  :: dK_2p1h(:)
  complex *16,allocatable  :: U1_2h1p(:)
  complex *16,allocatable  :: U1_2p1h(:)
  complex *16,allocatable  :: U2_2h1p(:)
  complex *16,allocatable  :: U2_2p1h(:)
  complex *16,allocatable  :: U3_2h1p(:)
  complex *16,allocatable  :: U3_2p1h(:)
  complex *16,allocatable  :: tmp_2h1p(:)
  complex *16,allocatable  :: tmp_2p1h(:)

! Output variables
  double precision,intent(out)  :: Re_SigC,Im_SigC
  double precision,intent(out)  :: Z

! Initialize 

  Re_SigC = 0d0
  Im_SigC = 0d0
  Z    = 0d0

! Dimension of the 2h1p and 2p1h subspaces

  n2h1p = nO*nO*nV
  n2p1h = nV*nV*nO

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
    
  allocate(tmp_2h1p(n2h1p))
  allocate(tmp_2p1h(n2p1h))

!---------------!
! Block C1_2h1p !
!---------------!
  
  C1_2h1p(:,:) = 0d0

  !$OMP PARALLEL DEFAULT(NONE)                                  &
  !$OMP          PRIVATE(i, nu, j, mu, inu, jmu, num, dem) &
  !$OMP          SHARED(nC, nO, nOrb, nR, nS, eta, eHF, rho, Om, C1_2h1p)
  !$OMP DO COLLAPSE(2)
  do i=nC+1,nO
    do nu=1,nS
      inu = (i - nC - 1) * nS + nu
  
      ! First-order terms
      jmu = 0
      do j=nC+1,nO
        do mu=1,nS
          jmu = jmu + 1
   
          do c=nO+1,nOrb-nR

            num = rho(i,c,mu)*rho(j,c,nu)
            dem = eHF(i) - eHF(c) + Om(mu) - im*eta
         
            C1_2h1p(inu,jmu) = C1_2h1p(inu,jmu) + num/dem
         
            num = rho(i,c,mu)*rho(j,c,nu)
            dem = eHF(j) - eHF(c) + Om(nu) - im*eta
         
            C1_2h1p(inu,jmu) = C1_2h1p(inu,jmu) + num/dem

          end do
  
        end do
      end do
  
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
 
!---------------!
! Block C1_2p1h !
!---------------!

  C1_2p1h(:,:) = 0d0

  !$OMP PARALLEL DEFAULT(NONE)                                  &
  !$OMP          PRIVATE(a, nu, b, mu, anu, bmu, num, dem) &
  !$OMP          SHARED(nC, nO, nOrb, nR, nS, eta, eHF, rho, Om, C1_2p1h)
  !$OMP DO COLLAPSE(2)
  do a=nO+1,nOrb-nR
    do nu=1,nS
      anu = (a - nO - 1) * nS + nu
 
      ! First-order terms
 
      bmu = 0
      do b=nO+1,nOrb-nR
        do mu=1,nS
          bmu = bmu + 1
   
          do k=nC+1,nO

            num = rho(k,a,mu)*rho(k,b,nu)
            dem = eHF(a) - eHF(k) - Om(mu) + im*eta
         
            C1_2p1h(anu,bmu) = C1_2p1h(anu,bmu) + num/dem
         
            num = rho(k,a,mu)*rho(k,b,nu)
            dem = eHF(b) - eHF(k) - Om(nu) + im*eta
         
            C1_2p1h(anu,bmu) = C1_2p1h(anu,bmu) + num/dem

          end do
  
        end do
      end do
  
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  
  ! Downfolding the 3h2p configurations
  
  !$OMP PARALLEL DEFAULT(NONE)                                  &
  !$OMP          PRIVATE(i,j,k,a,b, nu, mu, num1,num2,num3, dem1,dem2,dem3) &
  !$OMP          SHARED(p, w, nC, nO, nOrb, nR, nS, eta, eHF, rho, Om) &
  !$OMP REDUCTION(+:Re_SigC,Im_SigC,Z)
  !$OMP DO COLLAPSE(3)
  do i=nC+1,nO
     do nu=1,nS
        do mu=1,nS
           
           do j=nC+1,nO
              do k=nC+1,nO
                 
                 num1 = 2d0*rho(k,i,nu)*rho(p,k,mu)
                 dem1 = w - eHF(k) + Om(mu) - 2d0*im*eta
                 
                 dem2 = w - eHF(i) + Om(nu) + Om(mu) - 3d0*im*eta
                 
                 num3 = 2d0*rho(j,i,mu)*rho(p,j,nu)
                 dem3 = w - eHF(j) + Om(nu) - 2d0*im*eta
                 
                 Re_SigC = Re_SigC + Real(num1*num3/dem1/dem2/dem3)
                 Im_SigC = Im_SigC + Aimag(num1*num3/dem1/dem2/dem3)
                 
                 Z    = Z    - Real(num1*num3/dem1/dem2/dem3/dem1)
                 Z    = Z    - Real(num1*num3/dem1/dem2/dem3/dem2)
                 Z    = Z    - Real(num1*num3/dem1/dem2/dem3/dem3)
                 
              end do
           end do
           
           do k=nC+1,nO
              do a=nO+1,nOrb-nR
                 
                 num1 = 2d0*rho(k,i,nu)*rho(p,k,mu)
                 dem1 = w - eHF(k) + Om(mu) - 2d0*im*eta
                 
                 dem2 = w - eHF(i) + Om(nu) + Om(mu) - 3d0*im*eta
                 
                 num3 = 2d0*rho(a,i,mu)*rho(p,a,nu)
                 dem3 = eHF(i) - eHF(a) - Om(mu) + 3d0*im*eta
                 
                 Re_SigC = Re_SigC + Real(2d0*num1*num3/dem1/dem2/dem3)
                 Im_SigC = Im_SigC + Aimag(2d0*num1*num3/dem1/dem2/dem3)
                 
                 Z    = Z    - Real(2d0*num1*num3/dem1/dem2/dem3/dem1)
                 Z    = Z    - Real(2d0*num1*num3/dem1/dem2/dem3/dem2)
                 
              end do
           end do
             
           do b=nO+1,nOrb-nR
              do a=nO+1,nOrb-nR
                 
                 num1 = 2d0*rho(b,i,nu)*rho(p,b,mu)
                 dem1 = eHF(i) - eHF(b) - Om(nu) + 3d0*im*eta
                 
                 dem2 = w - eHF(i) + Om(nu) + Om(mu) - 3d0*im*eta
                
                 num3 = 2d0*rho(a,i,mu)*rho(p,a,nu)
                 dem3 = eHF(i) - eHF(a) - Om(mu) + 3d0*im*eta
                 
                 Re_SigC = Re_SigC + Real(num1*num3/dem1/dem2/dem3)
                 Im_SigC = Im_SigC + Aimag(num1*num3/dem1/dem2/dem3)
                 
                 Z    = Z    - Real(num1*num3/dem1/dem2/dem3/dem2)
                 
              end do
           end do
           
           
        end do
     end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  
  ! Downfolding the 3p2h configurations

    
  !$OMP PARALLEL DEFAULT(NONE)                                  &
  !$OMP          PRIVATE(a,b,c,i,j, nu, mu, num1,num2,num3, dem1,dem2,dem3) &
  !$OMP          SHARED(p, w, nC, nO, nOrb, nR, nS, eta, eHF, rho, Om) &
  !$OMP REDUCTION(+:Re_SigC,Im_SigC,Z)
  !$OMP DO COLLAPSE(3)
  do a=nO+1,nOrb-nR
     do nu=1,nS          
        do mu=1,nS
           
           do b=nO+1,nOrb-nR
              do c=nO+1,nOrb-nR
                 
                 num1 = 2d0*rho(a,c,nu)*rho(c,p,mu)
                 dem1 = w - eHF(c) - Om(mu) + 2d0*im*eta
                 
                 dem2 = w - eHF(a) - Om(nu) - Om(mu) + 3d0*im*eta
                 
                 num3 = 2d0*rho(a,b,mu)*rho(b,p,nu)
                 dem3 = w - eHF(b) - Om(nu)+ 2d0*im*eta

                 Re_SigC = Re_SigC + Real(num1*num3/dem1/dem2/dem3)
                 Im_SigC = Im_SigC + Aimag(num1*num3/dem1/dem2/dem3)
                 
                 Z    = Z    - Real(num1*num3/dem1/dem2/dem3/dem1)
                 Z    = Z    - Real(num1*num3/dem1/dem2/dem3/dem2)
                 Z    = Z    - Real(num1*num3/dem1/dem2/dem3/dem3)
                   
              end do
           end do
           
           do c=nO+1,nOrb-nR
              do i=nC+1,nO
                 
                 num1 = 2d0*rho(a,c,nu)*rho(c,p,mu)
                 dem1 = w - eHF(c) - Om(mu) + 2d0*im*eta
                 
                 dem2 = w - eHF(a) - Om(nu) - Om(mu) + 3d0*im*eta
                 
                 num3 = 2d0*rho(a,i,mu)*rho(i,p,nu)
                 dem3 = eHF(a) - eHF(i) + Om(mu) - 3d0*im*eta
                 
                 Re_SigC = Re_SigC + Real(num1*num3/dem1/dem2/dem3)
                 Im_SigC = Im_SigC + Aimag(num1*num3/dem1/dem2/dem3)
                 
                 Z    = Z    - Real(2d0*num1*num3/dem1/dem2/dem3/dem1)
                 Z    = Z    - Real(2d0*num1*num3/dem1/dem2/dem3/dem2)
                 
              end do
           end do
           
           do j=nC+1,nO
              do i=nC+1,nO
                 
                 num1 = 2d0*rho(a,j,nu)*rho(j,p,mu)
                 dem1 = eHF(a) - eHF(j) + Om(nu) - 3d0*im*eta
                 
                 dem2 = w - eHF(a) - Om(nu) - Om(mu) + 3d0*im*eta 
                 
                 num3 = 2d0*rho(a,i,mu)*rho(i,p,nu)
                 dem3 = eHF(a) - eHF(i) + Om(mu) - 3d0*im*eta
                 
                 Re_SigC = Re_SigC + Real(num1*num3/dem1/dem2/dem3)
                 Im_SigC = Im_SigC + Aimag(num1*num3/dem1/dem2/dem3)
                 
                 Z       = Z       - Real(num1*num3/dem1/dem2/dem3/dem2)
                 
              end do
           end do
           
           
        end do
     end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
     
  !---------------!
  ! Blocks U_2h1p !
  !---------------!
  
  U1_2h1p(:) = 0d0
  U2_2h1p(:) = 0d0
  U3_2h1p(:) = 0d0

  !$OMP PARALLEL DEFAULT(NONE)                                  &
  !$OMP          PRIVATE(a,b,c,i,j,k, nu, mu, inu, num,num1,num2,num3, dem,dem1,dem2,dem3) &
  !$OMP          SHARED(p, w, nC, nO, nOrb, nR, nS, eta, eHF, rho, Om, ERI, U1_2h1p,U2_2h1p,U3_2h1p)
  !$OMP DO COLLAPSE(2) 
  do i=nC+1,nO
     do nu=1,nS
        inu = (i - nC - 1) * nS + nu
        
        ! First-order terms
           
        U1_2h1p(inu) = sqrt(2d0)*rho(p,i,nu)
        
        ! Second-order terms
        
        do k=nC+1,nO
           do c=nO+1,nOrb-nR
              
              num = sqrt(2d0)*rho(k,c,nu)*ERI(i,k,c,p)
              dem = eHF(c) - eHF(k) - Om(nu) - im*eta
              
              U2_2h1p(inu) = U2_2h1p(inu) + num/dem
              
              num = sqrt(2d0)*rho(c,k,nu)*ERI(i,c,k,p)
              dem = eHF(c) - eHF(k) + Om(nu) - 3d0*im*eta
              
              U2_2h1p(inu) = U2_2h1p(inu) + num/dem
              
           end do
        end do
        
        ! Third-order terms
        
        do k=nC+1,nO
           do c=nO+1,nOrb-nR
              do mu=1,nS
                 
                 num = 2d0*sqrt(2d0)*rho(k,c,nu)*rho(c,i,mu)*rho(k,p,mu)
                 dem1 = eHF(c) - eHF(k) - Om(nu) - im*eta
                 dem2 = eHF(c) - eHF(i) + Om(mu) - 3d0*im*eta
                 
                 U3_2h1p(inu) = U3_2h1p(inu) - num/dem1/dem2
                 
                 num = 2d0*sqrt(2d0)*rho(k,c,nu)*rho(i,c,mu)*rho(p,k,mu)
                 dem1 = eHF(c) - eHF(k) - Om(nu) - im*eta
                 dem2 = eHF(c) - eHF(i) - Om(mu) - im*eta
                 
                 U3_2h1p(inu) = U3_2h1p(inu) + 0.5d0*num/dem1/dem2

                 num = 2d0*sqrt(2d0)*rho(k,i,mu)*rho(c,k,nu)*rho(c,p,mu)
                 dem1 = eHF(c) - eHF(i) + Om(mu) + Om(nu) - 4d0*im*eta
                 dem2 = eHF(c) - eHF(k) + Om(nu) - 3d0*im*eta
                 
                 U3_2h1p(inu) = U3_2h1p(inu) - num/dem1/dem2
                 
              end do
           end do
        end do
        
        do a=nO+1,nOrb-nR
           do b=nO+1,nOrb-nR
              do mu=1,nS
                 
                 num = 2d0*sqrt(2d0)*rho(a,i,mu)*rho(b,a,nu)*rho(b,p,mu)
                 dem1 = eHF(b) - eHF(i) + Om(mu) + Om(nu) - 4d0*im*eta
                 dem2 = eHF(a) - eHF(i) + Om(mu) - 3d0*im*eta
                 
                 U3_2h1p(inu) = U3_2h1p(inu) + num/dem1/dem2
                 
              end do
           end do
        end do

     end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  !---------------!
  ! Blocks U_2p1h !
  !---------------!
  
  U1_2p1h(:) = 0d0
  U2_2p1h(:) = 0d0
  U3_2p1h(:) = 0d0

  !$OMP PARALLEL DEFAULT(NONE)                                  &
  !$OMP          PRIVATE(a,b,c,i,j,k, nu, mu, anu, num,num1,num2,num3, dem,dem1,dem2,dem3) &
  !$OMP          SHARED(p, w, nC, nO, nOrb, nR, nS, eta, eHF, rho, Om, ERI, U1_2p1h,U2_2p1h,U3_2p1h)
  !$OMP DO COLLAPSE(2)
  do a=nO+1,nOrb-nR
     do nu=1,nS
        anu = (a - nO - 1) * nS + nu

        ! First-order terms
        
        U1_2p1h(anu) = sqrt(2d0)*rho(a,p,nu)
        
        ! Second-order terms
        
        do k=nC+1,nO
           do c=nO+1,nOrb-nR
              
              num = sqrt(2d0)*rho(k,c,nu)*ERI(a,c,k,p)
              dem = eHF(c) - eHF(k) - Om(nu) - im*eta
              
              U2_2p1h(anu) = U2_2p1h(anu) + num/dem
              
              num = sqrt(2d0)*rho(c,k,nu)*ERI(a,k,c,p)
              dem = eHF(c) - eHF(k) + Om(nu) - 3d0*im*eta
              
              U2_2p1h(anu) = U2_2p1h(anu) + num/dem
              
           end do
        end do
        
        ! Third-order terms
        
        do k=nC+1,nO
           do c=nO+1,nOrb-nR
              do mu=1,nS
                 
                 num = 2d0*sqrt(2d0)*rho(k,c,nu)*rho(a,k,mu)*rho(p,c,mu)
                 dem1 = eHF(c) - eHF(k) - Om(nu) - im*eta
                 dem2 = eHF(a) - eHF(k) + Om(mu) - 3d0*im*eta
                 
                 U3_2p1h(anu) = U3_2p1h(anu) - num/dem1/dem2
                 
                 num = 2d0*sqrt(2d0)*rho(k,c,nu)*rho(k,a,mu)*rho(c,p,mu)
                 dem1 = eHF(c) - eHF(k) - Om(nu) - im*eta
                 dem2 = eHF(a) - eHF(k) - Om(mu) - im*eta

                 U3_2p1h(anu) = U3_2p1h(anu) + 0.5d0*num/dem1/dem2
                 
                 num = 2d0*sqrt(2d0)*rho(a,c,mu)*rho(c,k,nu)*rho(p,k,mu)
                 dem1 = eHF(a) - eHF(k) + Om(mu) + Om(nu) - 4d0*im*eta
                 dem2 = eHF(c) - eHF(k) + Om(nu) - 3d0*im*eta
                 
                 U3_2p1h(anu) = U3_2p1h(anu) - num/dem1/dem2
                 
              end do
           end do
        end do
        
        do i=nC+1,nO
           do j=nC+1,nO
              do mu=1,nS
                 
                 num = 2d0*sqrt(2d0)*rho(a,j,mu)*rho(j,i,nu)*rho(p,i,mu)
                 dem1 = eHF(a) - eHF(i) + Om(mu) + Om(nu) - 4d0*im*eta
                 dem2 = eHF(a) - eHF(j) + Om(mu) - 3d0*im*eta
                 
                 U3_2p1h(anu) = U3_2p1h(anu) + num/dem1/dem2
                 
              end do
           end do
        end do
        
     end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  
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
        
        dem = w - eHF(i) + Om(nu) - 2d0*im*eta
        K_2h1p(inu)  =   1d0/dem
        dK_2h1p(inu) = - 1d0/dem/dem
        
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
        
        dem = w - eHF(a) - Om(nu) + 2d0*im*eta
        K_2p1h(anu)  =   1d0/dem
        dK_2p1h(anu) = - 1d0/dem/dem
        
     end do
  end do
  
  !----------------------------!
  !    Building self-energy    !
  !----------------------------!

  !$OMP PARALLEL DEFAULT(NONE)                &
  !$OMP          PRIVATE(i,j, nu,mu, inu,jmu) &
  !$OMP          SHARED(nC, nO, nOrb, nR, nS, K_2h1p,dK_2h1p, U1_2h1p,U2_2h1p,U3_2h1p, C1_2h1p) &
  !$OMP REDUCTION(+:Re_SigC,Im_SigC,Z)
  !$OMP DO COLLAPSE(2)
  do i=nC+1,nO
     do nu=1,nS
        inu = (i - nC - 1) * nS + nu
             
        Re_SigC = Re_SigC + Real(U1_2h1p(inu) *  K_2h1p(inu) * U1_2h1p(inu))
        Im_SigC = Im_SigC + Aimag(U1_2h1p(inu) *  K_2h1p(inu) * U1_2h1p(inu))
        Z    = Z    + Real(U1_2h1p(inu) * dK_2h1p(inu) * U1_2h1p(inu))
        
        Re_SigC = Re_SigC + Real(2d0*U2_2h1p(inu) *  K_2h1p(inu) * U1_2h1p(inu))
        Im_SigC = Im_SigC + Aimag(2d0*U2_2h1p(inu) *  K_2h1p(inu) * U1_2h1p(inu))
        Z    = Z    + Real(2d0*U2_2h1p(inu) * dK_2h1p(inu) * U1_2h1p(inu))
             
        ! SigC(p) = SigC(p) + U1_2h1p(inu) *  K_2h1p(inu) * U2_2h1p(inu)
        ! Z(p)    = Z(p)    + U1_2h1p(inu) * dK_2h1p(inu) * U2_2h1p(inu)
        
        Re_SigC = Re_SigC + Real(2d0*U3_2h1p(inu) *  K_2h1p(inu) * U1_2h1p(inu))
        Im_SigC = Im_SigC + Aimag(2d0*U3_2h1p(inu) *  K_2h1p(inu) * U1_2h1p(inu))
        Z    = Z    + Real(2d0*U3_2h1p(inu) * dK_2h1p(inu) * U1_2h1p(inu))
        
        ! SigC(p) = SigC(p) + U1_2h1p(inu) *  K_2h1p(inu) * U3_2h1p(inu)
        ! Z(p)    = Z(p)    + U1_2h1p(inu) * dK_2h1p(inu) * U3_2h1p(inu)

        jmu = 0
        do j=nC+1,nO
           do mu=1,nS
              jmu = jmu + 1
              
              Re_SigC = Re_SigC + Real(U1_2h1p(inu) *  K_2h1p(inu) * C1_2h1p(inu,jmu) *  K_2h1p(jmu) * U1_2h1p(jmu))
              Im_SigC = Im_SigC + Aimag(U1_2h1p(inu) *  K_2h1p(inu) * C1_2h1p(inu,jmu) *  K_2h1p(jmu) * U1_2h1p(jmu))
              Z    = Z    + Real(2d0*U1_2h1p(inu) * dK_2h1p(inu) * C1_2h1p(inu,jmu) *  K_2h1p(jmu) * U1_2h1p(jmu))
              ! Z(p)    = Z(p)    + U1_2h1p(inu) *  K_2h1p(inu) * C1_2h1p(inu,jmu) * dK_2h1p(jmu) * U1_2h1p(jmu)
              
           end do
        end do
        
     end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  !$OMP PARALLEL DEFAULT(NONE)                &
  !$OMP          PRIVATE(a,b, nu,mu, anu,bmu) &
  !$OMP          SHARED(nC, nO, nOrb, nR, nS, K_2p1h,dK_2p1h, U1_2p1h,U2_2p1h,U3_2p1h, C1_2p1h) &
  !$OMP REDUCTION(+:Re_SigC,Im_SigC,Z)
  !$OMP DO COLLAPSE(2)
  do a=nO+1,nOrb-nR
     do nu=1,nS
        anu = (a - nO - 1) * nS + nu
        
        Re_SigC = Re_SigC + Real(U1_2p1h(anu) *  K_2p1h(anu) * U1_2p1h(anu))
        Im_SigC = Im_SigC + Aimag(U1_2p1h(anu) *  K_2p1h(anu) * U1_2p1h(anu))
        Z       = Z       + Real(U1_2p1h(anu) * dK_2p1h(anu) * U1_2p1h(anu))
        
        Re_SigC = Re_SigC + Real(2d0*U2_2p1h(anu) *  K_2p1h(anu) * U1_2p1h(anu))
        Im_SigC = Im_SigC + Aimag(2d0*U2_2p1h(anu) *  K_2p1h(anu) * U1_2p1h(anu))
        Z       = Z       + Real(2d0*U2_2p1h(anu) * dK_2p1h(anu) * U1_2p1h(anu))
        
        ! SigC(p) = SigC(p) + U1_2p1h(anu) *  K_2p1h(anu) * U2_2p1h(anu)
        ! Z(p)    = Z(p)    + U1_2p1h(anu) * dK_2p1h(anu) * U2_2p1h(anu)
        
        Re_SigC = Re_SigC + Real(2d0*U3_2p1h(anu) *  K_2p1h(anu) * U1_2p1h(anu))
        Im_SigC = Im_SigC + Aimag(2d0*U3_2p1h(anu) *  K_2p1h(anu) * U1_2p1h(anu))
        Z       = Z       + Real(2d0*U3_2p1h(anu) * dK_2p1h(anu) * U1_2p1h(anu))
        
        ! SigC(p) = SigC(p) + U1_2p1h(anu) *  K_2p1h(anu) * U3_2p1h(anu)
        ! Z(p)    = Z(p)    + U1_2p1h(anu) * dK_2p1h(anu) * U3_2p1h(anu)
        
        bmu = 0
        do b=nO+1,nOrb-nR
           do mu=1,nS
              bmu = bmu + 1
              
              Re_SigC = Re_SigC + Real(U1_2p1h(anu) *  K_2p1h(anu) * C1_2p1h(anu,bmu) *  K_2p1h(bmu) * U1_2p1h(bmu))
              Im_SigC = Im_SigC + Aimag(U1_2p1h(anu) *  K_2p1h(anu) * C1_2p1h(anu,bmu) *  K_2p1h(bmu) * U1_2p1h(bmu))
              Z       = Z       + Real(2d0*U1_2p1h(anu) * dK_2p1h(anu) * C1_2p1h(anu,bmu) *  K_2p1h(bmu) * U1_2p1h(bmu))
              ! Z(p)    = Z(p)    + U1_2p1h(anu) *  K_2p1h(anu) * C1_2p1h(anu,bmu) * dK_2p1h(bmu) * U1_2p1h(bmu)
              
           end do
        end do
          
     end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  
  
  if (do_psd) then 
     inu = 0
     do i=nC+1,nO
        do nu=1,nS
           inu = inu + 1
           
           Re_SigC = Re_SigC + Real(U2_2h1p(inu) *  K_2h1p(inu) * U2_2h1p(inu))
           Im_SigC = Im_SigC + Aimag(U2_2h1p(inu) *  K_2h1p(inu) * U2_2h1p(inu))
           Z       = Z       + Real(U2_2h1p(inu) * dK_2h1p(inu) * U2_2h1p(inu))
           
           Re_SigC = Re_SigC + Real(U3_2h1p(inu) *  K_2h1p(inu) * U2_2h1p(inu))
           Im_SigC = Im_SigC + Aimag(U3_2h1p(inu) *  K_2h1p(inu) * U2_2h1p(inu))
           Z       = Z       + Real(U3_2h1p(inu) * dK_2h1p(inu) * U2_2h1p(inu))
           
           Re_SigC = Re_SigC + Real(U2_2h1p(inu) *  K_2h1p(inu) * U3_2h1p(inu))
           Im_SigC = Im_SigC + Aimag(U2_2h1p(inu) *  K_2h1p(inu) * U3_2h1p(inu))
           Z       = Z       + Real(U2_2h1p(inu) * dK_2h1p(inu) * U3_2h1p(inu))
           
           Re_SigC = Re_SigC + Real(U3_2h1p(inu) *  K_2h1p(inu) * U3_2h1p(inu))
           Im_SigC = Im_SigC + Aimag(U3_2h1p(inu) *  K_2h1p(inu) * U3_2h1p(inu))
           Z       = Z       + Real(U3_2h1p(inu) * dK_2h1p(inu) * U3_2h1p(inu))
           
           jmu = 0
           do j=nC+1,nO
              do mu=1,nS
                 jmu = jmu + 1
                 
                 Re_SigC = Re_SigC + Real(U2_2h1p(inu) *  K_2h1p(inu) * C1_2h1p(inu,jmu) *  K_2h1p(jmu) * U1_2h1p(jmu))
                 Im_SigC = Im_SigC + Aimag(U2_2h1p(inu) *  K_2h1p(inu) * C1_2h1p(inu,jmu) *  K_2h1p(jmu) * U1_2h1p(jmu))
                 Z    = Z    + Real(U2_2h1p(inu) * dK_2h1p(inu) * C1_2h1p(inu,jmu) *  K_2h1p(jmu) * U1_2h1p(jmu))
                 Z    = Z    + Real(U2_2h1p(inu) *  K_2h1p(inu) * C1_2h1p(inu,jmu) * dK_2h1p(jmu) * U1_2h1p(jmu))
                 
                 Re_SigC = Re_SigC + Real(U1_2h1p(inu) *  K_2h1p(inu) * C1_2h1p(inu,jmu) *  K_2h1p(jmu) * U2_2h1p(jmu))
                 Im_SigC = Im_SigC + Aimag(U1_2h1p(inu) *  K_2h1p(inu) * C1_2h1p(inu,jmu) *  K_2h1p(jmu) * U2_2h1p(jmu))
                 Z    = Z    + Real(U1_2h1p(inu) * dK_2h1p(inu) * C1_2h1p(inu,jmu) *  K_2h1p(jmu) * U2_2h1p(jmu))
                 Z    = Z    + Real(U1_2h1p(inu) *  K_2h1p(inu) * C1_2h1p(inu,jmu) * dK_2h1p(jmu) * U2_2h1p(jmu))
                 
                 Re_SigC = Re_SigC + Real(U3_2h1p(inu) *  K_2h1p(inu) * C1_2h1p(inu,jmu) *  K_2h1p(jmu) * U1_2h1p(jmu))
                 Im_SigC = Im_SigC + Aimag(U3_2h1p(inu) *  K_2h1p(inu) * C1_2h1p(inu,jmu) *  K_2h1p(jmu) * U1_2h1p(jmu))
                 Z    = Z    + Real(U3_2h1p(inu) * dK_2h1p(inu) * C1_2h1p(inu,jmu) *  K_2h1p(jmu) * U1_2h1p(jmu))
                 Z    = Z    + Real(U3_2h1p(inu) *  K_2h1p(inu) * C1_2h1p(inu,jmu) * dK_2h1p(jmu) * U1_2h1p(jmu))
                 
                 Re_SigC = Re_SigC + Real(U1_2h1p(inu) *  K_2h1p(inu) * C1_2h1p(inu,jmu) *  K_2h1p(jmu) * U3_2h1p(jmu))
                 Im_SigC = Im_SigC + Aimag(U1_2h1p(inu) *  K_2h1p(inu) * C1_2h1p(inu,jmu) *  K_2h1p(jmu) * U3_2h1p(jmu))
                 Z    = Z    + Real(U1_2h1p(inu) * dK_2h1p(inu) * C1_2h1p(inu,jmu) *  K_2h1p(jmu) * U3_2h1p(jmu))
                 Z    = Z    + Real(U1_2h1p(inu) *  K_2h1p(inu) * C1_2h1p(inu,jmu) * dK_2h1p(jmu) * U3_2h1p(jmu))

                 kla = 0
                 do k=nC+1,nO
                    do la=1,nS
                       kla = kla + 1

                       Re_SigC = Re_SigC + Real(U1_2h1p(inu) *  K_2h1p(inu) * C1_2h1p(inu,jmu) *  K_2h1p(jmu) * C1_2h1p(jmu,kla) *  K_2h1p(kla) * U1_2h1p(kla))
                       Im_SigC = Im_SigC + Aimag(U1_2h1p(inu) *  K_2h1p(inu) * C1_2h1p(inu,jmu) *  K_2h1p(jmu) * C1_2h1p(jmu,kla) *  K_2h1p(kla) * U1_2h1p(kla))
                       Z    = Z    + Real(3d0*U1_2h1p(inu) *  dK_2h1p(inu) * C1_2h1p(inu,jmu) *  K_2h1p(jmu)  * C1_2h1p(jmu,kla) *  K_2h1p(kla)  * U1_2h1p(kla))
                       ! Z    = Z    + U1_2h1p(inu) *  K_2h1p(inu)  * C1_2h1p(inu,jmu) *  dK_2h1p(jmu) * C1_2h1p(jmu,kla) *  K_2h1p(kla)  * U1_2h1p(kla)
                       ! Z    = Z    + U1_2h1p(inu) *  K_2h1p(inu)  * C1_2h1p(inu,jmu) *  K_2h1p(jmu)  * C1_2h1p(jmu,kla) *  dK_2h1p(kla) * U1_2h1p(kla)
                       
                    end do
                 end do
                 
              end do
           end do
           
        end do
     end do
     ! !$OMP END DO
     ! !$OMP END PARALLEL
     
     
     !$OMP PARALLEL DEFAULT(NONE)                &
     !$OMP          PRIVATE(a,b,c, nu,mu,la, anu,bmu,cla) &
     !$OMP          SHARED(nC, nO, nOrb, nR, nS, K_2p1h,dK_2p1h, U1_2p1h,U2_2p1h,U3_2p1h, C1_2p1h) &
     !$OMP REDUCTION(+:Re_SigC,Im_SigC,Z)
     !$OMP DO COLLAPSE(2)
     do a=nO+1,nOrb-nR
        do nu=1,nS
           anu = (a - nO - 1) * nS + nu           
           
           Re_SigC = Re_SigC + Real(U2_2p1h(anu) *  K_2p1h(anu) * U2_2p1h(anu))
           Im_SigC = Im_SigC + Aimag(U2_2p1h(anu) *  K_2p1h(anu) * U2_2p1h(anu))
           Z    = Z    + Real(U2_2p1h(anu) * dK_2p1h(anu) * U2_2p1h(anu))
           
           Re_SigC = Re_SigC + Real(U3_2p1h(anu) *  K_2p1h(anu) * U2_2p1h(anu))
           Im_SigC = Im_SigC + Aimag(U3_2p1h(anu) *  K_2p1h(anu) * U2_2p1h(anu))
           Z    = Z    + Real(U3_2p1h(anu) * dK_2p1h(anu) * U2_2p1h(anu))
             
           Re_SigC = Re_SigC + Real(U2_2p1h(anu) *  K_2p1h(anu) * U3_2p1h(anu))
           Im_SigC = Im_SigC + Aimag(U2_2p1h(anu) *  K_2p1h(anu) * U3_2p1h(anu))
           Z    = Z    + Real(U2_2p1h(anu) * dK_2p1h(anu) * U3_2p1h(anu))
           
           Re_SigC = Re_SigC + Real(U3_2p1h(anu) *  K_2p1h(anu) * U3_2p1h(anu))
           Im_SigC = Im_SigC + Aimag(U3_2p1h(anu) *  K_2p1h(anu) * U3_2p1h(anu))
           Z    = Z    + Real(U3_2p1h(anu) * dK_2p1h(anu) * U3_2p1h(anu))
           
           bmu = 0
           do b=nO+1,nOrb-nR
              do mu=1,nS
                 bmu = bmu + 1
                 
                 Re_SigC = Re_SigC + Real(U2_2p1h(anu) *  K_2p1h(anu) * C1_2p1h(anu,bmu) *  K_2p1h(bmu) * U1_2p1h(bmu))
                 Im_SigC = Im_SigC + Aimag(U2_2p1h(anu) *  K_2p1h(anu) * C1_2p1h(anu,bmu) *  K_2p1h(bmu) * U1_2p1h(bmu))
                 Z       = Z    + Real(U2_2p1h(anu) * dK_2p1h(anu) * C1_2p1h(anu,bmu) *  K_2p1h(bmu) * U1_2p1h(bmu))
                 Z       = Z    + Real(U2_2p1h(anu) *  K_2p1h(anu) * C1_2p1h(anu,bmu) * dK_2p1h(bmu) * U1_2p1h(bmu))
                 
                 Re_SigC = Re_SigC + Real(U1_2p1h(anu) *  K_2p1h(anu) * C1_2p1h(anu,bmu) *  K_2p1h(bmu) * U2_2p1h(bmu))
                 Im_SigC = Im_SigC + Aimag(U1_2p1h(anu) *  K_2p1h(anu) * C1_2p1h(anu,bmu) *  K_2p1h(bmu) * U2_2p1h(bmu))
                 Z       = Z    + Real(U1_2p1h(anu) * dK_2p1h(anu) * C1_2p1h(anu,bmu) *  K_2p1h(bmu) * U2_2p1h(bmu))
                 Z       = Z    + Real(U1_2p1h(anu) *  K_2p1h(anu) * C1_2p1h(anu,bmu) * dK_2p1h(bmu) * U2_2p1h(bmu))
                
                 Re_SigC = Re_SigC + Real(U3_2p1h(anu) *  K_2p1h(anu) * C1_2p1h(anu,bmu) *  K_2p1h(bmu) * U1_2p1h(bmu))
                 Im_SigC = Im_SigC + Aimag(U3_2p1h(anu) *  K_2p1h(anu) * C1_2p1h(anu,bmu) *  K_2p1h(bmu) * U1_2p1h(bmu))
                 Z       = Z    + Real(U3_2p1h(anu) * dK_2p1h(anu) * C1_2p1h(anu,bmu) *  K_2p1h(bmu) * U1_2p1h(bmu))
                 Z       = Z    + Real(U3_2p1h(anu) *  K_2p1h(anu) * C1_2p1h(anu,bmu) * dK_2p1h(bmu) * U1_2p1h(bmu))
                   
                 Re_SigC = Re_SigC + Real(U1_2p1h(anu) *  K_2p1h(anu) * C1_2p1h(anu,bmu) *  K_2p1h(bmu) * U3_2p1h(bmu))
                 Im_SigC = Im_SigC + Aimag(U1_2p1h(anu) *  K_2p1h(anu) * C1_2p1h(anu,bmu) *  K_2p1h(bmu) * U3_2p1h(bmu))
                 Z       = Z    + Real(U1_2p1h(anu) * dK_2p1h(anu) * C1_2p1h(anu,bmu) *  K_2p1h(bmu) * U3_2p1h(bmu))
                 Z       = Z    + Real(U1_2p1h(anu) *  K_2p1h(anu) * C1_2p1h(anu,bmu) * dK_2p1h(bmu) * U3_2p1h(bmu))

                 
                 cla = 0
                 do c=nO+1,nOrb-nR
                    do la=1,nS
                       cla = cla + 1

                       Re_SigC = Re_SigC + Real(U1_2p1h(anu) *  K_2p1h(anu) * C1_2p1h(anu,bmu) *  K_2p1h(bmu) * C1_2p1h(bmu,cla) *  K_2p1h(cla) * U1_2p1h(cla))
                       Im_SigC = Im_SigC + Aimag(U1_2p1h(anu) *  K_2p1h(anu) * C1_2p1h(anu,bmu) *  K_2p1h(bmu) * C1_2p1h(bmu,cla) *  K_2p1h(cla) * U1_2p1h(cla))
                       Z    = Z    + Real(3d0* U1_2p1h(anu) * dK_2p1h(anu) * C1_2p1h(anu,bmu) *  K_2p1h(bmu) * C1_2p1h(bmu,cla) *  K_2p1h(cla) * U1_2p1h(cla))
                       
                    end do
                 end do
                 
                 
              end do
           end do
           
        end do
     end do
     !$OMP END DO
     !$OMP END PARALLEL

  end if
  
end subroutine
