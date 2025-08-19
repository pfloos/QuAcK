subroutine G_Parquet_self_energy(eta,nOrb,nC,nO,nV,nR,nS,nOO,nVV,eQP,ERI,&
                                 eh_rho,eh_Om,ee_rho,ee_Om,hh_rho,hh_Om,SigC,Z)

! Compute correlation part of the self-energy coming from irreducible vertices contribution

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC,nO,nV,nR
  integer,intent(in)            :: nS,nOO,nVV
  double precision,intent(in)   :: eQP(nOrb)
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: eh_rho(nOrb,nOrb,nS)
  double precision,intent(in)   :: eh_Om(nS)
  double precision,intent(in)   :: ee_rho(nOrb,nOrb,nVV)
  double precision,intent(in)   :: ee_Om(nVV)
  double precision,intent(in)   :: hh_rho(nOrb,nOrb,nOO)
  double precision,intent(in)   :: hh_Om(nOO)

! Local variables
  integer                       :: i,j,k,a,b,c
  integer                       :: p,q,n
  double precision              :: eps_p,eps_q,eps_pq,reg
  double precision              :: dem1,dem2,reg1,reg2
  double precision              :: dem1_p,dem1_q,dem1_pq
  double precision              :: dem2_p,dem2_q,dem2_pq
  double precision              :: num
  double precision              :: start_t,end_t,t

  logical                       :: print_self_energy

! Output variables

  double precision,intent(out)  :: SigC(nOrb,nOrb)
  double precision,intent(out)  :: Z(nOrb)

! Initialize 

  SigC(:,:) = 0d0

! Memory allocation for self-energy decomposition
 
!-----------------------------------!
! 2nd-order part of the self-energy !
!-----------------------------------!

  call wall_time(start_t)

  do q=nC+1,nOrb-nR
     do p=nC+1,nOrb-nR

        ! 2h1p 
        do i=nC+1,nO
           do j=nC+1,nO
              do a=nO+1,nOrb-nR

                 eps_p = eQP(p) + eQP(a) - eQP(i) - eQP(j)
                 eps_q = eQP(q) + eQP(a) - eQP(i) - eQP(j)
                 eps_pq = eps_p * eps_p + eps_q * eps_q
                 reg = 1d0 - exp(- eta * eps_pq)
                 num = 0.5d0*(ERI(p,a,j,i) - ERI(p,a,i,j))*(ERI(q,a,j,i) - ERI(q,a,i,j))

                 SigC(p,q) = SigC(p,q) + num*reg*(eps_p + eps_q)/eps_pq
                 
              end do
           end do
        end do
     
        ! 2p1h
        do i=nC+1,nO
           do a=nO+1,nOrb-nR
              do b=nO+1,nOrb-nR
                 
                 eps_p = eQP(p) + eQP(i) - eQP(a) - eQP(b)
                 eps_q = eQP(q) + eQP(i) - eQP(a) - eQP(b)
                 eps_pq = eps_p * eps_p + eps_q * eps_q
                 reg = 1d0 - exp(- eta * eps_pq)
                 num = 0.5d0*(ERI(p,i,b,a) - ERI(p,i,a,b))*(ERI(q,i,b,a) - ERI(q,i,a,b))
                 
                 SigC(p,q) = SigC(p,q) + num*reg*(eps_p + eps_q)/eps_pq
                 
              end do
           end do
        end do
     end do
  end do
     
  call wall_time(end_t)
  t = end_t - start_t
  
  write(*,'(1X,A50,1X,F9.3,A8)') 'Wall time for building GF(2) self-energy =',t,' seconds'
  write(*,*)
  
!-----------------------------!
!  eh part of the self-energy !
!-----------------------------!
     
  call wall_time(start_t)
  
  !$OMP PARALLEL DEFAULT(NONE)    &
  !$OMP PRIVATE(p,q,i,a,j,b,n,num,dem1,dem1_p,dem1_q,dem1_pq,dem2_p,dem2_q,dem2_pq,reg1,reg2) &
  !$OMP SHARED(nC,nO,nOrb,nR,nS,eta,ERI,eQP,eh_rho,eh_Om,SigC)
  !$OMP DO
  do q=nC+1,nOrb-nR
     do p=nC+1,nOrb-nR
     
        do i=nC+1,nO
           do a=nO+1,nOrb-nR
           
              do n=1,nS

                 do j=nC+1,nO

                    ! 2h1p(d) * 2h1p 
                    num  = (ERI(p,a,i,j) - ERI(p,a,j,i)) * eh_rho(i,a,n) * eh_rho(q,j,n)
                    dem1_p = eQP(p) - eQP(j) + eh_Om(n)
                    dem1_q = eQP(q) - eQP(j) + eh_Om(n)
                    dem1_pq = dem1_p * dem1_p + dem1_q * dem1_q
                    reg1 = 1d0 - exp(- eta * dem1_pq)
                    dem2_p = eQP(p) - eQP(i) - eQP(j) + eQP(a)
                    dem2_q = eQP(q) - eQP(i) - eQP(j) + eQP(a)
                    dem2_pq = dem2_p * dem2_p + dem2_q * dem2_q
                    reg2 = 1d0 - exp(- eta * dem2_pq)
                    
                    SigC(p,q) = SigC(p,q) + num * reg1 * (dem1_p + dem1_q)/dem1_pq * reg2 * (dem2_p + dem2_q)/dem2_pq

                    
                    ! 2h2p(d) * 2h1p(d) 
                    num  = (ERI(p,i,a,j) - ERI(p,i,j,a)) * eh_rho(a,i,n) * eh_rho(q,j,n)
                    dem1 = eQP(a) - eQP(i) + eh_Om(n) 
                    reg1 = 1d0 - exp(- 2d0 * eta * dem1 * dem1)
                    dem2_p = eQP(p) - eQP(j) + eh_Om(n)
                    dem2_q = eQP(q) - eQP(j) + eh_Om(n)
                    dem2_pq = dem2_p * dem2_p + dem2_q * dem2_q
                    reg2 = 1d0 - exp(- eta * dem2_pq)
                    
                    SigC(p,q) = SigC(p,q) + num * (reg1/dem1) * reg2 * (dem2_p + dem2_q)/dem2_pq

                    
                    ! 2h2p(d) * 2h1p 
                    num  = (ERI(p,a,i,j) - ERI(p,a,j,i)) * eh_rho(a,i,n) * eh_rho(j,q,n) 
                    dem1 = eQP(a) - eQP(i) + eh_Om(n) 
                    reg1 = 1d0 - exp(- 2d0 * eta * dem1 * dem1)
                    dem2_p = eQP(p) - eQP(i) - eQP(j) + eQP(a)
                    dem2_q = eQP(q) - eQP(i) - eQP(j) + eQP(a)
                    dem2_pq = dem2_p * dem2_p + dem2_q * dem2_q
                    reg2 = 1d0 - exp(- eta * dem2_pq)
                    
                    SigC(p,q) = SigC(p,q) + num * (reg1/dem1) * reg2 * (dem2_p + dem2_q)/dem2_pq
                    
                 end do ! j
                 
                 do b=nO+1,nOrb-nR
                    
                    ! 2p1h(d) * 2p1h
                    num  = - (ERI(p,i,a,b) - ERI(p,i,b,a)) * eh_rho(i,a,n) * eh_rho(b,q,n) 
                    dem1_p = eQP(p) - eQP(b) - eh_Om(n)
                    dem1_q = eQP(q) - eQP(b) - eh_Om(n)
                    dem1_pq = dem1_p * dem1_p + dem1_q * dem1_q
                    reg1 = 1d0 - exp(- eta * dem1_pq)
                    dem2_p = eQP(p) + eQP(i) - eQP(a) - eQP(b)
                    dem2_q = eQP(q) + eQP(i) - eQP(a) - eQP(b)
                    dem2_pq = dem2_p * dem2_p + dem2_q * dem2_q
                    reg2 = 1d0 - exp(- eta * dem2_pq)
                    
                    SigC(p,q) = SigC(p,q) + num * reg1 * (dem1_p + dem1_q)/dem1_pq * reg2 * (dem2_p + dem2_q)/dem2_pq
                    
                    
                    ! 2h2p(d) * 2p1h(d)
                    num  = (ERI(p,a,i,b) - ERI(p,a,b,i)) * eh_rho(a,i,n) * eh_rho(b,q,n) 
                    dem1 = eQP(a) - eQP(i) + eh_Om(n)
                    reg1 = 1d0 - exp(- 2d0 * eta * dem1 * dem1)
                    dem2_p = eQP(p) - eQP(b) - eh_Om(n)
                    dem2_q = eQP(q) - eQP(b) - eh_Om(n)
                    dem2_pq = dem2_p * dem2_p + dem2_q * dem2_q
                    reg2 = 1d0 - exp(- eta * dem2_pq)
                    
                    SigC(p,q) = SigC(p,q) + num * (reg1/dem1) * reg2 * (dem2_p + dem2_q)/dem2_pq
                    
                    
                    ! 2h2p(d) * 2p1h
                    num  = (ERI(p,i,a,b) - ERI(p,i,b,a)) * eh_rho(a,i,n) * eh_rho(q,b,n) 
                    dem1 = eQP(a) - eQP(i) + eh_Om(n) 
                    reg1 = 1d0 - exp(- 2d0 * eta * dem1 * dem1)
                    dem2_p = eQP(p) + eQP(i) - eQP(a) - eQP(b)
                    dem2_q = eQP(q) + eQP(i) - eQP(a) - eQP(b)
                    dem2_pq = dem2_p * dem2_p + dem2_q * dem2_q
                    reg2 = 1d0 - exp(- eta * dem2_pq)
                    
                    SigC(p,q) = SigC(p,q) + num * (reg1/dem1) * reg2 * (dem2_p + dem2_q)/dem2_pq
                    
                 end do ! b
                 
              end do ! n
              
           end do ! a
        end do ! i
        
     end do ! p
  end do ! q
  !$OMP END DO
  !$OMP END PARALLEL

  call wall_time(end_t)
  t = end_t - start_t

  write(*,'(1X,A50,1X,F9.3,A8)') 'Wall time for building eh self-energy =',t,' seconds'
  write(*,*) 

!-----------------------------!
!  pp part of the self-energy !
!-----------------------------!

  call wall_time(start_t)

  !$OMP PARALLEL DEFAULT(NONE)    &
  !$OMP PRIVATE(p,q,i,j,k,c,n,num,dem1,dem1_p,dem1_q,dem1_pq,dem2,dem2_p,dem2_q,dem2_pq,reg1,reg2) &
  !$OMP SHARED(nC,nO,nOrb,nR,nOO,nVV,eta,ERI,eQP,ee_rho,ee_Om,hh_rho,hh_Om,SigC)
  !$OMP DO
  do q=nC+1,nOrb-nR
     do p=nC+1,nOrb-nR
     
        do i=nC+1,nO
           do j=nC+1,nO
              do n=1,nVV
                 
                 do k=nC+1,nO
                    
                    ! 2h2p(d) * 2p1h(d)
                    num  = - ERI(p,k,i,j) * ee_rho(i,j,n) * ee_rho(q,k,n)
                    dem1 = ee_Om(n) - eQP(i) - eQP(j)
                    reg1 = 1d0 - exp(- 2d0 * eta * dem1 * dem1)
                    dem2_p = eQP(p) + eQP(k) - ee_Om(n)
                    dem2_q = eQP(q) + eQP(k) - ee_Om(n)
                    dem2_pq = dem2_p * dem2_p + dem2_q * dem2_q
                    reg2 = 1d0 - exp(- eta * dem2_pq)
                    
                    SigC(p,q) = SigC(p,q) + num * (reg1/dem1) * reg2 * (dem2_p + dem2_q)/dem2_pq
                    
                 end do ! k
                 
                 do c=nO+1,nOrb-nR
                    
                    ! 2h2p(d) * 2h1p
                    num  = - ERI(p,c,i,j) * ee_rho(i,j,n) * ee_rho(q,c,n)
                    dem1 = ee_Om(n) - eQP(i) - eQP(j)
                    reg1 = 1d0 - exp(- 2d0 * eta * dem1 * dem1)
                    dem2_p = eQP(p) + eQP(c) - eQP(i) - eQP(j)
                    dem2_q = eQP(q) + eQP(c) - eQP(i) - eQP(j)
                    dem2_pq = dem2_p * dem2_p + dem2_q * dem2_q
                    reg2 = 1d0 - exp(- eta * dem2_pq)
                    
                    SigC(p,q) = SigC(p,q) + num * (reg1/dem1) * reg2 * (dem2_p + dem2_q)/dem2_pq
                    
                 end do ! a
              end do ! n
              do n=1,nOO
                 
                 do c=nO+1,nOrb-nR
                    
                    ! 2h1p(d) * 2h1p
                    num  = - ERI(p,c,i,j) * hh_rho(i,j,n) * hh_rho(q,c,n)
                    dem1_p = eQP(p) + eQP(c) - hh_Om(n)
                    dem1_q = eQP(q) + eQP(c) - hh_Om(n)
                    dem1_pq = dem1_p * dem1_p + dem1_q * dem1_q
                    reg1 = 1d0 - exp(- eta * dem1_pq)
                    dem2_p = eQP(p) + eQP(c) - eQP(i) - eQP(j)
                    dem2_q = eQP(q) + eQP(c) - eQP(i) - eQP(j)
                    dem2_pq = dem2_p * dem2_p + dem2_q * dem2_q
                    reg2 = 1d0 - exp(- eta * dem2_pq)
                    
                    SigC(p,q) = SigC(p,q) + num * reg1 * (dem1_p + dem1_q)/dem1_pq * reg2 * (dem2_p + dem2_q)/dem2_pq
                    
                 end do ! c
              end do ! n
           end do ! j
        end do ! i
        
     end do ! p
  end do ! q
  !$OMP END DO
  !$OMP END PARALLEL
  !$OMP PARALLEL DEFAULT(NONE)    &
  !$OMP PRIVATE(p,q,k,a,b,c,n,num,dem1,dem1_p,dem1_q,dem1_pq,dem2,dem2_p,dem2_q,dem2_pq,reg1,reg2) &
  !$OMP SHARED(nC,nO,nOrb,nR,nOO,nVV,eta,ERI,eQP,ee_rho,ee_Om,hh_rho,hh_Om,SigC)
  !$OMP DO
  do q=nC+1,nOrb-nR
     do p=nC+1,nOrb-nR
        do a=nO+1,nOrb-nR
           do b=nO+1,nOrb-nR
              do n=1,nOO
                 
                 do c=nO+1,nOrb-nR
                    
                    ! 2p2h(d) * 2h1p(d)   
                    num  = ERI(p,c,a,b) * hh_rho(a,b,n) * hh_rho(q,c,n)
                    dem1 = hh_Om(n) - eQP(a) - eQP(b)
                    reg1 = 1d0 - exp(- 2d0 * eta * dem1 * dem1)
                    dem2_p = eQP(p) + eQP(c) - hh_Om(n)
                    dem2_q = eQP(q) + eQP(c) - hh_Om(n)
                    dem2_pq = dem2_p * dem2_p + dem2_q * dem2_q
                    reg2 = 1d0 - exp(- eta * dem2_pq)
                    
                    SigC(p,q) = SigC(p,q) + num * (reg1/dem1) * reg2 * (dem2_p + dem2_q)/dem2_pq
                    
                 end do ! c
                 
                 do k=nC+1,nO
                    
                    ! 2p2h(d) * 2p1h
                    num  = ERI(p,k,a,b) * hh_rho(a,b,n) * hh_rho(q,k,n)
                    dem1 = hh_Om(n) - eQP(a) - eQP(b)
                    reg1 = 1d0 - exp(- 2d0 * eta * dem1 * dem1)
                    dem2_p = eQP(p) + eQP(k) - eQP(a) - eQP(b)
                    dem2_q = eQP(q) + eQP(k) - eQP(a) - eQP(b)
                    dem2_pq = dem2_p * dem2_p + dem2_q * dem2_q
                    reg2 = 1d0 - exp(- eta * dem2_pq)
                    
                    SigC(p,q) = SigC(p,q) + num * (reg1/dem1) * reg2 * (dem2_p + dem2_q)/dem2_pq
                    
                 end do ! k
              end do ! n
              do n=1,nVV
                 
                 do k=nC+1,nO
                    
                    ! 2p1h * 2p1h(d)
                    num  = ERI(p,k,a,b) * ee_rho(a,b,n) * ee_rho(q,k,n)
                    dem1_p = eQP(p) + eQP(k) - eQP(a) - eQP(b)
                    dem1_q = eQP(q) + eQP(k) - eQP(a) - eQP(b)
                    dem1_pq = dem1_p * dem1_p + dem1_q * dem1_q
                    reg1 = 1d0 - exp(- eta * dem1_pq)
                    dem2_p = eQP(p) + eQP(k) - ee_Om(n)
                    dem2_q = eQP(q) + eQP(k) - ee_Om(n)
                    dem2_pq = dem2_p * dem2_p + dem2_q * dem2_q
                    reg2 = 1d0 - exp(- eta * dem2_pq)
                    
                    SigC(p,q) = SigC(p,q) + num * reg1 * (dem1_p + dem1_q)/dem1_pq * reg2 * (dem2_p + dem2_q)/dem2_pq
                    
                 end do ! c
              end do ! n
           end do ! b
        end do ! a
        
     end do ! p
  end do ! q
  !$OMP END DO
  !$OMP END PARALLEL

  call wall_time(end_t)
  t = end_t - start_t

  write(*,'(1X,A50,1X,F9.3,A8)') 'Wall time for building pp self-energy =',t,' seconds'
  write(*,*)
 
!------------------------!
! Renormalization factor !
!------------------------!

  Z(:) = 1d0
  
end subroutine 
