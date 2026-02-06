subroutine R_Parquet_self_energy_diag(eta,nOrb,nC,nO,nV,nR,nS,nOOs,nVVs,nOOt,nVVt,eQP,ERI, &
                                      eh_sing_rho,eh_sing_Om,eh_trip_rho,eh_trip_Om,       &
                                      ee_sing_rho,ee_sing_Om,ee_trip_rho,ee_trip_Om,       &
                                      hh_sing_rho,hh_sing_Om,hh_trip_rho,hh_trip_Om,       &
                                      EcGM,SigC,Z)

! Compute correlation part of the self-energy with only irreducible vertices contribution

  implicit none
  include 'parameters.h'

! Input variables
  double precision,intent(in)   :: eta
  integer,intent(in)            :: nOrb,nC,nO,nV,nR
  integer,intent(in)            :: nS,nOOs,nVVs,nOOt,nVVt
  double precision,intent(in)   :: eQP(nOrb)
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: eh_sing_rho(nOrb,nOrb,nS)
  double precision,intent(in)   :: eh_sing_Om(nS)
  double precision,intent(in)   :: eh_trip_rho(nOrb,nOrb,nS)
  double precision,intent(in)   :: eh_trip_Om(nS)
  double precision,intent(in)   :: ee_sing_rho(nOrb,nOrb,nVVs)
  double precision,intent(in)   :: ee_sing_Om(nVVs)
  double precision,intent(in)   :: ee_trip_rho(nOrb,nOrb,nVVt)
  double precision,intent(in)   :: ee_trip_Om(nVVt)
  double precision,intent(in)   :: hh_sing_rho(nOrb,nOrb,nOOs)
  double precision,intent(in)   :: hh_sing_Om(nOOs)
  double precision,intent(in)   :: hh_trip_rho(nOrb,nOrb,nOOt)
  double precision,intent(in)   :: hh_trip_Om(nOOt)
  
! Local variables
  integer                       :: i,j,k,a,b,c
  integer                       :: p,n
  double precision              :: eps,dem1,dem2,reg,reg1,reg2
  double precision              :: num
  double precision              :: start_t,end_t,t

  logical                       :: print_self_energy = .false.

  double precision,allocatable  :: Sig_2d(:)
  double precision,allocatable  :: Sig_2x(:)
  double precision,allocatable  :: Sig_1eh(:)
  double precision,allocatable  :: Sig_3eh(:)
  double precision,allocatable  :: Sig_1pp(:)
  double precision,allocatable  :: Sig_3pp(:)

  double precision,allocatable  :: Z_2d(:)
  double precision,allocatable  :: Z_2x(:)
  double precision,allocatable  :: Z_1eh(:)
  double precision,allocatable  :: Z_3eh(:)
  double precision,allocatable  :: Z_1pp(:)
  double precision,allocatable  :: Z_3pp(:)
  
  logical                       :: do_2d_channel  = .true.
  logical                       :: do_2x_channel  = .true.
  logical                       :: do_1eh_channel = .true.
  logical                       :: do_3eh_channel = .true.
  logical                       :: do_1pp_channel = .true.
  logical                       :: do_3pp_channel = .true.
  double precision              :: Kx = 1d0

! Output variables
  double precision,intent(out)  :: EcGM
  double precision,intent(out)  :: SigC(nOrb)
  double precision,intent(out)  :: Z(nOrb)

! Initialize 

  SigC(:) = 0d0
  Z(:)    = 0d0
  EcGM    = 0d0

! Memory allocation for self-energy decomposition

  allocate(Sig_2d(nOrb))
  allocate(Sig_2x(nOrb))
  allocate(Sig_1eh(nOrb))
  allocate(Sig_3eh(nOrb))
  allocate(Sig_1pp(nOrb))
  allocate(Sig_3pp(nOrb))

  allocate(Z_2d(nOrb))
  allocate(Z_2x(nOrb))
  allocate(Z_1eh(nOrb))
  allocate(Z_3eh(nOrb))
  allocate(Z_1pp(nOrb))
  allocate(Z_3pp(nOrb))
  
  Sig_2x(:) = 0d0
  Sig_2d(:) = 0d0

  Z_2x(:) = 0d0
  Z_2d(:) = 0d0

!-----------------------------------!
! 2nd-order part of the self-energy !
!-----------------------------------!
  call wall_time(start_t)
  do p=nC+1,nOrb-nR
     ! 2h1p sum
     do i=nC+1,nO
        do j=nC+1,nO
           do a=nO+1,nOrb-nR

              eps = eQP(p) + eQP(a) - eQP(i) - eQP(j)
              reg = (1d0 - exp(- 2d0 * eta * eps * eps))
              num = 2d0*ERI(p,a,j,i)*ERI(j,i,p,a)

              Sig_2d(p) = Sig_2d(p) + num*reg/eps
              Z_2d(p)   = Z_2d(p)   - num*reg/eps**2

              num = - ERI(p,a,j,i)*ERI(j,i,a,p)

              Sig_2x(p) = Sig_2x(p) + num*reg/eps
              Z_2x(p)   = Z_2x(p)   - num*reg/eps**2

           end do
        end do
     end do
     ! 2p1h sum
     do i=nC+1,nO
        do a=nO+1,nOrb-nR
           do b=nO+1,nOrb-nR

              eps = eQP(p) + eQP(i) - eQP(a) - eQP(b)
              reg = (1d0 - exp(- 2d0 * eta * eps * eps))
              num = 2d0*ERI(p,i,b,a)*ERI(b,a,p,i)

              Sig_2d(p) = Sig_2d(p) + num*reg/eps
              Z_2d(p)   = Z_2d(p)   - num*reg/eps**2

              num = - ERI(p,i,b,a)*ERI(b,a,i,p)

              Sig_2x(p) = Sig_2x(p) + num*reg/eps
              Z_2x(p)   = Z_2x(p)   - num*reg/eps**2

           end do
        end do
     end do
  end do
  call wall_time(end_t)
  t = end_t - start_t

  write(*,'(1X,A50,1X,F9.3,A8)') 'Wall time for 2nd-order self-energy =',t,' seconds'
  write(*,*)

! Self-energy decomposition

  SigC(:) = Sig_2d(:) + Sig_2x(:)
  Z(:)    = Z_2d(:) + Z_2x(:)

!-------------------------------------!
!  singlet eh part of the self-energy !
!-------------------------------------!
  call wall_time(start_t)
  !$OMP PARALLEL DEFAULT(NONE)    &
  !$OMP PRIVATE(p,i,a,j,b,n,num,dem1,dem2,reg1,reg2) &
  !$OMP SHARED(nC,nO,nOrb,nR,nS,eta,ERI,eQP,eh_sing_rho,eh_sing_Om,SigC,Z,Kx)
  !$OMP DO
  do p=nC+1,nOrb-nR
     
     do i=nC+1,nO
        do a=nO+1,nOrb-nR
           
           do n=1,nS
              !3h2p
              do j=nC+1,nO
                 
                 num  = (Kx*0.5d0*ERI(p,a,i,j) - ERI(p,a,j,i))* &
                      eh_sing_rho(i,a,n) * eh_sing_rho(p,j,n)
                 
                 dem1 = eQP(a) - eQP(i) - eh_sing_Om(n)
                 reg1 = 1d0 - exp(- 2d0 * eta * dem1 * dem1)
                 dem2 = eQP(p) - eQP(j) + eh_sing_Om(n)
                 reg2 = 1d0 - exp(- 2d0 * eta * dem2 * dem2)
                 
                 SigC(p) = SigC(p) + num * (reg1/dem1) * (reg2/dem2)
                 Z(p)    = Z(p)    - num * (reg1/dem1) * (reg2/dem2/dem2)
                 
                 num  = - (Kx*0.5d0*ERI(p,a,i,j) - ERI(p,a,j,i))* &
                      eh_sing_rho(i,a,n) * eh_sing_rho(p,j,n)
                 
                 dem1 = eQP(a) - eQP(i) - eh_sing_Om(n) 
                 reg1 = 1d0 - exp(- 2d0 * eta * dem1 * dem1)
                 dem2 = eQP(p) - eQP(i) - eQP(j) + eQP(a)
                 reg2 = 1d0 - exp(- 2d0 * eta * dem2 * dem2)
                 
                 SigC(p) = SigC(p) + num * (reg1/dem1) * (reg2/dem2)
                 Z(p)    = Z(p)    - num * (reg1/dem1) * (reg2/dem2/dem2)
                                  
                 num  = (Kx*0.5d0*ERI(p,i,a,j) - ERI(p,i,j,a)) * &
                 eh_sing_rho(a,i,n) * eh_sing_rho(p,j,n)
                 
                 dem1 = eQP(a) - eQP(i) + eh_sing_Om(n) 
                 dem2 = eQP(p) - eQP(j) + eh_sing_Om(n)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 SigC(p) = SigC(p) + num * (reg1/dem1) * (reg2/dem2)
                 Z(p)    = Z(p)    - num * (reg1/dem1) * (reg2/dem2/dem2)
                 
                 num  = (Kx*0.5d0*ERI(p,a,i,j) - ERI(p,a,j,i))* &
                 eh_sing_rho(a,i,n) * eh_sing_rho(j,p,n)

                 dem1 = eQP(a) - eQP(i) + eh_sing_Om(n) 
                 dem2 = eQP(p) + eQP(a) - eQP(i) - eQP(j)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 SigC(p) = SigC(p) + num * (reg1/dem1) * (reg2/dem2)
                 Z(p)    = Z(p)    - num * (reg1/dem1) * (reg2/dem2/dem2)
                 
              end do ! j
              !3p2h
              do b=nO+1,nOrb-nR

                 num  = (Kx*0.5d0*ERI(p,i,a,b) - ERI(p,i,b,a)) * &
                      eh_sing_rho(i,a,n) * eh_sing_rho(b,p,n)
                 
                 dem1 = eQP(a) - eQP(i) - eh_sing_Om(n)
                 reg1 = 1d0 - exp(- 2d0 * eta * dem1 * dem1)
                 dem2 = eQP(p) - eQP(b) - eh_sing_Om(n)
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))
                 
                 SigC(p) = SigC(p) + num * (reg1/dem1) * (reg2/dem2)
                 Z(p)    = Z(p)    - num * (reg1/dem1) * (reg2/dem2/dem2)
                 
                 num  = - (Kx*0.5d0*ERI(p,i,a,b) - ERI(p,i,b,a)) * &
                      eh_sing_rho(i,a,n) * eh_sing_rho(b,p,n)
                 
                 dem1 = eQP(a) - eQP(i) - eh_sing_Om(n) 
                 reg1 = 1d0 - exp(- 2d0 * eta * dem1 * dem1)
                 dem2 = eQP(p) + eQP(i) - eQP(a) - eQP(b)
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))
                    
                 SigC(p) = SigC(p) + num * (reg1/dem1) * (reg2/dem2)
                 Z(p)    = Z(p)    - num * (reg1/dem1) * (reg2/dem2/dem2)
                 
                 num  = (Kx*0.5d0*ERI(p,a,i,b) - ERI(p,a,b,i)) * &
                 eh_sing_rho(a,i,n) * eh_sing_rho(b,p,n)

                 dem1 = eQP(a) - eQP(i) + eh_sing_Om(n) 
                 dem2 = eQP(p) - eQP(b) - eh_sing_Om(n)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 SigC(p) = SigC(p) + num * (reg1/dem1) * (reg2/dem2)
                 Z(p)    = Z(p)    - num * (reg1/dem1) * (reg2/dem2/dem2)
                 
                 num  = (Kx*0.5d0*ERI(p,i,a,b) - ERI(p,i,b,a)) * &
                 eh_sing_rho(a,i,n) * eh_sing_rho(p,b,n)

                 dem1 = eQP(a) - eQP(i) + eh_sing_Om(n) 
                 dem2 = eQP(p) + eQP(i) - eQP(a) - eQP(b)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 SigC(p) = SigC(p) + num * (reg1/dem1) * (reg2/dem2)
                 Z(p)    = Z(p)    - num * (reg1/dem1) * (reg2/dem2/dem2)

                 
              end do ! b
              
           end do ! n
           
        end do ! a
     end do ! i
     
  end do ! p
  !$OMP END DO
  !$OMP END PARALLEL
  call wall_time(end_t)
  t = end_t - start_t

  write(*,'(1X,A50,1X,F9.3,A8)') 'Wall time for singlet eh self-energy =',t,' seconds'
  write(*,*)

! Self-energy decomposition

  Sig_1eh(:) = SigC(:) - Sig_2d(:) - Sig_2x(:)
  Z_1eh(:)   = Z(:) - Z_2d(:) - Z_2x(:)

!-------------------------------------!
!  triplet eh part of the self-energy !
!-------------------------------------!
   call wall_time(start_t)
   !$OMP PARALLEL DEFAULT(NONE)    &
   !$OMP PRIVATE(p,i,a,j,b,n,num,dem1,dem2,reg1,reg2) &
   !$OMP SHARED(nC,nO,nOrb,nR,nS,eta,ERI,eQP,eh_trip_rho,eh_trip_Om,SigC,Z)
   !$OMP DO
   do p=nC+1,nOrb-nR
    
      do i=nC+1,nO
         do a=nO+1,nOrb-nR
          
            do n=1,nS
               !3h2p
               do j=nC+1,nO

                  num  = 1.5d0*ERI(p,a,i,j) * &
                       eh_trip_rho(i,a,n) * eh_trip_rho(p,j,n)

                  dem1 = eQP(a) - eQP(i) - eh_trip_Om(n)
                  reg1 = 1d0 - exp(- 2d0 * eta * dem1 * dem1)
                  dem2 = eQP(p) - eQP(j) + eh_trip_Om(n)
                  reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))
                     
                  SigC(p) = SigC(p) + num * (reg1/dem1) * (reg2/dem2)
                  Z(p)    = Z(p)    - num * (reg1/dem1) * (reg2/dem2/dem2)
                     
                  num  = - 1.5d0*ERI(p,a,i,j) * &
                       eh_trip_rho(i,a,n) * eh_trip_rho(p,j,n)
                     
                  dem1 = eQP(a) - eQP(i) - eh_trip_Om(n) 
                  reg1 = 1d0 - exp(- 2d0 * eta * dem1 * dem1)
                  dem2 = eQP(p) - eQP(i) - eQP(j) + eQP(a)
                  reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))
                     
                  SigC(p) = SigC(p) + num * (reg1/dem1) * (reg2/dem2)
                  Z(p)    = Z(p)    - num * (reg1/dem1) * (reg2/dem2/dem2)
                                 
                  num  = 1.5d0*ERI(p,i,a,j) * &
                  eh_trip_rho(a,i,n) * eh_trip_rho(p,j,n)
                
                  dem1 = eQP(a) - eQP(i) + eh_trip_Om(n) 
                  dem2 = eQP(p) - eQP(j) + eh_trip_Om(n)
                  reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                  reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                  SigC(p) = SigC(p) + num * (reg1/dem1) * (reg2/dem2)
                  Z(p)    = Z(p)    - num * (reg1/dem1) * (reg2/dem2/dem2)
                
                  num  = 1.5d0*ERI(p,a,i,j) * &
                  eh_trip_rho(a,i,n) * eh_trip_rho(j,p,n)

                  dem1 = eQP(a) - eQP(i) + eh_trip_Om(n) 
                  dem2 = eQP(p) + eQP(a) - eQP(i) - eQP(j)
                  reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                  reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                  SigC(p) = SigC(p) + num * (reg1/dem1) * (reg2/dem2)
                  Z(p)    = Z(p)    - num * (reg1/dem1) * (reg2/dem2/dem2)
                
               end do ! j
               !3p2h
               do b=nO+1,nOrb-nR

                  num  = 1.5d0*ERI(p,i,a,b) * &
                       eh_trip_rho(i,a,n) * eh_trip_rho(b,p,n)
                     
                  dem1 = eQP(a) - eQP(i) - eh_trip_Om(n)
                  reg1 = 1d0 - exp(- 2d0 * eta * dem1 * dem1)
                  dem2 = eQP(p) - eQP(b) - eh_trip_Om(n)
                  reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))
                     
                  SigC(p) = SigC(p) + num * (reg1/dem1) * (reg2/dem2)
                  Z(p)    = Z(p)    - num * (reg1/dem1) * (reg2/dem2/dem2)
                     
                  num  = - 1.5d0*ERI(p,i,a,b) * &
                       eh_trip_rho(i,a,n) * eh_trip_rho(b,p,n)
                     
                  dem1 = eQP(a) - eQP(i) - eh_trip_Om(n) 
                  reg1 = 1d0 - exp(- 2d0 * eta * dem1 * dem1)
                  dem2 = eQP(p) + eQP(i) - eQP(a) - eQP(b)
                  reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))
                    
                  SigC(p) = SigC(p) + num * (reg1/dem1) * (reg2/dem2)
                  Z(p)    = Z(p)    - num * (reg1/dem1) * (reg2/dem2/dem2)
                
                  num  = 1.5d0*ERI(p,a,i,b) * &
                  eh_trip_rho(a,i,n) * eh_trip_rho(b,p,n)

                  dem1 = eQP(a) - eQP(i) + eh_trip_Om(n) 
                  dem2 = eQP(p) - eQP(b) - eh_trip_Om(n)
                  reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                  reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                  SigC(p) = SigC(p) + num * (reg1/dem1) * (reg2/dem2)
                  Z(p)    = Z(p)    - num * (reg1/dem1) * (reg2/dem2/dem2)
                
                  num  = 1.5d0*ERI(p,i,a,b) * &
                  eh_trip_rho(a,i,n) * eh_trip_rho(p,b,n)

                  dem1 = eQP(a) - eQP(i) + eh_trip_Om(n) 
                  dem2 = eQP(p) + eQP(i) - eQP(a) - eQP(b)
                  reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                  reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                  SigC(p) = SigC(p) + num * (reg1/dem1) * (reg2/dem2)
                  Z(p)    = Z(p)    - num * (reg1/dem1) * (reg2/dem2/dem2)

                
               end do ! b
             
            end do ! n
          
         end do ! a
      end do ! i
    
   end do ! p
   !$OMP END DO
   !$OMP END PARALLEL
   call wall_time(end_t)
   t = end_t - start_t

   write(*,'(1X,A50,1X,F9.3,A8)') 'Wall time for triplet eh self-energy =',t,' seconds'
   write(*,*)

! Self-energy decomposition

  Sig_3eh(:) = SigC(:) - Sig_2d(:) - Sig_2x(:) - Sig_1eh(:)
  Z_3eh(:)   = Z(:) - Z_2d(:) - Z_2x(:) - Z_1eh(:)
  
!-------------------------------------!
!  singlet pp part of the self-energy !
!-------------------------------------!
  call wall_time(start_t)
  !$OMP PARALLEL DEFAULT(NONE)    &
  !$OMP PRIVATE(p,i,j,k,c,n,num,dem1,dem2,reg1,reg2) &
  !$OMP SHARED(nC,nO,nOrb,nR,nOOs,nVVs,eta,ERI,eQP,ee_sing_rho,ee_sing_Om,hh_sing_rho,hh_sing_Om,SigC,Z)
  !$OMP DO
  do p=nC+1,nOrb-nR
     
     do i=nC+1,nO
        do j=nC+1,nO
           do n=1,nVVs
              ! 4h1p
              do k=nC+1,nO
                 num  = - 0.5d0 * ERI(p,k,i,j) * ee_sing_rho(i,j,n) * ee_sing_rho(p,k,n)
                 dem1 = ee_sing_Om(n) - eQP(i) - eQP(j)
                 dem2 = eQP(p) + eQP(k) - ee_sing_Om(n)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 SigC(p) = SigC(p) + num * (reg1/dem1) * (reg2/dem2)
                 Z(p)    = Z(p)    - num * (reg1/dem1) * (reg2/dem2/dem2)

              end do ! k
              ! 3h2p
              do c=nO+1,nOrb-nR

                 num  = - 0.5d0*ERI(p,c,i,j) * ee_sing_rho(i,j,n) * ee_sing_rho(p,c,n)
                 dem1 = ee_sing_Om(n) - eQP(i) - eQP(j)
                 dem2 = eQP(p) + eQP(c) - eQP(i) - eQP(j)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 SigC(p) = SigC(p) + num * (reg1/dem1) * (reg2/dem2)
                 Z(p)    = Z(p)    - num * (reg1/dem1) * (reg2/dem2/dem2)
                 
              end do ! a
           end do ! n
           do n=1,nOOs
              ! 3h2p
              do a=nO+1,nOrb-nR

                 num  = 0.5d0*ERI(p,a,i,j) * hh_sing_rho(i,j,n) * hh_sing_rho(p,a,n)
                 
                 dem1 = eQP(i) + eQP(j) - hh_sing_Om(n)
                 reg1 = 1d0 - exp(- 2d0*eta * dem1 * dem1)
                 dem2 = eQP(p) + eQP(a) - hh_sing_Om(n)
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))
                    
                 SigC(p) = SigC(p) + num * (reg1/dem1) * (reg2/dem2)
                 Z(p)    = Z(p)    - num * (reg1/dem1) * (reg2/dem2/dem2)
                    
                 num  = - 0.5d0*ERI(p,a,i,j) * hh_sing_rho(i,j,n) * hh_sing_rho(p,a,n)

                 dem1 = eQP(i) + eQP(j) - hh_sing_Om(n)
                 reg1 = 1d0 - exp(- 2d0*eta * dem1 * dem1)
                 dem2 = eQP(p) + eQP(a) - eQP(i) - eQP(j)
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))
                    
                 SigC(p) = SigC(p) + num * (reg1/dem1) * (reg2/dem2)
                 Z(p)    = Z(p)    - num * (reg1/dem1) * (reg2/dem2/dem2)
                 
              end do ! a
           end do ! n
        end do ! j
     end do ! i
     
  end do ! p
  !$OMP END DO
  !$OMP END PARALLEL
  !$OMP PARALLEL DEFAULT(NONE)    &
  !$OMP PRIVATE(p,k,a,b,c,n,num,dem1,dem2,reg1,reg2) &
  !$OMP SHARED(nC,nO,nOrb,nR,nOOs,nVVs,eta,ERI,eQP,ee_sing_rho,ee_sing_Om,hh_sing_rho,hh_sing_Om,SigC,Z)
  !$OMP DO
  do p=nC+1,nOrb-nR
     do a=nO+1,nOrb-nR
        do b=nO+1,nOrb-nR
           do n=1,nOOs
              ! 4p1h
              do c=nO+1,nOrb-nR

                 num  = 0.5d0*ERI(p,c,a,b) * hh_sing_rho(a,b,n) * hh_sing_rho(p,c,n)
                 dem1 = hh_sing_Om(n) - eQP(a) - eQP(b)
                 dem2 = eQP(p) + eQP(c) - hh_sing_Om(n)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 SigC(p) = SigC(p) + num * (reg1/dem1) * (reg2/dem2)
                 Z(p)    = Z(p)    - num * (reg1/dem1) * (reg2/dem2/dem2)

              end do ! c
              ! 3p2h
              do k=nC+1,nO

                 num  = 0.5d0*ERI(p,k,a,b) * hh_sing_rho(a,b,n) * hh_sing_rho(p,k,n)
                 dem1 = hh_sing_Om(n) - eQP(a) - eQP(b)
                 dem2 = eQP(p) + eQP(k) - eQP(a) - eQP(b)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 SigC(p) = SigC(p) + num * (reg1/dem1) * (reg2/dem2)
                 Z(p)    = Z(p)    - num * (reg1/dem1) * (reg2/dem2/dem2)

              end do ! k
           end do ! n
           do n=1,nVVs
              ! 3p2h
              do i=nC+1,nO

                 num  = 0.5d0*ERI(p,i,a,b) * ee_sing_rho(a,b,n) * ee_sing_rho(p,i,n)

                 dem1 = eQP(a) + eQP(b) - ee_sing_Om(n)
                 reg1 = 1d0 - exp(- 2d0 * eta * dem1 * dem1)
                 dem2 = eQP(p) + eQP(i) - eQP(a) - eQP(b)
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))
                    
                 SigC(p) = SigC(p) + num * (reg1/dem1) * (reg2/dem2)
                 Z(p)    = Z(p)    - num * (reg1/dem1) * (reg2/dem2/dem2)

                 num  = - 0.5d0*ERI(p,i,a,b) * ee_sing_rho(a,b,n) * ee_sing_rho(p,i,n)

                 dem1 = eQP(a) + eQP(b) - ee_sing_Om(n)
                 reg1 = 1d0 - exp(- 2d0 * eta * dem1 * dem1)
                 dem2 = eQP(p) + eQP(i) - ee_sing_Om(n)
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))
                    
                 SigC(p) = SigC(p) + num * (reg1/dem1) * (reg2/dem2)
                 Z(p)    = Z(p)    - num * (reg1/dem1) * (reg2/dem2/dem2)
                 
              end do ! i
           end do ! n
        end do ! b
     end do ! a
     
  end do ! p
  !$OMP END DO
  !$OMP END PARALLEL
  call wall_time(end_t)
  t = end_t - start_t

  write(*,'(1X,A50,1X,F9.3,A8)') 'Wall time for singlet pp self-energy =',t,' seconds'
  write(*,*)

! Self-energy decomposition

  Sig_1pp(:) = SigC(:) - Sig_2d(:) - Sig_2x(:) - Sig_1eh(:) - Sig_3eh(:)
  Z_1pp(:)   = Z(:) - Z_2d(:) - Z_2x(:) - Z_1eh(:) - Z_3eh(:)

!-------------------------------------!
!  triplet pp part of the self-energy !
!-------------------------------------!
  call wall_time(start_t)
  !$OMP PARALLEL DEFAULT(NONE)    &
  !$OMP PRIVATE(p,i,j,k,c,n,num,dem1,dem2,reg1,reg2) &
  !$OMP SHARED(nC,nO,nOrb,nR,nOOt,nVVt,eta,ERI,eQP,ee_trip_rho,ee_trip_Om,hh_trip_rho,hh_trip_Om,SigC,Z)
  !$OMP DO
  do p=nC+1,nOrb-nR
    
     do i=nC+1,nO
        do j=nC+1,nO
           do n=1,nVVt
              ! 4h1p
              do k=nC+1,nO
                 num  = - 1.5d0 * ERI(p,k,i,j) * ee_trip_rho(i,j,n) * ee_trip_rho(p,k,n)
                 dem1 = ee_trip_Om(n) - eQP(i) - eQP(j)
                 dem2 = eQP(p) + eQP(k) - ee_trip_Om(n)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 SigC(p) = SigC(p) + num * (reg1/dem1) * (reg2/dem2)
                 Z(p)    = Z(p)    - num * (reg1/dem1) * (reg2/dem2/dem2)

              end do ! k
              ! 3h2p
              do c=nO+1,nOrb-nR

                 num  = - 1.5d0 * ERI(p,c,i,j) * ee_trip_rho(i,j,n) * ee_trip_rho(p,c,n)
                 dem1 = ee_trip_Om(n) - eQP(i) - eQP(j)
                 dem2 = eQP(p) + eQP(c) - eQP(i) - eQP(j)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 SigC(p) = SigC(p) + num * (reg1/dem1) * (reg2/dem2)
                 Z(p)    = Z(p)    - num * (reg1/dem1) * (reg2/dem2/dem2)
                 
              end do ! a
           end do ! n
           do n=1,nOOt
              ! 3h2p
              do a=nO+1,nOrb-nR
                 
                 num  = 1.5d0 * ERI(p,a,i,j) * hh_trip_rho(i,j,n) * hh_trip_rho(p,a,n)
                 
                 dem1 = eQP(i) + eQP(j) - hh_trip_Om(n)
                 reg1 = 1d0 - exp(- 2d0 * eta * dem1 * dem1)
                 dem2 = eQP(p) + eQP(a) - hh_trip_Om(n)
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))
                 
                 SigC(p) = SigC(p) + num * (reg1/dem1) * (reg2/dem2)
                 Z(p)    = Z(p)    - num * (reg1/dem1) * (reg2/dem2/dem2)
                 
                 num  = - 1.5d0 * ERI(p,a,i,j) * hh_trip_rho(i,j,n) * hh_trip_rho(p,a,n)
                 
                 dem1 = eQP(i) + eQP(j) - hh_trip_Om(n)
                 reg1 = 1d0 - exp(- 2d0 * eta * dem1 * dem1)
                 dem2 = eQP(p) + eQP(a) - eQP(i) - eQP(j)
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))
                 
                 SigC(p) = SigC(p) + num * (reg1/dem1) * (reg2/dem2)
                 Z(p)    = Z(p)    - num * (reg1/dem1) * (reg2/dem2/dem2)
                 
              end do ! a
           end do ! n
        end do ! j
     end do ! i
     
  end do ! p
  !$OMP END DO
  !$OMP END PARALLEL
  !$OMP PARALLEL DEFAULT(NONE)    &
  !$OMP PRIVATE(p,k,a,b,c,n,num,dem1,dem2,reg1,reg2) &
  !$OMP SHARED(nC,nO,nOrb,nR,nOOt,nVVt,eta,ERI,eQP,ee_trip_rho,ee_trip_Om,hh_trip_rho,hh_trip_Om,SigC,Z)
  !$OMP DO
  do p=nC+1,nOrb-nR
     do a=nO+1,nOrb-nR
        do b=nO+1,nOrb-nR
           do n=1,nOOt
              ! 4p1h
              do c=nO+1,nOrb-nR

                 num  = 1.5d0 * ERI(p,c,a,b) * hh_trip_rho(a,b,n) * hh_trip_rho(p,c,n)
                 dem1 = hh_trip_Om(n) - eQP(a) - eQP(b)
                 dem2 = eQP(p) + eQP(c) - hh_trip_Om(n)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 SigC(p) = SigC(p) + num * (reg1/dem1) * (reg2/dem2)
                 Z(p)    = Z(p)    - num * (reg1/dem1) * (reg2/dem2/dem2)

              end do ! c
              ! 3p2h
              do k=nC+1,nO

                 num  = 1.5d0 * ERI(p,k,a,b) * hh_trip_rho(a,b,n) * hh_trip_rho(p,k,n)
                 dem1 = hh_trip_Om(n) - eQP(a) - eQP(b)
                 dem2 = eQP(p) + eQP(k) - eQP(a) - eQP(b)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 SigC(p) = SigC(p) + num * (reg1/dem1) * (reg2/dem2)
                 Z(p)    = Z(p)    - num * (reg1/dem1) * (reg2/dem2/dem2)

              end do ! k
           end do ! n
           do n=1,nVVt
              ! 3p2h
              do i=nC+1,nO

                 num  = 1.5d0 * ERI(p,i,a,b) * ee_trip_rho(a,b,n) * ee_trip_rho(p,i,n)

                 dem1 = eQP(a) + eQP(b) - ee_trip_Om(n)
                 reg1 = 1d0 - exp(- 2d0 * eta * dem1 * dem1)
                 dem2 = eQP(p) + eQP(i) - eQP(a) - eQP(b)
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 SigC(p) = SigC(p) + num * (reg1/dem1) * (reg2/dem2)
                 Z(p)    = Z(p)    - num * (reg1/dem1) * (reg2/dem2/dem2)
                 
                 num  = - 1.5d0 * ERI(p,i,a,b) * ee_trip_rho(a,b,n) * ee_trip_rho(p,i,n)

                 dem1 = eQP(a) + eQP(b) - ee_trip_Om(n)
                 reg1 = 1d0 - exp(- 2d0 * eta * dem1 * dem1)
                 dem2 = eQP(p) + eQP(i) - ee_trip_Om(n)
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))
                 
                 SigC(p) = SigC(p) + num * (reg1/dem1) * (reg2/dem2)
                 Z(p)    = Z(p)    - num * (reg1/dem1) * (reg2/dem2/dem2)
                 
              end do ! i
           end do ! n
        end do ! b
     end do ! a
     
  end do ! p
  !$OMP END DO
  !$OMP END PARALLEL 
  call wall_time(end_t)
  t = end_t - start_t

  write(*,'(1X,A50,1X,F9.3,A8)') 'Wall time for triplet pp self-energy =',t,' seconds'
  write(*,*)
 
! Self-energy decomposition

  Sig_3pp(:) = SigC(:) - Sig_2d(:) - Sig_2x(:) - Sig_1eh(:) - Sig_3eh(:) - Sig_1pp(:)
  Z_3pp(:)   = Z(:) - Z_2d(:) - Z_2x(:) - Z_1eh(:) - Z_3eh(:) - Z_1pp(:)
 
!----------------------!
! Channel restrictions !
!----------------------!

  SigC(:) = 0d0

  if(do_2d_channel) SigC(:) = SigC(:) + Sig_2d(:) 
  if(do_2x_channel) SigC(:) = SigC(:) + Sig_2x(:)
  if(do_1eh_channel) SigC(:) = SigC(:) + Sig_1eh(:)
  if(do_3eh_channel) SigC(:) = SigC(:) + Sig_3eh(:)
  if(do_1pp_channel) SigC(:) = SigC(:) + Sig_1pp(:)
  if(do_3pp_channel) SigC(:) = SigC(:) + Sig_3pp(:)

!-----------------------------!
!   Renormalization factor    !
!-----------------------------!

  Z(:) = 0d0

  if(do_2d_channel) Z(:) = Z(:) + Z_2d(:)
  if(do_2x_channel) Z(:) = Z(:) + Z_2x(:)
  if(do_1eh_channel) Z(:) = Z(:) + Z_1eh(:)
  if(do_3eh_channel) Z(:) = Z(:) + Z_3eh(:)
  if(do_1pp_channel) Z(:) = Z(:) + Z_1pp(:)
  if(do_3pp_channel) Z(:) = Z(:) + Z_3pp(:)

  Z(:) = 1d0/(1d0 - Z(:))
  
!-------------------------------------!
! Galitskii-Migdal correlation energy !
!-------------------------------------!

  EcGM = 0d0
  
!---------------------------------!
! Print self-energy decomposition !
!---------------------------------!

  if(print_self_energy) &
    call dump_RParquet_self_energy(nOrb,nC,nO,nV,nR,Sig_2d,Sig_2x,Sig_1eh,Sig_3eh,Sig_1pp,Sig_3pp,SigC)

end subroutine 
