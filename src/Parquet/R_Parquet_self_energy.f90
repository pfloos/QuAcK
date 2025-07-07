subroutine R_Parquet_self_energy(eta,nOrb,nC,nO,nV,nR,nS,nOOs,nVVs,nOOt,nVVt,eQP,ERI, &
                                 eh_sing_rho,eh_sing_Om,eh_trip_rho,eh_trip_Om, &
                                 ee_sing_rho,ee_sing_Om,ee_trip_rho,ee_trip_Om, &
                                 hh_sing_rho,hh_sing_Om,hh_trip_rho,hh_trip_Om, &
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
  double precision,intent(in)   :: eh_sing_Om(nS+nS)
  double precision,intent(in)   :: eh_trip_rho(nOrb,nOrb,nS)
  double precision,intent(in)   :: eh_trip_Om(nS+nS)
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

  logical                       :: print_self_energy
  double precision,allocatable  :: Sig2d(:)
  double precision,allocatable  :: Sig2x(:)
  double precision,allocatable  :: Sig1eh(:)
  double precision,allocatable  :: Sig3eh(:)
  double precision,allocatable  :: Sig1pp(:)
  double precision,allocatable  :: Sig3pp(:)
  
! Output variables
  double precision,intent(out)  :: EcGM
  double precision,intent(out)  :: SigC(nOrb)
  double precision,intent(out)  :: Z(nOrb)

! Initialize 

  SigC(:) = 0d0
  Z(:)    = 0d0
  EcGM    = 0d0

! Memory allocation for self-energy decomposition

  allocate(Sig2d(nOrb))
  allocate(Sig2x(nOrb))
  allocate(Sig1eh(nOrb))
  allocate(Sig3eh(nOrb))
  allocate(Sig1pp(nOrb))
  allocate(Sig3pp(nOrb))
  
  Sig2x(:) = 0d0
  Sig2d(:) = 0d0

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

              Sig2d(p) = Sig2d(p) + num*reg/eps
              Z(p)     = Z(p)     - num*reg/eps**2

              num = - ERI(p,a,j,i)*ERI(j,i,a,p)

              Sig2x(p) = Sig2x(p) + num*reg/eps
              Z(p)     = Z(p)     - num*reg/eps**2

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

              Sig2d(p) = Sig2d(p) + num*reg/eps
              Z(p)     = Z(p)    - num*reg/eps**2

              num = - ERI(p,i,b,a)*ERI(b,a,i,p)

              Sig2x(p) = Sig2x(p) + num*reg/eps
              Z(p)    = Z(p)    - num*reg/eps**2

           end do
        end do
     end do
  end do
  call wall_time(end_t)
  t = end_t - start_t

  write(*,'(1X,A50,1X,F9.3,A8)') 'Wall time for 2nd-order self-energy =',t,' seconds'
  write(*,*)

! Self-energy decomposition

  SigC(:) = Sig2d(:) + Sig2x(:)

!-------------------------------------!
!  singlet eh part of the self-energy !
!-------------------------------------!
  call wall_time(start_t)
  ! !$OMP PARALLEL DEFAULT(NONE)    &
  ! !$OMP PRIVATE(p,i,a,j,b,n,num,dem1,dem2,reg1,reg2) &
  ! !$OMP SHARED(nC,nO,nOrb,nR,nS,eta,ERI,eQP,eh_sing_rho,eh_sing_Om,SigC,Z)
  ! !$OMP DO COLLAPSE(2)
  do p=nC+1,nOrb-nR
     
     do i=nC+1,nO
        do a=nO+1,nOrb-nR
           
           do n=1,nS
              !3h2p
              do j=nC+1,nO
                 num  = (0.5d0*ERI(p,a,i,j) - ERI(p,a,j,i))* &
                 eh_sing_rho(i,a,n) * eh_sing_rho(j,p,nS+n)

                 dem1 = eQP(p) - eQP(j) + eh_sing_Om(n)
                 dem2 = eQP(p) - eQP(j) + eQP(a) - eQP(i)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 SigC(p) = SigC(p) + num * (reg1/dem1) * (reg2/dem2)
                 Z(p)    = Z(p)    - num * (reg1/dem1) * (reg2/dem2/dem2) &
                                   - num * (reg1/dem1/dem1) * (reg2/dem2) 
                                  
                 num  = (0.5d0*ERI(p,i,a,j) - ERI(p,i,j,a)) * &
                 eh_sing_rho(a,i,n) * eh_sing_rho(j,p,nS+n)
                 
                 dem1 = eQP(a) - eQP(i) + eh_sing_Om(n) 
                 dem2 = eQP(p) - eQP(j) + eh_sing_Om(n)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 SigC(p) = SigC(p) + num * (reg1/dem1) * (reg2/dem2)
                 Z(p)    = Z(p)    - num * (reg1/dem1) * (reg2/dem2/dem2)
                 
                 num  = (0.5d0*ERI(p,a,i,j) - ERI(p,a,j,i))* &
                 eh_sing_rho(i,a,nS+n) * eh_sing_rho(j,p,n)

                 dem1 = eQP(a) - eQP(i) + eh_sing_Om(n) 
                 dem2 = eQP(p) + eQP(a) - eQP(i) - eQP(j)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 SigC(p) = SigC(p) + num * (reg1/dem1) * (reg2/dem2)
                 Z(p)    = Z(p)    - num * (reg1/dem1) * (reg2/dem2/dem2)
                 
              end do ! j
              !3p2h
              do b=nO+1,nOrb-nR
                 num  = -(0.5d0*ERI(p,i,a,b) - ERI(p,i,b,a)) * &
                 eh_sing_rho(a,i,nS+n) * eh_sing_rho(b,p,n)

                 dem1 = eQP(p) - eQP(b) - eh_sing_Om(n) 
                 dem2 = eQP(p) - eQP(b) - eQP(a) + eQP(i)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 SigC(p) = SigC(p) + num * (reg1/dem1) * (reg2/dem2)
                 Z(p)    = Z(p)    - num * (reg1/dem1) * (reg2/dem2/dem2) &
                                   - num * (reg1/dem1/dem1) * (reg2/dem2)
                 
                 num  = (0.5d0*ERI(p,a,i,b) - ERI(p,a,b,i)) * &
                 eh_sing_rho(i,a,nS+n) * eh_sing_rho(b,p,n)

                 dem1 = eQP(a) - eQP(i) + eh_sing_Om(n) 
                 dem2 = eQP(p) - eQP(b) - eh_sing_Om(n)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 SigC(p) = SigC(p) + num * (reg1/dem1) * (reg2/dem2)
                 Z(p)    = Z(p)    - num * (reg1/dem1) * (reg2/dem2/dem2)
                 
                 num  = (0.5d0*ERI(p,i,a,b) - ERI(p,i,b,a)) * &
                 eh_sing_rho(a,i,n) * eh_sing_rho(b,p,nS+n)

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
  ! !$OMP END DO
  ! !$OMP END PARALLEL
  call wall_time(end_t)
  t = end_t - start_t

  write(*,'(1X,A50,1X,F9.3,A8)') 'Wall time for singlet eh self-energy =',t,' seconds'
  write(*,*)

! Self-energy decomposition

  Sig1eh(:) = SigC(:) - Sig2d(:) - Sig2x(:)

!-------------------------------------!
!  triplet eh part of the self-energy !
!-------------------------------------!
   call wall_time(start_t)
   ! !$OMP PARALLEL DEFAULT(NONE)    &
   ! !$OMP PRIVATE(p,i,a,j,b,n,num,dem1,dem2,reg1,reg2) &
   ! !$OMP SHARED(nC,nO,nOrb,nR,nS,eta,ERI,eQP,eh_trip_rho,eh_trip_Om,SigC,Z)
   ! !$OMP DO COLLAPSE(2)
   do p=nC+1,nOrb-nR
    
      do i=nC+1,nO
         do a=nO+1,nOrb-nR
          
            do n=1,nS
               !3h2p
               do j=nC+1,nO
                  num  = (1.5d0*ERI(p,a,i,j) - 0d0*ERI(p,a,j,i))* &
                  eh_trip_rho(i,a,n) * eh_trip_rho(j,p,nS+n)

                  dem1 = eQP(p) - eQP(j) + eh_trip_Om(n)
                  dem2 = eQP(p) - eQP(j) + eQP(a) - eQP(i)
                  reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                  reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                  SigC(p) = SigC(p) + num * (reg1/dem1) * (reg2/dem2)
                  Z(p)    = Z(p)    - num * (reg1/dem1) * (reg2/dem2/dem2) &
                                    - num * (reg1/dem1/dem1) * (reg2/dem2) 
                                 
                  num  = (1.5d0*ERI(p,i,a,j) - 0d0*ERI(p,i,j,a)) * &
                  eh_trip_rho(a,i,n) * eh_trip_rho(j,p,nS+n)
                
                  dem1 = eQP(a) - eQP(i) + eh_trip_Om(n) 
                  dem2 = eQP(p) - eQP(j) + eh_trip_Om(n)
                  reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                  reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                  SigC(p) = SigC(p) + num * (reg1/dem1) * (reg2/dem2)
                  Z(p)    = Z(p)    - num * (reg1/dem1) * (reg2/dem2/dem2)
                
                  num  = (1.5d0*ERI(p,a,i,j) - 0d0*ERI(p,a,j,i))* &
                  eh_trip_rho(i,a,nS+n) * eh_trip_rho(j,p,n)

                  dem1 = eQP(a) - eQP(i) + eh_trip_Om(n) 
                  dem2 = eQP(p) + eQP(a) - eQP(i) - eQP(j)
                  reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                  reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                  SigC(p) = SigC(p) + num * (reg1/dem1) * (reg2/dem2)
                  Z(p)    = Z(p)    - num * (reg1/dem1) * (reg2/dem2/dem2)
                
               end do ! j
               !3p2h
               do b=nO+1,nOrb-nR
                  num  = -(1.5d0*ERI(p,i,a,b) - 0d0*ERI(p,i,b,a)) * &
                  eh_trip_rho(a,i,nS+n) * eh_trip_rho(b,p,n)

                  dem1 = eQP(p) - eQP(b) - eh_trip_Om(n) 
                  dem2 = eQP(p) - eQP(b) - eQP(a) + eQP(i)
                  reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                  reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                  SigC(p) = SigC(p) + num * (reg1/dem1) * (reg2/dem2)
                  Z(p)    = Z(p)    - num * (reg1/dem1) * (reg2/dem2/dem2) &
                                    - num * (reg1/dem1/dem1) * (reg2/dem2)
                
                  num  = (1.5d0*ERI(p,a,i,b) - 0d0*ERI(p,a,b,i)) * &
                  eh_trip_rho(i,a,nS+n) * eh_trip_rho(b,p,n)

                  dem1 = eQP(a) - eQP(i) + eh_trip_Om(n) 
                  dem2 = eQP(p) - eQP(b) - eh_trip_Om(n)
                  reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                  reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                  SigC(p) = SigC(p) + num * (reg1/dem1) * (reg2/dem2)
                  Z(p)    = Z(p)    - num * (reg1/dem1) * (reg2/dem2/dem2)
                
                  num  = (1.5d0*ERI(p,i,a,b) - 0d0*ERI(p,i,b,a)) * &
                  eh_trip_rho(a,i,n) * eh_trip_rho(b,p,nS+n)

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
   ! !$OMP END DO
   ! !$OMP END PARALLEL
   call wall_time(end_t)
   t = end_t - start_t

   write(*,'(1X,A50,1X,F9.3,A8)') 'Wall time for triplet eh self-energy =',t,' seconds'
   write(*,*)

! Self-energy decomposition

  Sig3eh(:) = SigC(:) - Sig2d(:) - Sig2x(:) - Sig1eh(:)
  
!-------------------------------------!
!  singlet pp part of the self-energy !
!-------------------------------------!
  call wall_time(start_t)
  ! !$OMP PARALLEL DEFAULT(NONE)    &
  ! !$OMP PRIVATE(p,i,j,k,c,n,num,dem1,dem2,reg1,reg2) &
  ! !$OMP SHARED(nC,nO,nOrb,nR,nOOs,nVVs,eta,ERI,eQP,ee_sing_rho,ee_sing_Om,hh_sing_rho,hh_sing_Om,SigC,Z)
  ! !$OMP DO COLLAPSE(2)
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
              do c=nO+1,nOrb-nR

                 num  = - 0.5d0*ERI(p,c,i,j) * hh_sing_rho(i,j,n) * hh_sing_rho(p,c,n)
                 dem1 = eQP(p) + eQP(c) - hh_sing_Om(n)
                 dem2 = eQP(p) + eQP(c) - eQP(i) - eQP(j)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 SigC(p) = SigC(p) + num * (reg1/dem1) * (reg2/dem2)
                 Z(p)    = Z(p)    - num * (reg1/dem1) * (reg2/dem2/dem2) &
                                   - num * (reg1/dem1/dem1) * (reg2/dem2)
                 
              end do ! c
           end do ! n
        end do ! j
     end do ! i
     
  end do ! p
  ! !$OMP END DO
  ! !$OMP END PARALLEL
  ! !$OMP PARALLEL DEFAULT(NONE)    &
  ! !$OMP PRIVATE(p,k,a,b,c,n,num,dem1,dem2,reg1,reg2) &
  ! !$OMP SHARED(nC,nO,nOrb,nR,nOOs,nVVs,eta,ERI,eQP,ee_sing_rho,ee_sing_Om,hh_sing_rho,hh_sing_Om,SigC,Z)
  ! !$OMP DO COLLAPSE(2)
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
              do k=nC+1,nO

                 num  = 0.5d0*ERI(p,k,a,b) * ee_sing_rho(a,b,n) * ee_sing_rho(p,k,n)
                 dem1 = eQP(p) + eQP(k) - eQP(a) - eQP(b)
                 dem2 = eQP(p) + eQP(k) - ee_sing_Om(n)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 SigC(p) = SigC(p) + num * (reg1/dem1) * (reg2/dem2)
                 Z(p)    = Z(p)    - num * (reg1/dem1) * (reg2/dem2/dem2) &
                                   - num * (reg1/dem1/dem1) * (reg2/dem2)
                 
              end do ! c
           end do ! n
        end do ! b
     end do ! a
     
  end do ! p
  ! !$OMP END DO
  ! !$OMP END PARALLEL
  call wall_time(end_t)
  t = end_t - start_t

  write(*,'(1X,A50,1X,F9.3,A8)') 'Wall time for singlet pp self-energy =',t,' seconds'
  write(*,*)

! Self-energy decomposition

  Sig1pp(:) = SigC(:) - Sig2d(:) - Sig2x(:) - Sig1eh(:) - Sig3eh(:)

!-------------------------------------!
!  triplet pp part of the self-energy !
!-------------------------------------!
  call wall_time(start_t)
  ! !$OMP PARALLEL DEFAULT(NONE)    &
  ! !$OMP PRIVATE(p,i,j,k,c,n,num,dem1,dem2,reg1,reg2) &
  ! !$OMP SHARED(nC,nO,nOrb,nR,nOOt,nVVt,eta,ERI,eQP,ee_trip_rho,ee_trip_Om,hh_trip_rho,hh_trip_Om,SigC,Z)
  ! !$OMP DO COLLAPSE(2)
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
              do c=nO+1,nOrb-nR

                 num  = - 1.5d0 * ERI(p,c,i,j) * hh_trip_rho(i,j,n) * hh_trip_rho(p,c,n)
                 dem1 = eQP(p) + eQP(c) - hh_trip_Om(n)
                 dem2 = eQP(p) + eQP(c) - eQP(i) - eQP(j)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 SigC(p) = SigC(p) + num * (reg1/dem1) * (reg2/dem2)
                 Z(p)    = Z(p)    - num * (reg1/dem1) * (reg2/dem2/dem2) &
                                   - num * (reg1/dem1/dem1) * (reg2/dem2)
                 
              end do ! c
           end do ! n
        end do ! j
     end do ! i
     
  end do ! p
  ! !$OMP END DO
  ! !$OMP END PARALLEL
  ! !$OMP PARALLEL DEFAULT(NONE)    &
  ! !$OMP PRIVATE(p,k,a,b,c,n,num,dem1,dem2,reg1,reg2) &
  ! !$OMP SHARED(nC,nO,nOrb,nR,nOOt,nVVt,eta,ERI,eQP,ee_trip_rho,ee_trip_Om,hh_trip_rho,hh_trip_Om,SigC,Z)
  ! !$OMP DO COLLAPSE(2)
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
              do k=nC+1,nO

                 num  = 1.5d0 * ERI(p,k,a,b) * ee_trip_rho(a,b,n) * ee_trip_rho(p,k,n)
                 dem1 = eQP(p) + eQP(k) - eQP(a) - eQP(b)
                 dem2 = eQP(p) + eQP(k) - ee_trip_Om(n)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 SigC(p) = SigC(p) + num * (reg1/dem1) * (reg2/dem2)
                 Z(p)    = Z(p)    - num * (reg1/dem1) * (reg2/dem2/dem2) &
                                   - num * (reg1/dem1/dem1) * (reg2/dem2)
                 
              end do ! c
           end do ! n
        end do ! b
     end do ! a
     
  end do ! p
  ! !$OMP END DO
  ! !$OMP END PARALLEL 
  call wall_time(end_t)
  t = end_t - start_t

  write(*,'(1X,A50,1X,F9.3,A8)') 'Wall time for triplet pp self-energy =',t,' seconds'
  write(*,*)
 
! Self-energy decomposition

  Sig3pp(:) = SigC(:) - Sig2d(:) - Sig2x(:) - Sig1eh(:) - Sig3eh(:) - Sig1pp(:)
 
!-----------------------------!
!   Renormalization factor    !
!-----------------------------!

  Z(:) = 1d0/(1d0 - Z(:))
  
!-------------------------------------!
! Galitskii-Migdal correlation energy !
!-------------------------------------!

  EcGM = 0d0
  
!---------------------------------!
! Print self-energy decomposition !
!---------------------------------!

  print_self_energy = .true.
  call dump_RParquet_self_energy(nOrb,nC,nO,nV,nR,Sig2d,Sig2x,Sig1eh,Sig3eh,Sig1pp,Sig3pp,SigC)

end subroutine 
