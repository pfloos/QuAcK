subroutine G_Parquet_self_energy(eta,nOrb,nC,nO,nV,nR,nS,nOO,nVV,eQP,ERI,&
                                 eh_rho,eh_Om,ee_rho,ee_Om,hh_rho,hh_Om,EcGM,SigC,Z)

! Compute correlation part of the self-energy coming from irreducible vertices contribution

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC, nO, nV, nR
  integer,intent(in)            :: nS, nOO, nVV
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
  integer                       :: p,n
  double precision              :: eps,dem1,dem2,reg,reg1,reg2
  double precision              :: num
  double precision              :: start_t,end_t,t

! Output variables

  double precision,intent(out)  :: SigC(nOrb)
  double precision,intent(out)  :: Z(nOrb)
  double precision,intent(out)  :: EcGM

! Initialize 

  SigC(:) = 0d0
  Z(:)    = 0d0
  EcGM    = 0d0
  
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
              num = 0.5d0*(ERI(p,a,j,i) - ERI(p,a,i,j))**2

              SigC(p) = SigC(p) + num*reg/eps
              Z(p)    = Z(p)    - num*reg/eps**2

           end do
        end do
     end do
     ! 2p1h sum
     do i=nC+1,nO
        do a=nO+1,nOrb-nR
           do b=nO+1,nOrb-nR

              eps = eQP(p) + eQP(i) - eQP(a) - eQP(b)
              reg = (1d0 - exp(- 2d0 * eta * eps * eps))
              num = 0.5d0*(ERI(p,i,b,a) - ERI(p,i,a,b))**2

              SigC(p) = SigC(p) + num*reg/eps
              Z(p)    = Z(p)    - num*reg/eps**2

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
  !$OMP PRIVATE(p,i,a,j,b,n,num,dem1,dem2,reg1,reg2) &
  !$OMP SHARED(nC,nO,nOrb,nR,nS,eta,ERI,eQP,eh_rho,eh_Om,SigC,Z)
  !$OMP DO COLLAPSE(2)
  do p=nC+1,nOrb-nR
     
     do i=nC+1,nO
        do a=nO+1,nOrb-nR
           
           do n=1,nS
              !3h2p
              do j=nC+1,nO
                 num  = ERI(p,a,j,i) * &
                 (eh_rho(j,p,n) * eh_rho(i,a,n) - eh_rho(a,p,n) * eh_rho(i,j,n))

                 dem1 = eQP(a) - eQP(i) - eh_Om(n) 
                 dem2 = eQP(p) - eQP(j) + eh_Om(n)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 SigC(p) = SigC(p) - num * (reg1/dem1) * (reg2/dem2)
                 Z(p)    = Z(p)    + num * (reg1/dem1) * (reg2/dem2/dem2)
                 
                 !num  = ERI(p,a,j,i) * &
                 !(eh_rho(j,p,n) * eh_rho(i,a,n) - eh_rho(a,p,n) * eh_rho(i,j,n))

                 !dem1 = eQP(a) - eQP(i) - eh_Om(n) 
                 dem2 = eQP(p) + eQP(a) - eQP(i) - eQP(j)
                 !reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 SigC(p) = SigC(p) + num * (reg1/dem1) * (reg2/dem2)
                 Z(p)    = Z(p)    - num * (reg1/dem1) * (reg2/dem2/dem2)
                                  
                 !num  = ERI(p,a,j,i) * &
                 !(eh_rho(j,p,n) * eh_rho(i,a,n) - eh_rho(a,p,n) * eh_rho(i,j,n))

                 dem1 = eQP(a) - eQP(i) + eh_Om(n) 
                 !dem2 = eQP(p) + eQP(a) - eQP(i) - eQP(j)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 !reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 SigC(p) = SigC(p) - num * (reg1/dem1) * (reg2/dem2)
                 Z(p)    = Z(p)    + num * (reg1/dem1) * (reg2/dem2/dem2)
                 
                 num  = ERI(p,i,j,a) * &
                 (eh_rho(j,p,n) * eh_rho(a,i,n) - eh_rho(i,p,n) * eh_rho(a,j,n))

                 !dem1 = eQP(a) - eQP(i) + eh_Om(n) 
                 dem2 = eQP(p) - eQP(j) + eh_Om(n)
                 !reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 SigC(p) = SigC(p) - num * (reg1/dem1) * (reg2/dem2)
                 Z(p)    = Z(p)    + num * (reg1/dem1) * (reg2/dem2/dem2)

                 
              end do ! j
              !3p2h
              do b=nO+1,nOrb-nR
                 num  = ERI(p,a,b,i) * &
                 (eh_rho(b,p,n) * eh_rho(i,a,n) - eh_rho(a,p,n) * eh_rho(i,b,n))

                 dem1 = eQP(a) - eQP(i) + eh_Om(n) 
                 dem2 = eQP(p) - eQP(b) - eh_Om(n)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 SigC(p) = SigC(p) - num * (reg1/dem1) * (reg2/dem2)
                 Z(p)    = Z(p)    + num * (reg1/dem1) * (reg2/dem2/dem2)
                 
                 num  = ERI(p,i,b,a) * &
                 (eh_rho(b,p,n) * eh_rho(a,i,n) - eh_rho(i,p,n) * eh_rho(a,b,n))

                 !dem1 = eQP(a) - eQP(i) + eh_Om(n) 
                 dem2 = eQP(p) + eQP(i) - eQP(a) - eQP(b)
                 !reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 SigC(p) = SigC(p) - num * (reg1/dem1) * (reg2/dem2)
                 Z(p)    = Z(p)    + num * (reg1/dem1) * (reg2/dem2/dem2)
                 
                 !num  = ERI(p,i,b,a) * &
                 !(eh_rho(b,p,n) * eh_rho(a,i,n) - eh_rho(b,i,n) * eh_rho(a,p,n))

                 dem1 = eQP(a) - eQP(i) - eh_Om(n) 
                 !dem2 = eQP(p) + eQP(i) - eQP(a) - eQP(b)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 !reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 SigC(p) = SigC(p) + num * (reg1/dem1) * (reg2/dem2)
                 Z(p)    = Z(p)    - num * (reg1/dem1) * (reg2/dem2/dem2)
                 
                 !num  = ERI(p,i,b,a) * &
                 !(eh_rho(b,p,n) * eh_rho(a,i,n) - eh_rho(b,i,n) * eh_rho(a,p,n))

                 !dem1 = eQP(a) - eQP(i) - eh_Om(n) 
                 dem2 = eQP(p) - eQP(b) - eh_Om(n)
                 !reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 SigC(p) = SigC(p) - num * (reg1/dem1) * (reg2/dem2)
                 Z(p)    = Z(p)    + num * (reg1/dem1) * (reg2/dem2/dem2)
                 
              end do ! b
              
           end do ! n
           
        end do ! a
     end do ! i
     
  end do ! p
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
  !$OMP PRIVATE(p,i,j,k,c,n,num,dem1,dem2,reg1,reg2) &
  !$OMP SHARED(nC,nO,nOrb,nR,nOO,nVV,eta,ERI,eQP,ee_rho,ee_Om,hh_rho,hh_Om,SigC,Z)
  !$OMP DO COLLAPSE(2)
  do p=nC+1,nOrb-nR
     
     do i=nC+1,nO
        do j=nC+1,nO
           do n=1,nVV
              ! 4h1p
              do k=nC+1,nO
                 num  = ERI(p,k,i,j) * ee_rho(i,j,n) * ee_rho(p,k,n)
                 dem1 = ee_Om(n) - eQP(i) - eQP(j)
                 dem2 = eQP(p) + eQP(k) - ee_Om(n)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 SigC(p) = SigC(p) + num * (reg1/dem1) * (reg2/dem2)
                 Z(p)    = Z(p)    - num * (reg1/dem1) * (reg2/dem2/dem2)

              end do ! k
              ! 3h2p
              do c=nO+1,nOrb-nR

                 num  = ERI(p,c,i,j) * ee_rho(i,j,n) * ee_rho(p,c,n)
                 !dem1 = ee_Om(n) - eQP(i) - eQP(j)
                 dem2 = eQP(p) + eQP(c) - eQP(i) - eQP(j)
                 !reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 SigC(p) = SigC(p) + num * (reg1/dem1) * (reg2/dem2)
                 Z(p)    = Z(p)    - num * (reg1/dem1) * (reg2/dem2/dem2)
                 
              end do ! a
           end do ! n
           do n=1,nOO
              ! 3h2p
              do c=nO+1,nOrb-nR

                 num  = ERI(p,c,i,j) * hh_rho(i,j,n) * hh_rho(p,c,n)
                 dem1 = hh_Om(n) - eQP(i) - eQP(j)
                 dem2 = eQP(p) + eQP(c) - hh_Om(n)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 SigC(p) = SigC(p) + num * (reg1/dem1) * (reg2/dem2)
                 Z(p)    = Z(p)    - num * (reg1/dem1) * (reg2/dem2/dem2)

!                num  = ERI(p,c,i,j) * hh_rho(i,j,n) * hh_rho(p,c,n)
                 !dem1 = hh_Om(n) - eQP(i) - eQP(j)
                 dem2 = eQP(p) + eQP(c) - eQP(i) - eQP(j)
                 !reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 SigC(p) = SigC(p) - num * (reg1/dem1) * (reg2/dem2)
                 Z(p)    = Z(p)    + num * (reg1/dem1) * (reg2/dem2/dem2)
                 
              end do ! c
           end do ! n
        end do ! j
     end do ! i
     
  end do ! p
  !$OMP END DO
  !$OMP END PARALLEL
  !$OMP PARALLEL DEFAULT(NONE)    &
  !$OMP PRIVATE(p,k,a,b,c,n,num,dem1,dem2,reg1,reg2) &
  !$OMP SHARED(nC,nO,nOrb,nR,nOO,nVV,eta,ERI,eQP,ee_rho,ee_Om,hh_rho,hh_Om,SigC,Z)
  !$OMP DO COLLAPSE(2)
  do p=nC+1,nOrb-nR
     do a=nO+1,nOrb-nR
        do b=nO+1,nOrb-nR
           do n=1,nOO
              ! 4p1h
              do c=nO+1,nOrb-nR
                 
                 num  = ERI(p,c,a,b) * hh_rho(a,b,n) * hh_rho(p,c,n)
                 dem1 = hh_Om(n) - eQP(a) - eQP(b)
                 dem2 = eQP(p) + eQP(c) - hh_Om(n)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 SigC(p) = SigC(p) - num * (reg1/dem1) * (reg2/dem2)
                 Z(p)    = Z(p)    + num * (reg1/dem1) * (reg2/dem2/dem2)

              end do ! c
              ! 3p2h
              do k=nC+1,nO

                 num  = ERI(p,k,a,b) * hh_rho(a,b,n) * hh_rho(p,k,n)
                 !dem1 = hh_Om(n) - eQP(a) - eQP(b)
                 dem2 = eQP(p) + eQP(k) - eQP(a) - eQP(b)
                 !reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 SigC(p) = SigC(p) - num * (reg1/dem1) * (reg2/dem2)
                 Z(p)    = Z(p)    + num * (reg1/dem1) * (reg2/dem2/dem2)
                 
              end do ! k
           end do ! n
           do n=1,nVV
              ! 3p2h
              do k=nC+1,nO

                 num  = ERI(p,k,a,b) * ee_rho(a,b,n) * ee_rho(p,k,n)
                 dem1 = ee_Om(n) - eQP(a) - eQP(b)
                 dem2 = eQP(p) + eQP(k) - ee_Om(n)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 SigC(p) = SigC(p) - num * (reg1/dem1) * (reg2/dem2)
                 Z(p)    = Z(p)    + num * (reg1/dem1) * (reg2/dem2/dem2)

!                num  = ERI(p,k,a,b) * ee_rho(a,b,n) * ee_rho(p,k,n)
                 !dem1 = ee_Om(n) - eQP(a) - eQP(b)
                 dem2 = eQP(p) + eQP(k) - eQP(a) - eQP(b)
                 !reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 SigC(p) = SigC(p) + num * (reg1/dem1) * (reg2/dem2)
                 Z(p)    = Z(p)    - num * (reg1/dem1) * (reg2/dem2/dem2)
                 
              end do ! c
           end do ! n
        end do ! b
     end do ! a
     
  end do ! p
  !$OMP END DO
  !$OMP END PARALLEL

  call wall_time(end_t)
  t = end_t - start_t

  write(*,'(1X,A50,1X,F9.3,A8)') 'Wall time for building pp self-energy =',t,' seconds'
  write(*,*)
  
!-----------------------------!
!   Renormalization factor    !
!-----------------------------!

  Z(:) = 1d0/(1d0 - Z(:))
  
!-------------------------------------!
! Galitskii-Migdal correlation energy !
!-------------------------------------!

  EcGM = 0d0
  
end subroutine 
