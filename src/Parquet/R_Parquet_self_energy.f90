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
  
! Output variables
  double precision,intent(out)  :: EcGM
  double precision,intent(out)  :: SigC(nOrb)
  double precision,intent(out)  :: Z(nOrb)

! Initialize 

  SigC(:) = 0d0
  Z(:)    = 0d0
  EcGM    = 0d0
  
!-----------------------------!
! GF2 part of the self-energy !
!-----------------------------!
  call wall_time(start_t)
  do p=nC+1,nOrb-nR
     ! 2h1p sum
     do i=nC+1,nO
        do j=nC+1,nO
           do a=nO+1,nOrb-nR

              eps = eQP(p) + eQP(a) - eQP(i) - eQP(j)
              reg = (1d0 - exp(- 2d0 * eta * eps * eps))
              num = ERI(p,a,j,i)*(2d0*ERI(j,i,p,a) - ERI(j,i,a,p))

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
              num = ERI(p,i,b,a)*(2d0*ERI(b,a,p,i) - ERI(b,a,i,p))

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
!   Renormalization factor    !
!-----------------------------!

  Z(:) = 1d0/(1d0 - Z(:))
  
!-------------------------------------!
! Galitskii-Migdal correlation energy !
!-------------------------------------!

  EcGM = 0d0

end subroutine 
