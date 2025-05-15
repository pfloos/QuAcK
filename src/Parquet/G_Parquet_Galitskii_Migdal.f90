subroutine G_Parquet_Galitskii_Migdal(eta,nOrb,nC,nO,nV,nR,nS,nOO,nVV,eQP,ERI,&
                                 eh_rho,eh_Om,ee_rho,ee_Om,hh_rho,hh_Om,EcGM)

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
  double precision,intent(in)   :: eh_rho(nOrb,nOrb,nS+nS)
  double precision,intent(in)   :: eh_Om(nS)
  double precision,intent(in)   :: ee_rho(nOrb,nOrb,nVV)
  double precision,intent(in)   :: ee_Om(nVV)
  double precision,intent(in)   :: hh_rho(nOrb,nOrb,nOO)
  double precision,intent(in)   :: hh_Om(nOO)

! Local variables
  integer                       :: i,j,k,l,a,b,c,d
  integer                       :: p,n
  double precision              :: eps,dem1,dem2,reg,reg1,reg2
  double precision              :: num
  double precision              :: start_t,end_t,t

! Output variables

  double precision,intent(out)  :: EcGM

! Initialize 
  EcGM    = 0d0
  
!-----------------------------------!
! 2nd-order part of the self-energy !
!-----------------------------------!

  call wall_time(start_t)

  do i=nC+1,nO
     do j=nC+1,nO
        do a=nO+1,nOrb-nR
           do b=nO+1,nOrb-nR

              eps = eQP(i) + eQP(j) - eQP(a) - eQP(b)
              reg = (1d0 - exp(- 2d0 * eta * eps * eps))
              num = 0.5d0*(ERI(a,b,i,j) - ERI(a,b,j,i))**2

              EcGM = EcGM + num*reg/eps

           end do
        end do
     end do
  end do

  call wall_time(end_t)
  t = end_t - start_t

  write(*,'(1X,A50,1X,F9.3,A8)') 'Wall time for computing GF(2) correlation energy =',t,' seconds'
  write(*,*)

!-----------------------------!
!  eh part of the self-energy !
!-----------------------------!

  ! call wall_time(start_t)

  ! ! !$OMP PARALLEL DEFAULT(NONE)    &
  ! ! !$OMP PRIVATE(p,i,a,j,b,n,num,dem1,dem2,reg1,reg2) &
  ! ! !$OMP SHARED(nC,nO,nOrb,nR,nS,eta,ERI,eQP,eh_rho,eh_Om,SigC,Z)
  ! ! !$OMP DO COLLAPSE(2)
  ! do i=nC+1,nO
  !    do j=nC+1,nO
  !       do a=nO+1,nOrb-nR
  !          do b=nO+1,nOrb-nR
  !             do n=1,nS

  !                num  = - 0.5d0 * (ERI(b,a,i,j) - ERI(b,a,j,i)) * eh_rho(i,a,n) * eh_rho(j,b,nS+n)
                 
  !                dem1 = eQP(b) - eQP(j) + eh_Om(n)
  !                dem2 = eQP(a) + eQP(b) - eQP(i) - eQP(j)
  !                reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
  !                reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

  !                EcGM = EcGM + num * (reg1/dem1) * (reg2/dem2)
                 
  !                num  = 0.5d0 * (ERI(j,i,a,b) - ERI(j,i,b,a)) * eh_rho(a,i,n) * eh_rho(b,j,nS+n)

  !                dem1 = eQP(a) - eQP(i) + eh_Om(n)
  !                dem2 = eQP(i) + eQP(j) - eQP(a) - eQP(b)
  !                reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
  !                reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))
 
  !                EcGM = EcGM + num * (reg1/dem1) * (reg2/dem2)
                 
  !                num  = 0.5d0 * (ERI(j,a,i,b) - ERI(j,a,b,i)) * eh_rho(i,a,nS+n) * eh_rho(b,j,n)

  !                dem1 = eQP(a) - eQP(i) + eh_Om(n) 
  !                dem2 = eQP(j) - eQP(b) - eh_Om(n)
  !                reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
  !                reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

  !                EcGM = EcGM + num * (reg1/dem1) * (reg2/dem2)

  !                num  = - 0.5d0 * (ERI(j,i,a,b) - ERI(j,i,b,a)) * eh_rho(a,i,nS+n) * eh_rho(b,j,n)

  !                dem1 = eQP(j) - eQP(b) - eh_Om(n) 
  !                dem2 = eQP(i) + eQP(j) - eQP(a) - eQP(b)
  !                reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
  !                reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

  !                EcGM = EcGM + num * (reg1/dem1) * (reg2/dem2)

  !                num  = 0.5d0 * (ERI(b,a,i,j) - ERI(b,a,j,i)) * eh_rho(i,a,nS+n) * eh_rho(j,b,n)

  !                dem1 = eQP(a) - eQP(i) + eh_Om(n) 
  !                dem2 = eQP(i) + eQP(j) - eQP(a) - eQP(b)
  !                reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
  !                reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

  !                EcGM = EcGM + num * (reg1/dem1) * (reg2/dem2)
                 
  !                num  = 0.5d0 * (ERI(b,i,a,j) - ERI(b,i,j,a)) * eh_rho(a,i,n) * eh_rho(j,b,nS+n)

  !                dem1 = eQP(a) - eQP(i) + eh_Om(n) 
  !                dem2 = eQP(j) - eQP(b) - eh_Om(n)
  !                reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
  !                reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

  !                EcGM = EcGM + num * (reg1/dem1) * (reg2/dem2)
                 
  !             end do ! n
  !          end do ! b
  !       end do ! a
  !    end do ! j
  ! end do ! i
  ! ! !$OMP END DO
  ! ! !$OMP END PARALLEL

  ! call wall_time(end_t)
  ! t = end_t - start_t

  ! write(*,'(1X,A50,1X,F9.3,A8)') 'Wall time for computing eh correlation energy =',t,' seconds'
  ! write(*,*) 
  
!-----------------------------!
!  pp part of the self-energy !
!-----------------------------!

  call wall_time(start_t)

  ! !$OMP PARALLEL DEFAULT(NONE)    &
  ! !$OMP PRIVATE(p,i,j,k,l,a,b,n,num,dem1,dem2,reg1,reg2) &
  ! !$OMP SHARED(nC,nO,nOrb,nR,nOO,nVV,eta,ERI,eQP,ee_rho,ee_Om,hh_rho,hh_Om,EcGM)
  ! !$OMP DO COLLAPSE(2)
  do i=nC+1,nO
     do j=nC+1,nO
        do n=1,nVV

           do k=nC+1,nO
              do l=nC+1,nO
                 num  = - 0.5d0 * ERI(k,l,i,j) * ee_rho(i,j,n) * ee_rho(k,l,n)
                 dem1 = ee_Om(n) - eQP(i) - eQP(j)
                 dem2 = eQP(k) + eQP(l) - ee_Om(n)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 EcGM = EcGM + num * (reg1/dem1) * (reg2/dem2)

              end do ! l
           end do ! k
           
           do a=nO+1,nOrb-nR
              do b=nO+1,nOrb-nR

                 num  = - 0.5d0 * ERI(a,b,i,j) * ee_rho(i,j,n) * ee_rho(a,b,n)
                 dem1 = ee_Om(n) - eQP(i) - eQP(j)
                 dem2 = eQP(i) + eQP(j) - eQP(a) - eQP(b)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 EcGM = EcGM + num * (reg1/dem1) * (reg2/dem2)
                 
              end do ! a
           end do ! a
           
        end do ! n
        do n=1,nOO
           
           do a=nO+1,nOrb-nR
              do b=nO+1,nOrb-nR

                 num  = - 0.5d0 * ERI(a,b,i,j) * hh_rho(i,j,n) * hh_rho(a,b,n)
                 dem1 = eQP(a) + eQP(b) - hh_Om(n)
                 dem2 = eQP(i) + eQP(j) - eQP(a) - eQP(b)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 EcGM = EcGM + num * (reg1/dem1) * (reg2/dem2)
                 
              end do ! a
           end do ! a
           
        end do ! n
     end do ! j
  end do ! i

  ! !$OMP END DO
  ! !$OMP END PARALLEL
  ! !$OMP PARALLEL DEFAULT(NONE)    &
  ! !$OMP PRIVATE(p,k,a,b,c,n,num,dem1,dem2,reg1,reg2) &
  ! !$OMP SHARED(nC,nO,nOrb,nR,nOO,nVV,eta,ERI,eQP,ee_rho,ee_Om,hh_rho,hh_Om)
  ! !$OMP DO COLLAPSE(2)
  
  do a=nO+1,nOrb-nR
     do b=nO+1,nOrb-nR
        do n=1,nOO
           
           do c=nO+1,nOrb-nR
              do d=nO+1,nOrb-nR
                 
                 num  = 0.5d0 * ERI(c,d,a,b) * hh_rho(a,b,n) * hh_rho(c,d,n)
                 dem1 = hh_Om(n) - eQP(a) - eQP(b)
                 dem2 = hh_Om(n) - eQP(c) - eQP(d)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 EcGM = EcGM + num * (reg1/dem1) * (reg2/dem2)

              end do ! d
           end do ! c
           
           do i=nC+1,nO
              do j=nC+1,nO

                 num  = 0.5d0 * ERI(i,j,a,b) * hh_rho(a,b,n) * hh_rho(i,j,n)
                 dem1 = hh_Om(n) - eQP(a) - eQP(b)
                 dem2 = eQP(i) + eQP(j) - eQP(a) - eQP(b)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 EcGM = EcGM + num * (reg1/dem1) * (reg2/dem2)
                 
              end do ! j
           end do ! i
           
        end do ! n
        do n=1,nVV

           do i=nC+1,nO
              do j=nC+1,nO

                 num  = 0.5d0 * ERI(i,j,a,b) * ee_rho(a,b,n) * ee_rho(i,j,n)
                 dem1 = eQP(i) + eQP(j) - eQP(a) - eQP(b)
                 dem2 = eQP(i) + eQP(j) - ee_Om(n)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 EcGM = EcGM + num * (reg1/dem1) * (reg2/dem2)
                 
              end do ! j
           end do ! i
              
        end do ! n
     end do ! b
  end do ! a

  ! !$OMP END DO
  ! !$OMP END PARALLEL

  call wall_time(end_t)
  t = end_t - start_t

  write(*,'(1X,A50,1X,F9.3,A8)') 'Wall time for computing pp correlation energy =',t,' seconds'
  write(*,*)

  write(*,*) 'The correlation energy is', EcGM
  write(*,*)
  
end subroutine 
