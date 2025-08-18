subroutine R_Parquet_self_energy_interm(eta,nOrb,nC,nO,nV,nR,nS,nOOs,nVVs,nOOt,nVVt,eQP,ERI, &
                                 eh_sing_rho,eh_sing_Om,eh_trip_rho,eh_trip_Om, &
                                 ee_sing_rho,ee_sing_Om,ee_trip_rho,ee_trip_Om, &
                                 hh_sing_rho,hh_sing_Om,hh_trip_rho,hh_trip_Om, &
                                 SigC,Z)

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
  integer                       :: p,q,n
  double precision              :: eps_p,eps_q,eps_pq,reg
  double precision              :: dem1,dem2,reg1,reg2
  double precision              :: dem2_p,dem2_q,dem2_pq
  double precision              :: num
  double precision              :: start_t,end_t,t

  double precision,allocatable  :: int3(:,:,:)
  double precision,allocatable  :: int4(:,:,:,:)
  
! Output variables
  double precision,intent(out)  :: SigC(nOrb,nOrb)
  double precision,intent(out)  :: Z(nOrb)

! Initialize 

  SigC(:,:) = 0d0

!-----------------------------------!
! 2nd-order part of the self-energy !
!-----------------------------------!
  call wall_time(start_t)
  do q=nC+1,nOrb-nR
    do p=nC+1,nOrb-nR
      ! 2h1p sum
      do i=nC+1,nO
        do j=nC+1,nO
          do a=nO+1,nOrb-nR
            num  = ERI(p,a,j,i) * ( 2d0*ERI(j,i,q,a) - ERI(j,i,a,q) )

            eps_p = eQP(p) + eQP(a) - eQP(i) - eQP(j)
            eps_q = eQP(q) + eQP(a) - eQP(i) - eQP(j)
            eps_pq = eps_p * eps_p + eps_q * eps_q
            reg = 1d0 - exp(- eta * eps_pq)
            
            SigC(p,q) = SigC(p,q) + num*reg*(eps_p + eps_q)/eps_pq
            
          end do ! a
        end do ! j
      end do ! i
      ! 2p1h sum
      do i=nC+1,nO
        do a=nO+1,nOrb-nR
          do b=nO+1,nOrb-nR
            num  = ERI(p,i,b,a) * ( 2d0*ERI(b,a,q,i) - ERI(b,a,i,q) )

            eps_p = eQP(p) + eQP(i) - eQP(a) - eQP(b)
            eps_q = eQP(q) + eQP(i) - eQP(a) - eQP(b)
            eps_pq = eps_p * eps_p + eps_q * eps_q
            reg = 1d0 - exp(- eta * eps_pq)
            
            SigC(p,q) = SigC(p,q) + num*reg*(eps_p + eps_q)/eps_pq
                 
          end do ! b
        end do ! a
      end do ! i
        
    end do ! p
  end do ! q
  call wall_time(end_t)
  t = end_t - start_t
  
  write(*,'(1X,A50,1X,F9.3,A8)') 'Wall time for 2nd-order self-energy =',t,' seconds'
  write(*,*)

!-------------------------------------!
!  singlet eh part of the self-energy !
!-------------------------------------!
  call wall_time(start_t)
!-------------------------------------!
! Intermediates for singlet eh part 1 !
!-------------------------------------!

  allocate(int3(nOrb,nOrb,nS))
  int3(:,:,:) = 0d0

  do p=nC+1,nOrb-nR
    do j=nC+1,nO
      do n=1,nS

        do i=nC+1,nO
          do a=nO+1,nOrb-nR

            num  = (0.5d0*ERI(p,a,i,j) - ERI(p,a,j,i)) * eh_sing_rho(i,a,n)
            dem1 = eQP(a) - eQP(i) - eh_sing_Om(n)
            reg1 = 1d0 - exp(- 2d0 * eta * dem1 * dem1)

            int3(p,j,n) = int3(p,j,n) + num*(reg1/dem1) 

            num  = (0.5d0*ERI(p,i,a,j) - ERI(p,i,j,a)) * eh_sing_rho(a,i,n)
            dem1 = eQP(a) - eQP(i) + eh_sing_Om(n) 
            reg1 = 1d0 - exp(- 2d0 * eta * dem1 * dem1)

            int3(p,j,n) = int3(p,j,n) + num*(reg1/dem1) 

          end do
        end do

      end do
    end do
  end do

  do q=nC+1,nOrb-nR
    do p=nC+1,nOrb-nR
      do j=nC+1,nO
        do n=1,nS

          num  = eh_sing_rho(q,j,n)
          dem2_p = eQP(p) - eQP(j) + eh_sing_Om(n)
          dem2_q = eQP(q) - eQP(j) + eh_sing_Om(n)
          dem2_pq = dem2_p * dem2_p + dem2_q * dem2_q
          reg2 = 1d0 - exp(- eta * dem2_pq)
             
          SigC(p,q) = SigC(p,q) + num * int3(p,j,n) * reg2 * (dem2_p + dem2_q)/dem2_pq

        end do
      end do
    end do
  end do

  deallocate(int3)

!-------------------------------------!
! Intermediates for singlet eh part 2 !
!-------------------------------------!

  allocate(int4(nOrb,nOrb,nOrb,nOrb))
  int4(:,:,:,:) = 0d0

  do q=nC+1,nOrb-nR
    do i=nC+1,nO
      do a=nO+1,nOrb-nR
        do j=nC+1,nO

          do n=1,nS

            num  = - eh_sing_rho(i,a,n) * eh_sing_rho(q,j,n)
            dem1 = eQP(a) - eQP(i) - eh_sing_Om(n)
            reg1 = 1d0 - exp(- 2d0 * eta * dem1 * dem1)

            int4(q,i,j,a) = int4(q,i,j,a) + num*(reg1/dem1) 

            num  = eh_sing_rho(a,i,n) * eh_sing_rho(j,q,n)
            dem1 = eQP(a) - eQP(i) + eh_sing_Om(n)
            reg1 = 1d0 - exp(- 2d0 * eta * dem1 * dem1)

            int4(q,i,j,a) = int4(q,i,j,a) + num*(reg1/dem1) 

          end do

        end do
      end do
    end do
  end do

  do q=nC+1,nOrb-nR
    do p=nC+1,nOrb-nR
      do i=nC+1,nO
        do a=nO+1,nOrb-nR
          do j=nC+1,nO

            num  = 0.5d0*ERI(p,a,i,j) - ERI(p,a,j,i)
            dem2_p = eQP(p) - eQP(i) - eQP(j) + eQP(a)
            dem2_q = eQP(q) - eQP(i) - eQP(j) + eQP(a)
            dem2_pq = dem2_p * dem2_p + dem2_q * dem2_q
            reg2 = 1d0 - exp(- eta * dem2_pq)
               
            SigC(p,q) = SigC(p,q) + num * int4(q,i,j,a) * reg2 * (dem2_p + dem2_q)/dem2_pq

          end do
        end do
      end do
    end do
  end do

  deallocate(int4)

!-------------------------------------!
! Intermediates for singlet eh part 3 !
!-------------------------------------!

  allocate(int3(nOrb,nOrb,nS))
  int3(:,:,:) = 0d0

  do p=nC+1,nOrb-nR
    do b=nO+1,nOrb-nR
      do n=1,nS

        do i=nC+1,nO
          do a=nO+1,nOrb-nR

            num  = (0.5d0*ERI(p,i,a,b) - ERI(p,i,b,a)) * eh_sing_rho(i,a,n)
            dem1 = eQP(a) - eQP(i) - eh_sing_Om(n)
            reg1 = 1d0 - exp(- 2d0 * eta * dem1 * dem1)

            int3(p,b,n) = int3(p,b,n) + num*(reg1/dem1) 

            num  = (0.5d0*ERI(p,a,i,b) - ERI(p,a,b,i)) * eh_sing_rho(a,i,n) 
            dem1 = eQP(a) - eQP(i) + eh_sing_Om(n)
            reg1 = 1d0 - exp(- 2d0 * eta * dem1 * dem1)

            int3(p,b,n) = int3(p,b,n) + num*(reg1/dem1) 

          end do
        end do

      end do
    end do
  end do

  do q=nC+1,nOrb-nR
    do p=nC+1,nOrb-nR
      do b=nO+1,nOrb-nR
        do n=1,nS

          num  = eh_sing_rho(b,q,n)
          dem2_p = eQP(p) - eQP(b) - eh_sing_Om(n)
          dem2_q = eQP(q) - eQP(b) - eh_sing_Om(n)
          dem2_pq = dem2_p * dem2_p + dem2_q * dem2_q
          reg2 = 1d0 - exp(- eta * dem2_pq)

          SigC(p,q) = SigC(p,q) + num * int3(p,b,n) * reg2 * (dem2_p + dem2_q)/dem2_pq

        end do
      end do
    end do
  end do

  deallocate(int3)

!-------------------------------------!
! Intermediates for singlet eh part 4 !
!-------------------------------------!

  allocate(int4(nOrb,nOrb,nOrb,nOrb))
  int4(:,:,:,:) = 0d0

  do q=nC+1,nOrb-nR
    do i=nC+1,nO
      do a=nO+1,nOrb-nR
        do b=nO+1,nOrb-nR

          do n=1,nS

            num  = - eh_sing_rho(i,a,n) * eh_sing_rho(b,q,n)
            dem1 = eQP(a) - eQP(i) - eh_sing_Om(n)
            reg1 = 1d0 - exp(- 2d0 * eta * dem1 * dem1)

            int4(q,i,a,b) = int4(q,i,a,b) + num*(reg1/dem1) 

            num  = eh_sing_rho(a,i,n) * eh_sing_rho(q,b,n)
            dem1 = eQP(a) - eQP(i) + eh_sing_Om(n)
            reg1 = 1d0 - exp(- 2d0 * eta * dem1 * dem1)

            int4(q,i,a,b) = int4(q,i,a,b) + num*(reg1/dem1) 

          end do

        end do
      end do
    end do
  end do

  do q=nC+1,nOrb-nR
    do p=nC+1,nOrb-nR
      do i=nC+1,nO
        do a=nO+1,nOrb-nR
          do b=nO+1,nOrb-nR

            num  = 0.5d0*ERI(p,i,a,b) - ERI(p,i,b,a)
            dem2_p = eQP(p) + eQP(i) - eQP(a) - eQP(b)
            dem2_q = eQP(q) + eQP(i) - eQP(a) - eQP(b)
            dem2_pq = dem2_p * dem2_p + dem2_q * dem2_q
            reg2 = 1d0 - exp(- eta * dem2_pq)
               
            SigC(p,q) = SigC(p,q) + num * int4(q,i,a,b) * reg2 * (dem2_p + dem2_q)/dem2_pq

          end do
        end do
      end do
    end do
  end do

  deallocate(int4)

  call wall_time(end_t)
  t = end_t - start_t

  write(*,'(1X,A50,1X,F9.3,A8)') 'Wall time for singlet eh self-energy =',t,' seconds'
  write(*,*)

!-------------------------------------!
!  triplet eh part of the self-energy !
!-------------------------------------!
  call wall_time(start_t)
!-------------------------------------!
! Intermediates for triplet eh part 1 !
!-------------------------------------!

  allocate(int3(nOrb,nOrb,nS))
  int3(:,:,:) = 0d0

  do p=nC+1,nOrb-nR
    do j=nC+1,nO
      do n=1,nS

        do i=nC+1,nO
          do a=nO+1,nOrb-nR

            num  = 1.5d0*ERI(p,a,i,j) * eh_trip_rho(i,a,n)
            dem1 = eQP(a) - eQP(i) - eh_trip_Om(n)
            reg1 = 1d0 - exp(- 2d0 * eta * dem1 * dem1)

            int3(p,j,n) = int3(p,j,n) + num*(reg1/dem1) 

            num  = 1.5d0*ERI(p,i,a,j) * eh_trip_rho(a,i,n)
            dem1 = eQP(a) - eQP(i) + eh_trip_Om(n) 
            reg1 = 1d0 - exp(- 2d0 * eta * dem1 * dem1)

            int3(p,j,n) = int3(p,j,n) + num*(reg1/dem1) 

          end do
        end do

      end do
    end do
  end do

  do q=nC+1,nOrb-nR
    do p=nC+1,nOrb-nR
      do j=nC+1,nO
        do n=1,nS

          num  = eh_trip_rho(q,j,n)
          dem2_p = eQP(p) - eQP(j) + eh_trip_Om(n)
          dem2_q = eQP(q) - eQP(j) + eh_trip_Om(n)
          dem2_pq = dem2_p * dem2_p + dem2_q * dem2_q
          reg2 = 1d0 - exp(- eta * dem2_pq)
             
          SigC(p,q) = SigC(p,q) + num * int3(p,j,n) * reg2 * (dem2_p + dem2_q)/dem2_pq

        end do
      end do
    end do
  end do

  deallocate(int3)

!-------------------------------------!
! Intermediates for triplet eh part 2 !
!-------------------------------------!

  allocate(int4(nOrb,nOrb,nOrb,nOrb))
  int4(:,:,:,:) = 0d0

  do q=nC+1,nOrb-nR
    do i=nC+1,nO
      do a=nO+1,nOrb-nR
        do j=nC+1,nO

          do n=1,nS

            num  = - eh_trip_rho(i,a,n) * eh_trip_rho(q,j,n)
            dem1 = eQP(a) - eQP(i) - eh_trip_Om(n)
            reg1 = 1d0 - exp(- 2d0 * eta * dem1 * dem1)

            int4(q,i,j,a) = int4(q,i,j,a) + num*(reg1/dem1) 

            num  = eh_trip_rho(a,i,n) * eh_trip_rho(j,q,n)
            dem1 = eQP(a) - eQP(i) + eh_trip_Om(n)
            reg1 = 1d0 - exp(- 2d0 * eta * dem1 * dem1)

            int4(q,i,j,a) = int4(q,i,j,a) + num*(reg1/dem1) 

          end do

        end do
      end do
    end do
  end do

  do q=nC+1,nOrb-nR
    do p=nC+1,nOrb-nR
      do i=nC+1,nO
        do a=nO+1,nOrb-nR
          do j=nC+1,nO

            num  = 1.5d0*ERI(p,a,i,j)
            dem2_p = eQP(p) - eQP(i) - eQP(j) + eQP(a)
            dem2_q = eQP(q) - eQP(i) - eQP(j) + eQP(a)
            dem2_pq = dem2_p * dem2_p + dem2_q * dem2_q
            reg2 = 1d0 - exp(- eta * dem2_pq)
               
            SigC(p,q) = SigC(p,q) + num * int4(q,i,j,a) * reg2 * (dem2_p + dem2_q)/dem2_pq

          end do
        end do
      end do
    end do
  end do

  deallocate(int4)

!-------------------------------------!
! Intermediates for triplet eh part 3 !
!-------------------------------------!

  allocate(int3(nOrb,nOrb,nS))
  int3(:,:,:) = 0d0

  do p=nC+1,nOrb-nR
    do b=nO+1,nOrb-nR
      do n=1,nS

        do i=nC+1,nO
          do a=nO+1,nOrb-nR

            num  = 1.5d0*ERI(p,i,a,b) * eh_trip_rho(i,a,n)
            dem1 = eQP(a) - eQP(i) - eh_trip_Om(n)
            reg1 = 1d0 - exp(- 2d0 * eta * dem1 * dem1)

            int3(p,b,n) = int3(p,b,n) + num*(reg1/dem1) 

            num  = 1.5d0*ERI(p,a,i,b) * eh_trip_rho(a,i,n) 
            dem1 = eQP(a) - eQP(i) + eh_trip_Om(n)
            reg1 = 1d0 - exp(- 2d0 * eta * dem1 * dem1)

            int3(p,b,n) = int3(p,b,n) + num*(reg1/dem1) 

          end do
        end do

      end do
    end do
  end do

  do q=nC+1,nOrb-nR
    do p=nC+1,nOrb-nR
      do b=nO+1,nOrb-nR
        do n=1,nS

          num  = eh_trip_rho(b,q,n)
          dem2_p = eQP(p) - eQP(b) - eh_trip_Om(n)
          dem2_q = eQP(q) - eQP(b) - eh_trip_Om(n)
          dem2_pq = dem2_p * dem2_p + dem2_q * dem2_q
          reg2 = 1d0 - exp(- eta * dem2_pq)

          SigC(p,q) = SigC(p,q) + num * int3(p,b,n) * reg2 * (dem2_p + dem2_q)/dem2_pq

        end do
      end do
    end do
  end do

  deallocate(int3)

!-------------------------------------!
! Intermediates for triplet eh part 4 !
!-------------------------------------!

  allocate(int4(nOrb,nOrb,nOrb,nOrb))
  int4(:,:,:,:) = 0d0

  do q=nC+1,nOrb-nR
    do i=nC+1,nO
      do a=nO+1,nOrb-nR
        do b=nO+1,nOrb-nR

          do n=1,nS

            num  = - eh_trip_rho(i,a,n) * eh_trip_rho(b,q,n)
            dem1 = eQP(a) - eQP(i) - eh_trip_Om(n)
            reg1 = 1d0 - exp(- 2d0 * eta * dem1 * dem1)

            int4(q,i,a,b) = int4(q,i,a,b) + num*(reg1/dem1) 

            num  = eh_trip_rho(a,i,n) * eh_trip_rho(q,b,n)
            dem1 = eQP(a) - eQP(i) + eh_trip_Om(n)
            reg1 = 1d0 - exp(- 2d0 * eta * dem1 * dem1)

            int4(q,i,a,b) = int4(q,i,a,b) + num*(reg1/dem1) 

          end do

        end do
      end do
    end do
  end do

  do q=nC+1,nOrb-nR
    do p=nC+1,nOrb-nR
      do i=nC+1,nO
        do a=nO+1,nOrb-nR
          do b=nO+1,nOrb-nR

            num  = 1.5d0*ERI(p,i,a,b)
            dem2_p = eQP(p) + eQP(i) - eQP(a) - eQP(b)
            dem2_q = eQP(q) + eQP(i) - eQP(a) - eQP(b)
            dem2_pq = dem2_p * dem2_p + dem2_q * dem2_q
            reg2 = 1d0 - exp(- eta * dem2_pq)
               
            SigC(p,q) = SigC(p,q) + num * int4(q,i,a,b) * reg2 * (dem2_p + dem2_q)/dem2_pq

          end do
        end do
      end do
    end do
  end do

  deallocate(int4)

  call wall_time(end_t)
  t = end_t - start_t

  write(*,'(1X,A50,1X,F9.3,A8)') 'Wall time for triplet eh self-energy =',t,' seconds'
  write(*,*)
  
!-------------------------------------!
!  singlet pp part of the self-energy !
!-------------------------------------!
  call wall_time(start_t)
!-------------------------------------!
! Intermediates for singlet pp part 1 !
!-------------------------------------!

  allocate(int3(nOrb,nOrb,nVVs))
  int3(:,:,:) = 0d0

  do p=nC+1,nOrb-nR
    do k=nC+1,nO
      do n=1,nVVs

        do i=nC+1,nO
          do j=nC+1,nO

            num  = - 0.5d0 * ERI(p,k,i,j) * ee_sing_rho(i,j,n)
            dem1 = ee_sing_Om(n) - eQP(i) - eQP(j)
            reg1 = 1d0 - exp(- 2d0 * eta * dem1 * dem1)

            int3(p,k,n) = int3(p,k,n) + num*(reg1/dem1)

          end do
        end do

        do a=nO+1,nOrb-nR
          do b=nO+1,nOrb-nR

            num  = - 0.5d0 * ERI(p,k,a,b) * ee_sing_rho(a,b,n)
            dem1 = eQP(a) + eQP(b) - ee_sing_Om(n)
            reg1 = 1d0 - exp(- 2d0 * eta * dem1 * dem1)

            int3(p,k,n) = int3(p,k,n) + num*(reg1/dem1)

          end do
        end do

      end do
    end do
  end do

  do q=nC+1,nOrb-nR
    do p=nC+1,nOrb-nR
      do k=nC+1,nO
        do n=1,nVVs

          num  = ee_sing_rho(q,k,n)
          dem2_p = eQP(p) + eQP(k) - ee_sing_Om(n)
          dem2_q = eQP(q) + eQP(k) - ee_sing_Om(n)
          dem2_pq = dem2_p * dem2_p + dem2_q * dem2_q
          reg2 = 1d0 - exp(- eta * dem2_pq)

          SigC(p,q) = SigC(p,q) + num * int3(p,k,n) * reg2 * (dem2_p + dem2_q)/dem2_pq

        end do
      end do
    end do
  end do

  deallocate(int3)

!-------------------------------------!
! Intermediates for singlet pp part 2 !
!-------------------------------------!

  allocate(int4(nOrb,nOrb,nOrb,nOrb))
  int4(:,:,:,:) = 0d0

  do q=nC+1,nOrb-nR
    do i=nC+1,nO
      do j=nC+1,nO
        do a=nO+1,nOrb-nR

          do n=1,nVVs

            num  = ee_sing_rho(i,j,n) * ee_sing_rho(q,a,n)
            dem1 = ee_sing_Om(n) - eQP(i) - eQP(j)
            reg1 = 1d0 - exp(- 2d0 * eta * dem1 * dem1)

            int4(q,i,j,a) = int4(q,i,j,a) + num*(reg1/dem1)

          end do

          do n=1,nOOs

            num  = hh_sing_rho(i,j,n) * hh_sing_rho(q,a,n)
            dem1 = eQP(i) + eQP(j) - hh_sing_Om(n)
            reg1 = 1d0 - exp(- 2d0*eta * dem1 * dem1)

            int4(q,i,j,a) = int4(q,i,j,a) + num*(reg1/dem1)

          end do

        end do
      end do
    end do
  end do

  do q=nC+1,nOrb-nR
    do p=nC+1,nOrb-nR
      do i=nC+1,nO
        do j=nC+1,nO
          do a=nO+1,nOrb-nR

            num  = - 0.5d0*ERI(p,a,i,j)
            dem2_p = eQP(p) + eQP(a) - eQP(i) - eQP(j)
            dem2_q = eQP(q) + eQP(a) - eQP(i) - eQP(j)
            dem2_pq = dem2_p * dem2_p + dem2_q * dem2_q
            reg2 = 1d0 - exp(- eta * dem2_pq)

            SigC(p,q) = SigC(p,q) + num * int4(q,i,j,a) * reg2 * (dem2_p + dem2_q)/dem2_pq

          end do
        end do
      end do
    end do
  end do

  deallocate(int4)

!-------------------------------------!
! Intermediates for singlet pp part 3 !
!-------------------------------------!

  allocate(int3(nOrb,nOrb,nOOs))
  int3(:,:,:) = 0d0

  do p=nC+1,nOrb-nR
    do c=nO+1,nOrb-nR
      do n=1,nOOs

        do i=nC+1,nO
          do j=nC+1,nO

            num  = 0.5d0*ERI(p,c,i,j) * hh_sing_rho(i,j,n)
            dem1 = eQP(i) + eQP(j) - hh_sing_Om(n)
            reg1 = 1d0 - exp(- 2d0*eta * dem1 * dem1)

            int3(p,c,n) = int3(p,c,n) + num*(reg1/dem1)

          end do
        end do

        do a=nO+1,nOrb-nR
          do b=nO+1,nOrb-nR

            num  = 0.5d0*ERI(p,c,a,b) * hh_sing_rho(a,b,n)
            dem1 = hh_sing_Om(n) - eQP(a) - eQP(b)
            reg1 = 1d0 - exp(- 2d0*eta * dem1 * dem1)

            int3(p,c,n) = int3(p,c,n) + num*(reg1/dem1)


          end do
        end do

      end do
    end do
  end do

  do q=nC+1,nOrb-nR
    do p=nC+1,nOrb-nR
      do c=nO+1,nOrb-nR
        do n=1,nOOs

          num  = hh_sing_rho(q,c,n)
          dem2_p = eQP(p) + eQP(c) - hh_sing_Om(n)
          dem2_q = eQP(q) + eQP(c) - hh_sing_Om(n)
          dem2_pq = dem2_p * dem2_p + dem2_q * dem2_q
          reg2 = 1d0 - exp(- eta * dem2_pq)

          SigC(p,q) = SigC(p,q) + num * int3(p,c,n) * reg2 * (dem2_p + dem2_q)/dem2_pq

        end do
      end do
    end do
  end do

  deallocate(int3)

!-------------------------------------!
! Intermediates for singlet pp part 4 !
!-------------------------------------!

  allocate(int4(nOrb,nOrb,nOrb,nOrb))
  int4(:,:,:,:) = 0d0

  do q=nC+1,nOrb-nR
    do i=nC+1,nO
      do a=nO+1,nOrb-nR
        do b=nO+1,nOrb-nR

          do n=1,nOOs

            num  = hh_sing_rho(a,b,n) * hh_sing_rho(q,i,n)
            dem1 = hh_sing_Om(n) - eQP(a) - eQP(b)
            reg1 = 1d0 - exp(- 2d0 * eta * dem1 * dem1)

            int4(q,i,a,b) = int4(q,i,a,b) + num*(reg1/dem1)

          end do

          do n=1,nVVs

            num  = ee_sing_rho(a,b,n) * ee_sing_rho(q,i,n)
            dem1 = eQP(a) + eQP(b) - ee_sing_Om(n)
            reg1 = 1d0 - exp(- 2d0 * eta * dem1 * dem1)

            int4(q,i,a,b) = int4(q,i,a,b) + num*(reg1/dem1)

          end do

        end do
      end do
    end do
  end do

  do q=nC+1,nOrb-nR
    do p=nC+1,nOrb-nR
      do i=nC+1,nO
        do a=nO+1,nOrb-nR
          do b=nO+1,nOrb-nR

            num  = 0.5d0*ERI(p,i,a,b)
            dem2_p = eQP(p) + eQP(i) - eQP(a) - eQP(b)
            dem2_q = eQP(q) + eQP(i) - eQP(a) - eQP(b)
            dem2_pq = dem2_p * dem2_p + dem2_q * dem2_q
            reg2 = 1d0 - exp(- eta * dem2_pq)

            SigC(p,q) = SigC(p,q) + num * int4(q,i,a,b) * reg2 * (dem2_p + dem2_q)/dem2_pq

          end do
        end do
      end do
    end do
  end do

  deallocate(int4)

  call wall_time(end_t)
  t = end_t - start_t

  write(*,'(1X,A50,1X,F9.3,A8)') 'Wall time for singlet pp self-energy =',t,' seconds'
  write(*,*)

!-------------------------------------!
!  triplet pp part of the self-energy !
!-------------------------------------!
  call wall_time(start_t)
!-------------------------------------!
! Intermediates for triplet pp part 1 !
!-------------------------------------!

  allocate(int3(nOrb,nOrb,nVVt))
  int3(:,:,:) = 0d0

  do p=nC+1,nOrb-nR
    do k=nC+1,nO
      do n=1,nVVt

        do i=nC+1,nO
          do j=nC+1,nO

            num  = - 1.5d0 * ERI(p,k,i,j) * ee_trip_rho(i,j,n)
            dem1 = ee_trip_Om(n) - eQP(i) - eQP(j)
            reg1 = 1d0 - exp(- 2d0 * eta * dem1 * dem1)

            int3(p,k,n) = int3(p,k,n) + num*(reg1/dem1)

          end do
        end do

        do a=nO+1,nOrb-nR
          do b=nO+1,nOrb-nR

            num  = - 1.5d0 * ERI(p,k,a,b) * ee_trip_rho(a,b,n)
            dem1 = eQP(a) + eQP(b) - ee_trip_Om(n)
            reg1 = 1d0 - exp(- 2d0 * eta * dem1 * dem1)

            int3(p,k,n) = int3(p,k,n) + num*(reg1/dem1)

          end do
        end do

      end do
    end do
  end do

  do q=nC+1,nOrb-nR
    do p=nC+1,nOrb-nR
      do k=nC+1,nO
        do n=1,nVVt

          num  = ee_trip_rho(q,k,n)
          dem2_p = eQP(p) + eQP(k) - ee_trip_Om(n)
          dem2_q = eQP(q) + eQP(k) - ee_trip_Om(n)
          dem2_pq = dem2_p * dem2_p + dem2_q * dem2_q
          reg2 = 1d0 - exp(- eta * dem2_pq)

          SigC(p,q) = SigC(p,q) + num * int3(p,k,n) * reg2 * (dem2_p + dem2_q)/dem2_pq

        end do
      end do
    end do
  end do

  deallocate(int3)

!-------------------------------------!
! Intermediates for triplet pp part 2 !
!-------------------------------------!

  allocate(int4(nOrb,nOrb,nOrb,nOrb))
  int4(:,:,:,:) = 0d0

  do q=nC+1,nOrb-nR
    do i=nC+1,nO
      do j=nC+1,nO
        do a=nO+1,nOrb-nR

          do n=1,nVVt

            num  = ee_trip_rho(i,j,n) * ee_trip_rho(q,a,n)
            dem1 = ee_trip_Om(n) - eQP(i) - eQP(j)
            reg1 = 1d0 - exp(- 2d0 * eta * dem1 * dem1)

            int4(q,i,j,a) = int4(q,i,j,a) + num*(reg1/dem1)

          end do

          do n=1,nOOt

            num  = hh_trip_rho(i,j,n) * hh_trip_rho(q,a,n)
            dem1 = eQP(i) + eQP(j) - hh_trip_Om(n)
            reg1 = 1d0 - exp(- 2d0*eta * dem1 * dem1)

            int4(q,i,j,a) = int4(q,i,j,a) + num*(reg1/dem1)

          end do

        end do
      end do
    end do
  end do

  do q=nC+1,nOrb-nR
    do p=nC+1,nOrb-nR
      do i=nC+1,nO
        do j=nC+1,nO
          do a=nO+1,nOrb-nR

            num  = - 1.5d0*ERI(p,a,i,j)
            dem2_p = eQP(p) + eQP(a) - eQP(i) - eQP(j)
            dem2_q = eQP(q) + eQP(a) - eQP(i) - eQP(j)
            dem2_pq = dem2_p * dem2_p + dem2_q * dem2_q
            reg2 = 1d0 - exp(- eta * dem2_pq)

            SigC(p,q) = SigC(p,q) + num * int4(q,i,j,a) * reg2 * (dem2_p + dem2_q)/dem2_pq

          end do
        end do
      end do
    end do
  end do

  deallocate(int4)

!-------------------------------------!
! Intermediates for triplet pp part 3 !
!-------------------------------------!

  allocate(int3(nOrb,nOrb,nOOt))
  int3(:,:,:) = 0d0

  do p=nC+1,nOrb-nR
    do c=nO+1,nOrb-nR
      do n=1,nOOt

        do i=nC+1,nO
          do j=nC+1,nO

            num  = 1.5d0*ERI(p,c,i,j) * hh_trip_rho(i,j,n)
            dem1 = eQP(i) + eQP(j) - hh_trip_Om(n)
            reg1 = 1d0 - exp(- 2d0*eta * dem1 * dem1)

            int3(p,c,n) = int3(p,c,n) + num*(reg1/dem1)

          end do
        end do

        do a=nO+1,nOrb-nR
          do b=nO+1,nOrb-nR

            num  = 1.5d0*ERI(p,c,a,b) * hh_trip_rho(a,b,n)
            dem1 = hh_trip_Om(n) - eQP(a) - eQP(b)
            reg1 = 1d0 - exp(- 2d0*eta * dem1 * dem1)

            int3(p,c,n) = int3(p,c,n) + num*(reg1/dem1)


          end do
        end do

      end do
    end do
  end do

  do q=nC+1,nOrb-nR
    do p=nC+1,nOrb-nR
      do c=nO+1,nOrb-nR
        do n=1,nOOt

          num  = hh_trip_rho(q,c,n)
          dem2_p = eQP(p) + eQP(c) - hh_trip_Om(n)
          dem2_q = eQP(q) + eQP(c) - hh_trip_Om(n)
          dem2_pq = dem2_p * dem2_p + dem2_q * dem2_q
          reg2 = 1d0 - exp(- eta * dem2_pq)

          SigC(p,q) = SigC(p,q) + num * int3(p,c,n) * reg2 * (dem2_p + dem2_q)/dem2_pq

        end do
      end do
    end do
  end do

  deallocate(int3)

!-------------------------------------!
! Intermediates for triplet pp part 4 !
!-------------------------------------!

  allocate(int4(nOrb,nOrb,nOrb,nOrb))
  int4(:,:,:,:) = 0d0

  do q=nC+1,nOrb-nR
    do i=nC+1,nO
      do a=nO+1,nOrb-nR
        do b=nO+1,nOrb-nR

          do n=1,nOOt

            num  = hh_trip_rho(a,b,n) * hh_trip_rho(q,i,n)
            dem1 = hh_trip_Om(n) - eQP(a) - eQP(b)
            reg1 = 1d0 - exp(- 2d0 * eta * dem1 * dem1)

            int4(q,i,a,b) = int4(q,i,a,b) + num*(reg1/dem1)

          end do

          do n=1,nVVt

            num  = ee_trip_rho(a,b,n) * ee_trip_rho(q,i,n)
            dem1 = eQP(a) + eQP(b) - ee_trip_Om(n)
            reg1 = 1d0 - exp(- 2d0 * eta * dem1 * dem1)

            int4(q,i,a,b) = int4(q,i,a,b) + num*(reg1/dem1)

          end do

        end do
      end do
    end do
  end do

  do q=nC+1,nOrb-nR
    do p=nC+1,nOrb-nR
      do i=nC+1,nO
        do a=nO+1,nOrb-nR
          do b=nO+1,nOrb-nR

            num  = 1.5d0*ERI(p,i,a,b)
            dem2_p = eQP(p) + eQP(i) - eQP(a) - eQP(b)
            dem2_q = eQP(q) + eQP(i) - eQP(a) - eQP(b)
            dem2_pq = dem2_p * dem2_p + dem2_q * dem2_q
            reg2 = 1d0 - exp(- eta * dem2_pq)

            SigC(p,q) = SigC(p,q) + num * int4(q,i,a,b) * reg2 * (dem2_p + dem2_q)/dem2_pq

          end do
        end do
      end do
    end do
  end do

  deallocate(int4)

  call wall_time(end_t)
  t = end_t - start_t

  write(*,'(1X,A50,1X,F9.3,A8)') 'Wall time for triplet pp self-energy =',t,' seconds'
  write(*,*)
 
!-----------------------------!
!   Renormalization factor    !
!-----------------------------!

  Z(:) = 1d0

end subroutine 
