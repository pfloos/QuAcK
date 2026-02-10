subroutine R_Parquet_self_energy_diag_psd(eta,nOrb,nC,nO,nV,nR,nS,nOOs,nVVs,nOOt,nVVt,eQP,ERI, &
                                          eh_sing_rho,eh_sing_Om,eh_trip_rho,eh_trip_Om, &
                                          ee_sing_rho,ee_sing_Om,ee_trip_rho,ee_trip_Om, &
                                          hh_sing_rho,hh_sing_Om,hh_trip_rho,hh_trip_Om, &
                                          EcGM,SigC,Z)

! Compute correlation part of the parquet self-energy and make it positive semi-definite
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
  double precision              :: num,dem,reg
  double precision              :: start_t,end_t,t

  double precision,allocatable  :: int3(:,:,:)
  double precision,allocatable  :: int4(:,:,:,:)
  
! Output variables
  double precision,intent(out)  :: EcGM
  double precision,intent(out)  :: SigC(nOrb)
  double precision,intent(out)  :: Z(nOrb)

! Initialize 

  SigC(:) = 0d0
  Z(:) = 0d0

!-----------------------------------!
! 2nd-order part of the self-energy !
!-----------------------------------!
  call wall_time(start_t)

  do p=nC+1,nOrb-nR
    ! 2h1p sum
    do i=nC+1,nO
      do j=nC+1,nO
        do a=nO+1,nOrb-nR

          num = ERI(p,a,j,i) * ( 2d0*ERI(j,i,p,a) - ERI(j,i,a,p) )
          dem = eQP(p) + eQP(a) - eQP(i) - eQP(j)
          reg = 1d0 - exp(- 2d0 * eta * dem * dem)/dem
          
          SigC(p) = SigC(p) + num * reg
          Z(p) = Z(p) - num * reg / dem
          
        end do ! a
      end do ! j
    end do ! i
    ! 2p1h sum
    do i=nC+1,nO
      do a=nO+1,nOrb-nR
        do b=nO+1,nOrb-nR

          num  = ERI(p,i,b,a) * ( 2d0*ERI(b,a,p,i) - ERI(b,a,i,p) )
          dem = eQP(p) + eQP(i) - eQP(a) - eQP(b)
          reg = 1d0 - exp(- 2d0 * eta * dem * dem)/dem
          
          SigC(p) = SigC(p) + num * reg
          Z(p) = Z(p) - num * reg / dem
               
        end do ! b
      end do ! a
    end do ! i
      
  end do ! p

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

  !$OMP PARALLEL DEFAULT(NONE)    &                                       
  !$OMP PRIVATE(p,i,a,j,n,num,dem,reg) & 
  !$OMP SHARED(nC,nO,nOrb,nR,nS,eta,ERI,eQP,eh_sing_rho,eh_sing_Om,int3)  
  !$OMP DO
  do n=1,nS
     do j=nC+1,nO
        do p=nC+1,nOrb-nR
           
           do i=nC+1,nO
              do a=nO+1,nOrb-nR
           
                 num = (0.5d0*ERI(p,a,i,j) - ERI(p,a,j,i)) * eh_sing_rho(i,a,n)
                 dem = eQP(a) - eQP(i) - eh_sing_Om(n)
                 reg = 1d0 - exp(- 2d0 * eta * dem * dem)/dem
                 
                 int3(p,j,n) = int3(p,j,n) + num*reg
                 
                 num = (0.5d0*ERI(p,i,a,j) - ERI(p,i,j,a)) * eh_sing_rho(a,i,n)
                 dem = eQP(a) - eQP(i) + eh_sing_Om(n) 
                 reg = 1d0 - exp(- 2d0 * eta * dem * dem)/dem
                 
                 int3(p,j,n) = int3(p,j,n) + num*reg
                 
              end do
           end do
           
        end do
     end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL 

  !$OMP PARALLEL DEFAULT(NONE)    &                                       
  !$OMP PRIVATE(p,j,n,num,dem,reg) & 
  !$OMP SHARED(nC,nO,nOrb,nR,nS,eta,eQP,eh_sing_rho,eh_sing_Om,int3,SigC)  
  !$OMP DO
  do p=nC+1,nOrb-nR
    do j=nC+1,nO
      do n=1,nS

        num = eh_sing_rho(p,j,n) + int3(p,j,n)
        dem = eQP(p) - eQP(j) + eh_sing_Om(n)
        reg = 1d0 - exp(- 2d0*eta * dem * dem)
           
        SigC(p) = SigC(p) + num * num * reg
        Z(p) = Z(p) - num * num * reg / dem

      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL 

  deallocate(int3)

!-------------------------------------!
! Intermediates for singlet eh part 2 !
!-------------------------------------!

  allocate(int4(nOrb,nOrb,nOrb,nOrb))
  int4(:,:,:,:) = 0d0

  !$OMP PARALLEL DEFAULT(NONE)    &                                       
  !$OMP PRIVATE(p,i,a,j,n,num,dem,reg) & 
  !$OMP SHARED(nC,nO,nOrb,nR,nS,eta,eQP,eh_sing_rho,eh_sing_Om,int4)  
  !$OMP DO
  do a=nO+1,nOrb-nR
    do j=nC+1,nO
      do i=nC+1,nO
        do p=nC+1,nOrb-nR

          do n=1,nS

            num = - eh_sing_rho(i,a,n) * eh_sing_rho(p,j,n)
            dem = eQP(a) - eQP(i) - eh_sing_Om(n)
            reg = 1d0 - exp(- 2d0 * eta * dem * dem)/dem

            int4(p,i,j,a) = int4(p,i,j,a) + num*reg

            num  = eh_sing_rho(a,i,n) * eh_sing_rho(j,p,n)
            dem = eQP(a) - eQP(i) + eh_sing_Om(n)
            reg = 1d0 - exp(- 2d0 * eta * dem * dem)/dem

            int4(p,i,j,a) = int4(p,i,j,a) + num*reg

         end do
         

        end do 
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL 

  !$OMP PARALLEL DEFAULT(NONE)    &                                       
  !$OMP PRIVATE(p,i,j,a,num,dem,reg) & 
  !$OMP SHARED(nC,nO,nOrb,nR,nS,eta,eQP,ERI,int4,SigC)  
  !$OMP DO
  do p=nC+1,nOrb-nR
    do i=nC+1,nO
      do j=nC+1,nO
        do a=nO+1,nOrb-nR

          num = 0.5d0*ERI(p,a,i,j) - ERI(p,a,j,i) + int4(p,i,j,a)
          dem = eQP(p) - eQP(i) - eQP(j) + eQP(a)
          reg = 1d0 - exp(- 2d0 * eta * dem * dem)/dem
             
          SigC(p) = SigC(p) + num * num * reg
          Z(p) = Z(p) - num * num * reg / dem

        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL 

  deallocate(int4)

!-------------------------------------!
! Intermediates for singlet eh part 3 !
!-------------------------------------!

  allocate(int3(nOrb,nOrb,nS))
  int3(:,:,:) = 0d0

  !$OMP PARALLEL DEFAULT(NONE)    &                                       
  !$OMP PRIVATE(p,i,a,b,n,num,dem,reg) & 
  !$OMP SHARED(nC,nO,nOrb,nR,nS,eta,ERI,eQP,eh_sing_rho,eh_sing_Om,int3)  
  !$OMP DO
  do n=1,nS
    do b=nO+1,nOrb-nR
      do p=nC+1,nOrb-nR

        do i=nC+1,nO
          do a=nO+1,nOrb-nR

            num  = (0.5d0*ERI(p,i,a,b) - ERI(p,i,b,a)) * eh_sing_rho(i,a,n)
            dem = eQP(a) - eQP(i) - eh_sing_Om(n)
            reg = 1d0 - exp(- 2d0 * eta * dem * dem)/dem

            int3(p,b,n) = int3(p,b,n) + num*reg

            num  = (0.5d0*ERI(p,a,i,b) - ERI(p,a,b,i)) * eh_sing_rho(a,i,n) 
            dem = eQP(a) - eQP(i) + eh_sing_Om(n)
            reg = 1d0 - exp(- 2d0 * eta * dem * dem)/dem

            int3(p,b,n) = int3(p,b,n) + num*reg

          end do
        end do

      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL 

  !$OMP PARALLEL DEFAULT(NONE)    &                                       
  !$OMP PRIVATE(p,b,n,num,dem,reg) & 
  !$OMP SHARED(nC,nO,nOrb,nR,nS,eta,eQP,eh_sing_rho,eh_sing_Om,int3,SigC)  
  !$OMP DO
  do p=nC+1,nOrb-nR
    do b=nO+1,nOrb-nR
      do n=1,nS

        num = eh_sing_rho(b,p,n) + int3(p,b,n)
        dem = eQP(p) - eQP(b) - eh_sing_Om(n)
        reg = 1d0 - exp(- 2d0 * eta * dem * dem)/dem

        SigC(p) = SigC(p) + num * num * reg
        Z(p) = Z(p) - num * num * reg / dem

      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL 

  deallocate(int3)

!-------------------------------------!
! Intermediates for singlet eh part 4 !
!-------------------------------------!

  allocate(int4(nOrb,nOrb,nOrb,nOrb))
  int4(:,:,:,:) = 0d0

  !$OMP PARALLEL DEFAULT(NONE)    &                                       
  !$OMP PRIVATE(p,i,a,j,n,num,dem,reg) & 
  !$OMP SHARED(nC,nO,nOrb,nR,nS,eta,eQP,eh_sing_rho,eh_sing_Om,int4)  
  !$OMP DO
  do b=nO+1,nOrb-nR
    do a=nO+1,nOrb-nR
      do i=nC+1,nO
        do p=nC+1,nOrb-nR

          do n=1,nS

            num = - eh_sing_rho(i,a,n) * eh_sing_rho(b,p,n)
            dem = eQP(a) - eQP(i) - eh_sing_Om(n)
            reg = 1d0 - exp(- 2d0 * eta * dem * dem)/dem

            int4(p,i,a,b) = int4(p,i,a,b) + num*reg

            num = eh_sing_rho(a,i,n) * eh_sing_rho(p,b,n)
            dem = eQP(a) - eQP(i) + eh_sing_Om(n)
            reg = 1d0 - exp(- 2d0 * eta * dem * dem)/dem

            int4(p,i,a,b) = int4(p,i,a,b) + num*reg

          end do

        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL 

  !$OMP PARALLEL DEFAULT(NONE)    &                                       
  !$OMP PRIVATE(p,i,a,b,num,dem,reg) & 
  !$OMP SHARED(nC,nO,nOrb,nR,nS,eta,eQP,ERI,int4,SigC)  
  !$OMP DO
  do p=nC+1,nOrb-nR
    do i=nC+1,nO
      do a=nO+1,nOrb-nR
        do b=nO+1,nOrb-nR

          num  = 0.5d0*ERI(p,i,a,b) - ERI(p,i,b,a) + int4(p,i,a,b)
          dem = eQP(p) + eQP(i) - eQP(a) - eQP(b)
          reg = 1d0 - exp(- 2d0 * eta * dem * dem)/dem
             
          SigC(p) = SigC(p) + num * num * reg
          Z(p) = Z(p) - num * num * reg / dem

        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL 

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
  
  !$OMP PARALLEL DEFAULT(NONE)    &                                       
  !$OMP PRIVATE(p,i,a,j,n,num,dem,reg) & 
  !$OMP SHARED(nC,nO,nOrb,nR,nS,eta,ERI,eQP,eh_trip_rho,eh_trip_Om,int3)  
  !$OMP DO
  do n=1,nS
    do j=nC+1,nO
      do p=nC+1,nOrb-nR

        do i=nC+1,nO
          do a=nO+1,nOrb-nR

            num = 1.5d0*ERI(p,a,i,j) * eh_trip_rho(i,a,n)
            dem = eQP(a) - eQP(i) - eh_trip_Om(n)
            reg = 1d0 - exp(- 2d0 * eta * dem * dem)/dem

            int3(p,j,n) = int3(p,j,n) + num*reg

            num = 1.5d0*ERI(p,i,a,j) * eh_trip_rho(a,i,n)
            dem = eQP(a) - eQP(i) + eh_trip_Om(n) 
            reg = 1d0 - exp(- 2d0 * eta * dem * dem)/dem

            int3(p,j,n) = int3(p,j,n) + num*reg

          end do
        end do

      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL 

  !$OMP PARALLEL DEFAULT(NONE)    &                                       
  !$OMP PRIVATE(p,j,n,num,dem,reg) & 
  !$OMP SHARED(nC,nO,nOrb,nR,nS,eta,eQP,eh_trip_rho,eh_trip_Om,int3,SigC)  
  !$OMP DO
  do p=nC+1,nOrb-nR
    do j=nC+1,nO
      do n=1,nS

        num = eh_trip_rho(p,j,n) + int3(p,j,n)
        dem = eQP(p) - eQP(j) + eh_trip_Om(n)
        reg = 1d0 - exp(- 2d0 * eta * dem * dem)/dem
           
        SigC(p) = SigC(p) + num * num * reg
        Z(p) = Z(p) - num * num * reg / dem

      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL 

  deallocate(int3)

!-------------------------------------!
! Intermediates for triplet eh part 2 !
!-------------------------------------!

  allocate(int4(nOrb,nOrb,nOrb,nOrb))
  int4(:,:,:,:) = 0d0

  !$OMP PARALLEL DEFAULT(NONE)    &                                       
  !$OMP PRIVATE(p,i,a,j,n,num,dem,reg) & 
  !$OMP SHARED(nC,nO,nOrb,nR,nS,eta,eQP,eh_trip_rho,eh_trip_Om,int4)  
  !$OMP DO
  do a=nO+1,nOrb-nR
    do j=nC+1,nO
      do i=nC+1,nO
        do p=nC+1,nOrb-nR

          do n=1,nS

            num = - eh_trip_rho(i,a,n) * eh_trip_rho(p,j,n)
            dem = eQP(a) - eQP(i) - eh_trip_Om(n)
            reg = 1d0 - exp(- 2d0 * eta * dem * dem)/dem

            int4(p,i,j,a) = int4(p,i,j,a) + num*reg

            num = eh_trip_rho(a,i,n) * eh_trip_rho(j,p,n)
            dem = eQP(a) - eQP(i) + eh_trip_Om(n)
            reg = 1d0 - exp(- 2d0 * eta * dem * dem)/dem

            int4(p,i,j,a) = int4(p,i,j,a) + num*reg

          end do

        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL 

  !$OMP PARALLEL DEFAULT(NONE)    &                                       
  !$OMP PRIVATE(p,i,j,a,num,dem,reg) & 
  !$OMP SHARED(nC,nO,nOrb,nR,nS,eta,eQP,ERI,int4,SigC)  
  !$OMP DO
  do p=nC+1,nOrb-nR
    do i=nC+1,nO
      do j=nC+1,nO
        do a=nO+1,nOrb-nR

          num = 1.5d0*ERI(p,a,i,j) + int4(p,i,j,a)
          dem = eQP(p) - eQP(i) - eQP(j) + eQP(a)
          reg = 1d0 - exp(- 2d0 * eta * dem * dem)/dem
             
          SigC(p) = SigC(p) + num * num * reg
          Z(p) = Z(p) - num * num * reg / dem

        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL 

  deallocate(int4)

!-------------------------------------!
! Intermediates for triplet eh part 3 !
!-------------------------------------!

  allocate(int3(nOrb,nOrb,nS))
  int3(:,:,:) = 0d0

  !$OMP PARALLEL DEFAULT(NONE)    &                                       
  !$OMP PRIVATE(p,i,a,b,n,num,dem,reg) & 
  !$OMP SHARED(nC,nO,nOrb,nR,nS,eta,ERI,eQP,eh_trip_rho,eh_trip_Om,int3)  
  !$OMP DO
  do n=1,nS
    do b=nO+1,nOrb-nR
      do p=nC+1,nOrb-nR

        do i=nC+1,nO
          do a=nO+1,nOrb-nR

            num  = 1.5d0*ERI(p,i,a,b) * eh_trip_rho(i,a,n)
            dem = eQP(a) - eQP(i) - eh_trip_Om(n)
            reg = 1d0 - exp(- 2d0 * eta * dem * dem)/dem

            int3(p,b,n) = int3(p,b,n) + num*reg

            num  = 1.5d0*ERI(p,a,i,b) * eh_trip_rho(a,i,n) 
            dem = eQP(a) - eQP(i) + eh_trip_Om(n)
            reg = 1d0 - exp(- 2d0 * eta * dem * dem)/dem

            int3(p,b,n) = int3(p,b,n) + num*reg

          end do
        end do

      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL 

  !$OMP PARALLEL DEFAULT(NONE)    &                                       
  !$OMP PRIVATE(p,b,n,num,dem,reg) & 
  !$OMP SHARED(nC,nO,nOrb,nR,nS,eta,eQP,eh_trip_rho,eh_trip_Om,int3,SigC)  
  !$OMP DO
  do p=nC+1,nOrb-nR
    do b=nO+1,nOrb-nR
      do n=1,nS

        num = eh_trip_rho(b,p,n) + int3(p,b,n)
        dem = eQP(p) - eQP(b) - eh_trip_Om(n)
        reg = 1d0 - exp(- 2d0 * eta * dem * dem)/dem

        SigC(p) = SigC(p) + num * num * reg
        Z(p) = Z(p) - num * num * reg / dem

      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL 

  deallocate(int3)

!-------------------------------------!
! Intermediates for triplet eh part 4 !
!-------------------------------------!

  allocate(int4(nOrb,nOrb,nOrb,nOrb))
  int4(:,:,:,:) = 0d0

  !$OMP PARALLEL DEFAULT(NONE)    &                                       
  !$OMP PRIVATE(p,i,a,j,n,num,dem,reg) & 
  !$OMP SHARED(nC,nO,nOrb,nR,nS,eta,eQP,eh_trip_rho,eh_trip_Om,int4)  
  !$OMP DO
  do b=nO+1,nOrb-nR
    do a=nO+1,nOrb-nR
      do i=nC+1,nO
        do p=nC+1,nOrb-nR

          do n=1,nS

            num  = - eh_trip_rho(i,a,n) * eh_trip_rho(b,p,n)
            dem = eQP(a) - eQP(i) - eh_trip_Om(n)
            reg = 1d0 - exp(- 2d0 * eta * dem * dem)/dem

            int4(p,i,a,b) = int4(p,i,a,b) + num*reg

            num  = eh_trip_rho(a,i,n) * eh_trip_rho(p,b,n)
            dem = eQP(a) - eQP(i) + eh_trip_Om(n)
            reg = 1d0 - exp(- 2d0 * eta * dem * dem)/dem

            int4(p,i,a,b) = int4(p,i,a,b) + num*reg

          end do

        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL 

  !$OMP PARALLEL DEFAULT(NONE)    &                                       
  !$OMP PRIVATE(p,i,a,b,num,dem,reg) & 
  !$OMP SHARED(nC,nO,nOrb,nR,nS,eta,eQP,ERI,int4,SigC)  
  !$OMP DO
  do p=nC+1,nOrb-nR
    do i=nC+1,nO
      do a=nO+1,nOrb-nR
        do b=nO+1,nOrb-nR

          num = 1.5d0*ERI(p,i,a,b) + int4(p,i,a,b)
          dem = eQP(p) + eQP(i) - eQP(a) - eQP(b)
          reg = 1d0 - exp(- 2d0 * eta * dem * dem)/dem
             
          SigC(p) = SigC(p) + num * num * reg
          Z(p) = Z(p) - num * num * reg / dem

        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL 

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

  !$OMP PARALLEL DEFAULT(NONE)    &                                       
  !$OMP PRIVATE(p,i,j,k,a,b,n,num,dem,reg) & 
  !$OMP SHARED(nC,nO,nOrb,nR,nVVs,eta,eQP,ERI,ee_sing_rho,ee_sing_Om,int3)  
  !$OMP DO
  do n=1,nVVs
    do k=nC+1,nO
      do p=nC+1,nOrb-nR

        do i=nC+1,nO
          do j=nC+1,nO

            num  = - 0.5d0 * ERI(p,k,i,j) * ee_sing_rho(i,j,n)
            dem = ee_sing_Om(n) - eQP(i) - eQP(j)
            reg = 1d0 - exp(- 2d0 * eta * dem * dem)/dem

            int3(p,k,n) = int3(p,k,n) + num*reg

          end do
        end do

        do a=nO+1,nOrb-nR
          do b=nO+1,nOrb-nR

            num  = - 0.5d0 * ERI(p,k,a,b) * ee_sing_rho(a,b,n)
            dem = eQP(a) + eQP(b) - ee_sing_Om(n)
            reg = 1d0 - exp(- 2d0 * eta * dem * dem)/dem

            int3(p,k,n) = int3(p,k,n) + num*reg

          end do
        end do

      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL 

  !$OMP PARALLEL DEFAULT(NONE)    &                                       
  !$OMP PRIVATE(p,k,n,num,dem,reg) & 
  !$OMP SHARED(nC,nO,nOrb,nR,nVVs,eta,eQP,ee_sing_rho,ee_sing_Om,int3,SigC)  
  !$OMP DO
  do p=nC+1,nOrb-nR
    do k=nC+1,nO
      do n=1,nVVs

        num = ee_sing_rho(p,k,n) + int3(p,k,n)
        dem = eQP(p) + eQP(k) - ee_sing_Om(n)
        reg = 1d0 - exp(- 2d0 * eta * dem * dem)/dem

        SigC(p) = SigC(p) + num * num * reg
        Z(p) = Z(p) - num * num * reg / dem

      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL 

  deallocate(int3)

!-------------------------------------!
! Intermediates for singlet pp part 2 !
!-------------------------------------!

  allocate(int4(nOrb,nOrb,nOrb,nOrb))
  int4(:,:,:,:) = 0d0

  !$OMP PARALLEL DEFAULT(NONE)    &                                       
  !$OMP PRIVATE(a,i,j,p,n,num,dem,reg) & 
  !$OMP SHARED(nC,nO,nOrb,nR,nVVs,nOOs,eta,eQP,ee_sing_rho,ee_sing_Om,hh_sing_rho,hh_sing_Om,int4)  
  !$OMP DO
  do a=nO+1,nOrb-nR
    do j=nC+1,nO
      do i=nC+1,nO
        do p=nC+1,nOrb-nR

          do n=1,nVVs

            num  = ee_sing_rho(i,j,n) * ee_sing_rho(p,a,n)
            dem = ee_sing_Om(n) - eQP(i) - eQP(j)
            reg = 1d0 - exp(- 2d0 * eta * dem * dem)/dem

            int4(p,i,j,a) = int4(p,i,j,a) + num*reg

          end do

          do n=1,nOOs

            num  = hh_sing_rho(i,j,n) * hh_sing_rho(p,a,n)
            dem = eQP(i) + eQP(j) - hh_sing_Om(n)
            reg = 1d0 - exp(- 2d0*eta * dem * dem)/dem

            int4(p,i,j,a) = int4(p,i,j,a) + num*reg

          end do

        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL 
 
  !$OMP PARALLEL DEFAULT(NONE)    &                                       
  !$OMP PRIVATE(p,i,j,a,num,dem,reg) & 
  !$OMP SHARED(nC,nO,nOrb,nR,eta,eQP,ERI,int4,SigC)  
  !$OMP DO
  do p=nC+1,nOrb-nR
    do i=nC+1,nO
      do j=nC+1,nO
        do a=nO+1,nOrb-nR

          num = - 0.5d0*ERI(p,a,i,j) + int4(p,i,j,a)
          dem = eQP(p) + eQP(a) - eQP(i) - eQP(j)
          reg = 1d0 - exp(- 2d0*eta * dem * dem)/dem

          SigC(p) = SigC(p) + num * num * reg
          Z(p) = Z(p) - num * num * reg / dem

        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL 

  deallocate(int4)

!-------------------------------------!
! Intermediates for singlet pp part 3 !
!-------------------------------------!

  allocate(int3(nOrb,nOrb,nOOs))
  int3(:,:,:) = 0d0

  !$OMP PARALLEL DEFAULT(NONE)    &                                       
  !$OMP PRIVATE(p,i,j,a,b,c,n,num,dem,reg) & 
  !$OMP SHARED(nC,nO,nOrb,nR,nOOs,eta,eQP,ERI,hh_sing_rho,hh_sing_Om,int3)  
  !$OMP DO
  do n=1,nOOs
    do c=nO+1,nOrb-nR
      do p=nC+1,nOrb-nR

        do i=nC+1,nO
          do j=nC+1,nO

            num  = 0.5d0*ERI(p,c,i,j) * hh_sing_rho(i,j,n)
            dem = eQP(i) + eQP(j) - hh_sing_Om(n)
            reg = 1d0 - exp(- 2d0*eta * dem * dem)/dem

            int3(p,c,n) = int3(p,c,n) + num*reg

          end do
        end do

        do a=nO+1,nOrb-nR
          do b=nO+1,nOrb-nR

            num  = 0.5d0*ERI(p,c,a,b) * hh_sing_rho(a,b,n)
            dem = hh_sing_Om(n) - eQP(a) - eQP(b)
            reg = 1d0 - exp(- 2d0*eta * dem * dem)/dem

            int3(p,c,n) = int3(p,c,n) + num*reg


          end do
        end do

      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL 

  !$OMP PARALLEL DEFAULT(NONE)    &                                       
  !$OMP PRIVATE(p,c,n,num,dem,reg) & 
  !$OMP SHARED(nC,nO,nOrb,nR,nOOs,eta,eQP,hh_sing_rho,hh_sing_Om,int3,SigC)  
  !$OMP DO
  do p=nC+1,nOrb-nR
    do c=nO+1,nOrb-nR
      do n=1,nOOs

        num  = hh_sing_rho(p,c,n) + int3(p,c,n)
        dem = eQP(p) + eQP(c) - hh_sing_Om(n)
        reg = 1d0 - exp(- 2d0*eta * dem * dem)/dem

        SigC(p) = SigC(p) + num * num * reg
        Z(p) = Z(p) - num * num * reg / dem

      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL 

  deallocate(int3)

!-------------------------------------!
! Intermediates for singlet pp part 4 !
!-------------------------------------!

  allocate(int4(nOrb,nOrb,nOrb,nOrb))
  int4(:,:,:,:) = 0d0

  !$OMP PARALLEL DEFAULT(NONE)    &                                       
  !$OMP PRIVATE(b,a,i,p,n,num,dem,reg) & 
  !$OMP SHARED(nC,nO,nOrb,nR,nVVs,nOOs,eta,eQP,ee_sing_rho,ee_sing_Om,hh_sing_rho,hh_sing_Om,int4)  
  !$OMP DO
  do b=nO+1,nOrb-nR
    do a=nO+1,nOrb-nR
      do i=nC+1,nO
        do p=nC+1,nOrb-nR

          do n=1,nOOs

            num  = hh_sing_rho(a,b,n) * hh_sing_rho(p,i,n)
            dem = hh_sing_Om(n) - eQP(a) - eQP(b)
            reg = 1d0 - exp(- 2d0 * eta * dem * dem)/dem

            int4(p,i,a,b) = int4(p,i,a,b) + num*reg

          end do

          do n=1,nVVs

            num  = ee_sing_rho(a,b,n) * ee_sing_rho(p,i,n)
            dem = eQP(a) + eQP(b) - ee_sing_Om(n)
            reg = 1d0 - exp(- 2d0 * eta * dem * dem)/dem

            int4(p,i,a,b) = int4(p,i,a,b) + num*reg

          end do

        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL 

  !$OMP PARALLEL DEFAULT(NONE)    &                                       
  !$OMP PRIVATE(p,i,a,b,num,dem,reg) & 
  !$OMP SHARED(nC,nO,nOrb,nR,eta,eQP,ERI,int4,SigC)  
  !$OMP DO
  do p=nC+1,nOrb-nR
    do i=nC+1,nO
      do a=nO+1,nOrb-nR
        do b=nO+1,nOrb-nR

          num = 0.5d0*ERI(p,i,a,b) + int4(p,i,a,b)
          dem = eQP(p) + eQP(i) - eQP(a) - eQP(b)
          reg = 1d0 - exp(- 2d0 * eta * dem * dem)/dem

          SigC(p) = SigC(p) + num * num * reg
          Z(p) = Z(p) - num * num * reg / dem

        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL 

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

  !$OMP PARALLEL DEFAULT(NONE)    &                                       
  !$OMP PRIVATE(p,i,j,k,a,b,n,num,dem,reg) & 
  !$OMP SHARED(nC,nO,nOrb,nR,nVVt,eta,eQP,ERI,ee_trip_rho,ee_trip_Om,int3)  
  !$OMP DO
  do n=1,nVVt
    do k=nC+1,nO
      do p=nC+1,nOrb-nR

        do i=nC+1,nO
          do j=nC+1,nO

            num  = - 1.5d0 * ERI(p,k,i,j) * ee_trip_rho(i,j,n)
            dem = ee_trip_Om(n) - eQP(i) - eQP(j)
            reg = 1d0 - exp(- 2d0 * eta * dem * dem)/dem

            int3(p,k,n) = int3(p,k,n) + num*reg

          end do
        end do

        do a=nO+1,nOrb-nR
          do b=nO+1,nOrb-nR

            num  = - 1.5d0 * ERI(p,k,a,b) * ee_trip_rho(a,b,n)
            dem = eQP(a) + eQP(b) - ee_trip_Om(n)
            reg = 1d0 - exp(- 2d0 * eta * dem * dem)/dem

            int3(p,k,n) = int3(p,k,n) + num*reg

          end do
        end do

      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL 

  !$OMP PARALLEL DEFAULT(NONE)    &                                       
  !$OMP PRIVATE(p,k,n,num,dem,reg) & 
  !$OMP SHARED(nC,nO,nOrb,nR,nVVt,eta,eQP,ee_trip_rho,ee_trip_Om,int3,SigC)  
  !$OMP DO
  do p=nC+1,nOrb-nR
    do k=nC+1,nO
      do n=1,nVVt

        num = ee_trip_rho(p,k,n) + int3(p,k,n)
        dem = eQP(p) + eQP(k) - ee_trip_Om(n)
        reg = 1d0 - exp(- 2d0 * eta * dem * dem)/dem

        SigC(p) = SigC(p) + num * num * reg
        Z(p) = Z(p) - num * num * reg / dem

      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL 

  deallocate(int3)

!-------------------------------------!
! Intermediates for triplet pp part 2 !
!-------------------------------------!

  allocate(int4(nOrb,nOrb,nOrb,nOrb))
  int4(:,:,:,:) = 0d0

  !$OMP PARALLEL DEFAULT(NONE)    &                                       
  !$OMP PRIVATE(a,i,j,p,n,num,dem,reg) & 
  !$OMP SHARED(nC,nO,nOrb,nR,nVVt,nOOt,eta,eQP,ee_trip_rho,ee_trip_Om,hh_trip_rho,hh_trip_Om,int4)  
  !$OMP DO
  do a=nO+1,nOrb-nR
    do j=nC+1,nO
      do i=nC+1,nO
        do p=nC+1,nOrb-nR

          do n=1,nVVt

            num  = ee_trip_rho(i,j,n) * ee_trip_rho(p,a,n)
            dem = ee_trip_Om(n) - eQP(i) - eQP(j)
            reg = 1d0 - exp(- 2d0 * eta * dem * dem)/dem

            int4(p,i,j,a) = int4(p,i,j,a) + num*reg

          end do

          do n=1,nOOt

            num  = hh_trip_rho(i,j,n) * hh_trip_rho(p,a,n)
            dem = eQP(i) + eQP(j) - hh_trip_Om(n)
            reg = 1d0 - exp(- 2d0*eta * dem * dem)/dem

            int4(p,i,j,a) = int4(p,i,j,a) + num*reg

          end do

        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL 

  !$OMP PARALLEL DEFAULT(NONE)    &                                       
  !$OMP PRIVATE(p,i,j,a,num,dem,reg) & 
  !$OMP SHARED(nC,nO,nOrb,nR,eta,eQP,ERI,int4,SigC)  
  !$OMP DO
  do p=nC+1,nOrb-nR
    do i=nC+1,nO
      do j=nC+1,nO
        do a=nO+1,nOrb-nR

          num  = - 1.5d0*ERI(p,a,i,j) + int4(p,i,j,a)
          dem = eQP(p) + eQP(a) - eQP(i) - eQP(j)
          reg = 1d0 - exp(- 2d0 * eta * dem * dem)/dem

          SigC(p) = SigC(p) + num * num * reg
          Z(p) = Z(p) - num * num * reg / dem

        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL 

  deallocate(int4)

!-------------------------------------!
! Intermediates for triplet pp part 3 !
!-------------------------------------!

  allocate(int3(nOrb,nOrb,nOOt))
  int3(:,:,:) = 0d0
  
  !$OMP PARALLEL DEFAULT(NONE)    &                                       
  !$OMP PRIVATE(p,i,j,a,b,c,n,num,dem,reg) & 
  !$OMP SHARED(nC,nO,nOrb,nR,nOOt,eta,eQP,ERI,hh_trip_rho,hh_trip_Om,int3)  
  !$OMP DO
  do n=1,nOOt
    do c=nO+1,nOrb-nR
      do p=nC+1,nOrb-nR

        do i=nC+1,nO
          do j=nC+1,nO

            num  = 1.5d0*ERI(p,c,i,j) * hh_trip_rho(i,j,n)
            dem = eQP(i) + eQP(j) - hh_trip_Om(n)
            reg = 1d0 - exp(- 2d0*eta * dem * dem)/dem

            int3(p,c,n) = int3(p,c,n) + num*reg

          end do
        end do

        do a=nO+1,nOrb-nR
          do b=nO+1,nOrb-nR

            num  = 1.5d0*ERI(p,c,a,b) * hh_trip_rho(a,b,n)
            dem = hh_trip_Om(n) - eQP(a) - eQP(b)
            reg = 1d0 - exp(- 2d0*eta * dem * dem)/dem

            int3(p,c,n) = int3(p,c,n) + num*reg


          end do
        end do

      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL 

  !$OMP PARALLEL DEFAULT(NONE)    &                                       
  !$OMP PRIVATE(p,c,n,num,dem,reg) & 
  !$OMP SHARED(nC,nO,nOrb,nR,nOOt,eta,eQP,hh_trip_rho,hh_trip_Om,int3,SigC)  
  !$OMP DO
  do p=nC+1,nOrb-nR
    do c=nO+1,nOrb-nR
      do n=1,nOOt

        num  = hh_trip_rho(p,c,n) + int3(p,c,n)
        dem = eQP(p) + eQP(c) - hh_trip_Om(n)
        reg = 1d0 - exp(- 2d0*eta * dem * dem)/dem

        SigC(p) = SigC(p) + num * num * reg
        Z(p) = Z(p) - num * num * reg / dem

      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL 

  deallocate(int3)

!-------------------------------------!
! Intermediates for triplet pp part 4 !
!-------------------------------------!

  allocate(int4(nOrb,nOrb,nOrb,nOrb))
  int4(:,:,:,:) = 0d0

  !$OMP PARALLEL DEFAULT(NONE)    &                                       
  !$OMP PRIVATE(b,a,i,p,n,num,dem,reg) & 
  !$OMP SHARED(nC,nO,nOrb,nR,nVVt,nOOt,eta,eQP,ee_trip_rho,ee_trip_Om,hh_trip_rho,hh_trip_Om,int4)  
  !$OMP DO
  do b=nO+1,nOrb-nR
    do a=nO+1,nOrb-nR
      do i=nC+1,nO
        do p=nC+1,nOrb-nR

          do n=1,nOOt

            num  = hh_trip_rho(a,b,n) * hh_trip_rho(p,i,n)
            dem = hh_trip_Om(n) - eQP(a) - eQP(b)
            reg = 1d0 - exp(- 2d0 * eta * dem * dem)/dem

            int4(p,i,a,b) = int4(p,i,a,b) + num*reg

          end do

          do n=1,nVVt

            num  = ee_trip_rho(a,b,n) * ee_trip_rho(p,i,n)
            dem = eQP(a) + eQP(b) - ee_trip_Om(n)
            reg = 1d0 - exp(- 2d0 * eta * dem * dem)/dem

            int4(p,i,a,b) = int4(p,i,a,b) + num*reg

          end do

        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL 

  !$OMP PARALLEL DEFAULT(NONE)    &                                       
  !$OMP PRIVATE(p,i,a,b,num,dem,reg) & 
  !$OMP SHARED(nC,nO,nOrb,nR,eta,eQP,ERI,int4,SigC)  
  !$OMP DO
  do p=nC+1,nOrb-nR
    do i=nC+1,nO
      do a=nO+1,nOrb-nR
        do b=nO+1,nOrb-nR

          num  = 1.5d0*ERI(p,i,a,b) + int4(p,i,a,b)
          dem = eQP(p) + eQP(i) - eQP(a) - eQP(b)
          reg = 1d0 - exp(- 2d0 * eta * dem * dem)/dem

          SigC(p) = SigC(p) + num * num * reg
          Z(p) = Z(p) - num * num * reg / dem

        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL 

  deallocate(int4)

  call wall_time(end_t)
  t = end_t - start_t

  write(*,'(1X,A50,1X,F9.3,A8)') 'Wall time for triplet pp self-energy =',t,' seconds'
  write(*,*)
 
end subroutine 
