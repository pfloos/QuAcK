subroutine G_Parquet_self_energy_diag(eta,nOrb,nC,nO,nV,nR,nS,nOO,nVV,eQP,ERI,&
                                 eh_rho,eh_Om,ee_rho,ee_Om,hh_rho,hh_Om,EcGM,SigC,Z)

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
  integer                       :: p,n
  double precision              :: eps,dem1,dem2,reg,reg1,reg2
  double precision              :: num
  double precision              :: start_t,end_t,t

  logical                       :: print_self_energy
  double precision,allocatable  :: Sig_2(:)
  double precision,allocatable  :: Sig_eh(:,:)
  double precision,allocatable  :: Sig_pp(:,:)
  double precision,allocatable  :: Z_2(:)
  double precision,allocatable  :: Z_eh(:,:)
  double precision,allocatable  :: Z_pp(:,:)

! Output variables

  double precision,intent(out)  :: SigC(nOrb)
  double precision,intent(out)  :: Z(nOrb)
  double precision,intent(out)  :: EcGM

! Initialize 

  SigC(:) = 0d0
  Z(:)    = 0d0
  EcGM    = 0d0

! Memory allocation for self-energy decomposition

  allocate(Sig_2(nOrb))
  allocate(Sig_eh(nOrb,6))
  allocate(Sig_pp(nOrb,6))
 
  allocate(Z_2(nOrb))
  allocate(Z_eh(nOrb,6))
  allocate(Z_pp(nOrb,6))
 
!-----------------------------------!
! 2nd-order part of the self-energy !
!-----------------------------------!

  call wall_time(start_t)

  Sig_2(:) = 0d0
  Z_2(:)   = 0d0

  do p=nC+1,nOrb-nR

     ! 2h1p 
     do i=nC+1,nO
        do j=nC+1,nO
           do a=nO+1,nOrb-nR

              eps = eQP(p) + eQP(a) - eQP(i) - eQP(j)
              reg = (1d0 - exp(- 2d0 * eta * eps * eps))
              num = 0.5d0*(ERI(p,a,j,i) - ERI(p,a,i,j))**2

              Sig_2(p) = Sig_2(p) + num*reg/eps
              Z_2(p)   = Z_2(p)   - num*reg/eps**2

           end do
        end do
     end do

     ! 2p1h
     do i=nC+1,nO
        do a=nO+1,nOrb-nR
           do b=nO+1,nOrb-nR

              eps = eQP(p) + eQP(i) - eQP(a) - eQP(b)
              reg = (1d0 - exp(- 2d0 * eta * eps * eps))
              num = 0.5d0*(ERI(p,i,b,a) - ERI(p,i,a,b))**2

              Sig_2(p) = Sig_2(p) + num*reg/eps
              Z_2(p)   = Z_2(p)   - num*reg/eps**2

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

  Sig_eh(:,:) = 0d0
  Z_eh(:,:)   = 0d0

  !$OMP PARALLEL DEFAULT(NONE)    &
  !$OMP PRIVATE(p,i,a,j,b,n,num,dem1,dem2,reg1,reg2) &
  !$OMP SHARED(nC,nO,nOrb,nR,nS,eta,ERI,eQP,eh_rho,eh_Om,Sig_eh,Z_eh)
  !$OMP DO
  do p=nC+1,nOrb-nR
     
     do i=nC+1,nO
        do a=nO+1,nOrb-nR
           
           do n=1,nS

              do j=nC+1,nO

                 ! 2h1p(d) * 2h1p 
                 num  = (ERI(p,a,i,j) - ERI(p,a,j,i)) * eh_rho(i,a,n) * eh_rho(p,j,n)
                 ! num  = ERI(p,a,i,j) * ( eh_rho(i,a,n) * eh_rho(j,p,nS+n) - eh_rho(j,a,n) * eh_rho(i,p,nS+n) )
                 ! num  = 0.5d0 * (ERI(p,a,i,j) - ERI(p,a,j,i)) * &
                 !      ( eh_rho(i,a,n) * eh_rho(j,p,nS+n) - eh_rho(j,a,n) * eh_rho(i,p,nS+n) )
                 dem1 = eQP(p) - eQP(j) + eh_Om(n)
                 dem2 = eQP(p) - eQP(i) - eQP(j) + eQP(a)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 Sig_eh(p,1) = Sig_eh(p,1) + num * (reg1/dem1) * (reg2/dem2)
                 Z_eh(p,1)   = Z_eh(p,1)   - num * (reg1/dem1) * (reg2/dem2/dem2) &
                                           - num * (reg1/dem1/dem1) * (reg2/dem2) 
                 
                 ! 2h2p(d) * 2h1p(d) 
                 num  = (ERI(p,i,a,j) - ERI(p,i,j,a)) * eh_rho(a,i,n) * eh_rho(p,j,n)
                 ! num  = ERI(p,i,a,j) * ( eh_rho(a,i,n) * eh_rho(j,p,nS+n) -  eh_rho(j,i,n) * eh_rho(a,p,nS+n) )
                 ! num  = 0.5d0 * (ERI(p,i,a,j) - ERI(p,i,j,a)) * &
                 !      ( eh_rho(a,i,n) * eh_rho(j,p,nS+n) -  eh_rho(j,i,n) * eh_rho(a,p,nS+n) )
                 dem1 = eQP(a) - eQP(i) + eh_Om(n) 
                 dem2 = eQP(p) - eQP(j) + eh_Om(n)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))
 
                 Sig_eh(p,2) = Sig_eh(p,2) + num * (reg1/dem1) * (reg2/dem2)
                 Z_eh(p,2)   = Z_eh(p,2)   - num * (reg1/dem1) * (reg2/dem2/dem2)
                 
                 ! 2h2p(d) * 2h1p 
                 num  = (ERI(p,a,i,j) - ERI(p,a,j,i)) * eh_rho(a,i,n) * eh_rho(j,p,n) 
                 ! num  = ERI(p,a,i,j) * ( eh_rho(i,a,nS+n) * eh_rho(j,p,n) - eh_rho(j,a,nS+n) * eh_rho(i,p,n) ) 
                 ! num  =  0.5d0 * (ERI(p,a,i,j) - ERI(p,a,j,i)) * &
                 !      ( eh_rho(i,a,nS+n) * eh_rho(j,p,n) - eh_rho(j,a,nS+n) * eh_rho(i,p,n) ) 
                 dem1 = eQP(a) - eQP(i) + eh_Om(n) 
                 dem2 = eQP(p) - eQP(i) - eQP(j) + eQP(a)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 Sig_eh(p,3) = Sig_eh(p,3) + num * (reg1/dem1) * (reg2/dem2)
                 Z_eh(p,3)   = Z_eh(p,3)   - num * (reg1/dem1) * (reg2/dem2/dem2)

              end do ! j

              do b=nO+1,nOrb-nR

                 ! 2p1h(d) * 2p1h
                 num  = - (ERI(p,i,a,b) - ERI(p,i,b,a)) * eh_rho(i,a,n) * eh_rho(b,p,n) 
                 ! num  = - ERI(p,i,a,b) * ( eh_rho(a,i,nS+n) * eh_rho(b,p,n) - eh_rho(b,i,nS+n) * eh_rho(a,p,n) )
                 ! num  = - 0.5d0 * (ERI(p,i,a,b) - ERI(p,i,b,a)) * &
                 !      ( eh_rho(a,i,nS+n) * eh_rho(b,p,n) - eh_rho(b,i,nS+n) * eh_rho(a,p,n) )
                 dem1 = eQP(p) - eQP(b) - eh_Om(n) 
                 dem2 = eQP(p) + eQP(i) - eQP(a) - eQP(b) 
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 Sig_eh(p,4) = Sig_eh(p,4) + num * (reg1/dem1) * (reg2/dem2)
                 Z_eh(p,4)   = Z_eh(p,4)   - num * (reg1/dem1) * (reg2/dem2/dem2) &
                                           - num * (reg1/dem1/dem1) * (reg2/dem2)

                 ! 2h2p(d) * 2p1h(d)
                 num  = (ERI(p,a,i,b) - ERI(p,a,b,i)) * eh_rho(a,i,n) * eh_rho(b,p,n) 
                 ! num  = ERI(p,a,i,b) * ( eh_rho(i,a,nS+n) * eh_rho(b,p,n) - eh_rho(b,a,nS+n) * eh_rho(i,p,n) )
                 ! num  = 0.5d0 * (ERI(p,a,i,b) - ERI(p,a,b,i)) * &
                 !      ( eh_rho(i,a,nS+n) * eh_rho(b,p,n) - eh_rho(b,a,nS+n) * eh_rho(i,p,n) )
                 dem1 = eQP(a) - eQP(i) + eh_Om(n) 
                 dem2 = eQP(p) - eQP(b) - eh_Om(n)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 Sig_eh(p,5) = Sig_eh(p,5) + num * (reg1/dem1) * (reg2/dem2)
                 Z_eh(p,5)   = Z_eh(p,5)   - num * (reg1/dem1) * (reg2/dem2/dem2)
                 
                 ! 2h2p(d) * 2p1h
                 num  = (ERI(p,i,a,b) - ERI(p,i,b,a)) * eh_rho(a,i,n) * eh_rho(p,b,n) 
                 ! num  = ERI(p,i,a,b) * ( eh_rho(a,i,n) * eh_rho(b,p,nS+n) - eh_rho(b,i,n) * eh_rho(a,p,nS+n) )
                 ! num  = 0.5d0 * (ERI(p,i,a,b) - ERI(p,i,b,a)) * &
                 !      ( eh_rho(a,i,n) * eh_rho(b,p,nS+n) - eh_rho(b,i,n) * eh_rho(a,p,nS+n) )
                 dem1 = eQP(a) - eQP(i) + eh_Om(n) 
                 dem2 = eQP(p) + eQP(i) - eQP(a) - eQP(b)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 Sig_eh(p,6) = Sig_eh(p,6) + num * (reg1/dem1) * (reg2/dem2)
                 Z_eh(p,6)   = Z_eh(p,6)   - num * (reg1/dem1) * (reg2/dem2/dem2)
                 
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

  Sig_pp(:,:) = 0d0
  Z_pp(:,:)   = 0d0

  !$OMP PARALLEL DEFAULT(NONE)    &
  !$OMP PRIVATE(p,i,j,k,c,n,num,dem1,dem2,reg1,reg2) &
  !$OMP SHARED(nC,nO,nOrb,nR,nOO,nVV,eta,ERI,eQP,ee_rho,ee_Om,hh_rho,hh_Om,Sig_pp,Z_pp)
  !$OMP DO
  do p=nC+1,nOrb-nR
     
     do i=nC+1,nO
        do j=nC+1,nO
           do n=1,nVV

              do k=nC+1,nO

                 ! 2h2p(d) * 2p1h(d)
                 num  = - ERI(p,k,i,j) * ee_rho(i,j,n) * ee_rho(p,k,n)
                 dem1 = ee_Om(n) - eQP(i) - eQP(j)
                 dem2 = eQP(p) + eQP(k) - ee_Om(n)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 Sig_pp(p,5) = Sig_pp(p,5) + num * (reg1/dem1) * (reg2/dem2)
                 Z_pp(p,5)   = Z_pp(p,5)   - num * (reg1/dem1) * (reg2/dem2/dem2)

              end do ! k

              do c=nO+1,nOrb-nR

                 ! 2h2p(d) * 2h1p
                 num  = - ERI(p,c,i,j) * ee_rho(i,j,n) * ee_rho(p,c,n)
                 dem1 = ee_Om(n) - eQP(i) - eQP(j)
                 dem2 = eQP(p) + eQP(c) - eQP(i) - eQP(j)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 Sig_pp(p,3) = Sig_pp(p,3) + num * (reg1/dem1) * (reg2/dem2)
                 Z_pp(p,3)   = Z_pp(p,3)   - num * (reg1/dem1) * (reg2/dem2/dem2)
                 
              end do ! a
           end do ! n
           do n=1,nOO

              do c=nO+1,nOrb-nR

                 ! 2h1p(d) * 2h1p
                 num  = - ERI(p,c,i,j) * hh_rho(i,j,n) * hh_rho(p,c,n)
                 dem1 = eQP(p) + eQP(c) - hh_Om(n)
                 dem2 = eQP(p) + eQP(c) - eQP(i) - eQP(j)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 Sig_pp(p,1) = Sig_pp(p,1) + num * (reg1/dem1) * (reg2/dem2)
                 Z_pp(p,1)   = Z_pp(p,1)   - num * (reg1/dem1) * (reg2/dem2/dem2) &
                                           - num * (reg1/dem1/dem1) * (reg2/dem2)
                 
              end do ! c
           end do ! n
        end do ! j
     end do ! i
     
  end do ! p
  !$OMP END DO
  !$OMP END PARALLEL
  !$OMP PARALLEL DEFAULT(NONE)    &
  !$OMP PRIVATE(p,k,a,b,c,n,num,dem1,dem2,reg1,reg2) &
  !$OMP SHARED(nC,nO,nOrb,nR,nOO,nVV,eta,ERI,eQP,ee_rho,ee_Om,hh_rho,hh_Om,Sig_pp,Z_pp)
  !$OMP DO
  do p=nC+1,nOrb-nR
     do a=nO+1,nOrb-nR
        do b=nO+1,nOrb-nR
           do n=1,nOO
             
              do c=nO+1,nOrb-nR
              
                 ! 2p2h(d) * 2h1p(d)   
                 num  = ERI(p,c,a,b) * hh_rho(a,b,n) * hh_rho(p,c,n)
                 dem1 = hh_Om(n) - eQP(a) - eQP(b)
                 dem2 = eQP(p) + eQP(c) - hh_Om(n)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 Sig_pp(p,2) = Sig_pp(p,2) + num * (reg1/dem1) * (reg2/dem2)
                 Z_pp(p,2)   = Z_pp(p,2)   - num * (reg1/dem1) * (reg2/dem2/dem2)

              end do ! c
              
              do k=nC+1,nO

                 ! 2p2h(d) * 2p1h
                 num  = ERI(p,k,a,b) * hh_rho(a,b,n) * hh_rho(p,k,n)
                 dem1 = hh_Om(n) - eQP(a) - eQP(b)
                 dem2 = eQP(p) + eQP(k) - eQP(a) - eQP(b)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 Sig_pp(p,6) = Sig_pp(p,6) + num * (reg1/dem1) * (reg2/dem2)
                 Z_pp(p,6)   = Z_pp(p,6)   - num * (reg1/dem1) * (reg2/dem2/dem2)
                 
              end do ! k
           end do ! n
           do n=1,nVV
              
              do k=nC+1,nO

                 ! 2p1h * 2p1h(d)
                 num  = ERI(p,k,a,b) * ee_rho(a,b,n) * ee_rho(p,k,n)
                 dem1 = eQP(p) + eQP(k) - eQP(a) - eQP(b)
                 dem2 = eQP(p) + eQP(k) - ee_Om(n)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 Sig_pp(p,4) = Sig_pp(p,4) + num * (reg1/dem1) * (reg2/dem2)
                 Z_pp(p,4)   = Z_pp(p,4)   - num * (reg1/dem1) * (reg2/dem2/dem2) &
                                           - num * (reg1/dem1/dem1) * (reg2/dem2)
                 
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
 
!---------------------------!
! Self-energy decomposition !
!---------------------------!

  SigC(:) = Sig_2(:) &
          + Sig_eh(:,1) + Sig_eh(:,2) + Sig_eh(:,3) + Sig_eh(:,4)+ Sig_eh(:,5)+ Sig_eh(:,6) &
          + Sig_pp(:,1) + Sig_pp(:,2) + Sig_pp(:,3) + Sig_pp(:,4)+ Sig_pp(:,5)+ Sig_pp(:,6)
 
!------------------------!
! Renormalization factor !
!------------------------!

  Z(:) = Z_2(:) &
       + Z_eh(:,1) + Z_eh(:,2) + Z_eh(:,3) + Z_eh(:,4) + Z_eh(:,5) + Z_eh(:,6) &
       + Z_pp(:,1) + Z_pp(:,2) + Z_pp(:,3) + Z_pp(:,4) + Z_pp(:,5) + Z_pp(:,6)  

  Z(:) = 1d0/(1d0 - Z(:))

!---------------------------------!
! Print self-energy decomposition !
!---------------------------------!

  print_self_energy = .true.
  call dump_GParquet_self_energy(nOrb,nC,nO,nV,nR,Sig_2,Sig_eh,Sig_pp,SigC,Z_2,Z_eh,Z_pp,Z)
  
end subroutine 
