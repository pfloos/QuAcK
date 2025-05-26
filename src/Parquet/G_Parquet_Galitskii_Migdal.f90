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
  integer                       :: n
  double precision              :: eps,dem1,dem2,reg,reg1,reg2
  double precision              :: num
  double precision              :: start_t,end_t,t

  double precision              :: Ec_2
  double precision              :: Ec_eh(6)
  double precision              :: Ec_pp(6)

! Output variables

  double precision,intent(out)  :: EcGM
  
!-----------------------------------!
! 2nd-order part of the self-energy !
!-----------------------------------!

  call wall_time(start_t)

  Ec_2 = 0d0

  do i=nC+1,nO
     do j=nC+1,nO
        do a=nO+1,nOrb-nR
           do b=nO+1,nOrb-nR

              eps = eQP(i) + eQP(j) - eQP(a) - eQP(b)
              reg = (1d0 - exp(- 2d0 * eta * eps * eps))
              num = 0.5d0*(ERI(a,b,i,j) - ERI(a,b,j,i))**2

              Ec_2 = Ec_2 + num*reg/eps

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

  call wall_time(start_t)

  Ec_eh(:) = 0d0


  ! !$OMP SHARED(nC,nO,nOrb,nR,nS,eta,ERI,eQP,eh_rho,eh_Om,SigC,Z)
  ! !$OMP DO COLLAPSE(2)
  do i=nC+1,nO
     do j=nC+1,nO
        do a=nO+1,nOrb-nR
           do b=nO+1,nOrb-nR
              do n=1,nS

                 ! 2h2p(d) * 2p2h
                 num  = 0.5d0 * (ERI(b,a,i,j) - ERI(b,a,j,i)) * eh_rho(i,a,n) * eh_rho(j,b,nS+n)
                 dem1 = eQP(b) - eQP(j) + eh_Om(n)
                 dem2 = eQP(i) + eQP(j) - eQP(a) - eQP(b) 
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 Ec_eh(1) = Ec_eh(1) + num * (reg1/dem1) * (reg2/dem2)

                 ! 2h2p(d) * 2p2h(d)
                 num  = 0.5d0 * (ERI(b,i,a,j) - ERI(b,i,j,a)) * eh_rho(a,i,n) * eh_rho(j,b,nS+n)
                 dem1 = eQP(a) - eQP(i) + eh_Om(n) 
                 dem2 = eQP(j) - eQP(b) - eh_Om(n)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 Ec_eh(2) = Ec_eh(2) + num * (reg1/dem1) * (reg2/dem2)

                 ! 2h2p(d) * 2p2h
                 num  = 0.5d0 * (ERI(b,a,i,j) - ERI(b,a,j,i)) * eh_rho(i,a,nS+n) * eh_rho(j,b,n)
                 dem1 = eQP(a) - eQP(i) + eh_Om(n) 
                 dem2 = eQP(i) + eQP(j) - eQP(a) - eQP(b)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 Ec_eh(3) = Ec_eh(3) + num * (reg1/dem1) * (reg2/dem2)
                
                 ! 2p2h(d) * 2h2p
                 num  = 0.5d0 * (ERI(j,i,a,b) - ERI(j,i,b,a)) * eh_rho(a,i,nS+n) * eh_rho(b,j,n)
                 dem1 = eQP(j) - eQP(b) - eh_Om(n) 
                 dem2 = eQP(a) + eQP(b) - eQP(i) - eQP(j)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 Ec_eh(4) = Ec_eh(4) + num * (reg1/dem1) * (reg2/dem2)

                 ! 2h2p(d) * 2p2h(d)
                 num  = 0.5d0 * (ERI(j,a,i,b) - ERI(j,a,b,i)) * eh_rho(i,a,nS+n) * eh_rho(b,j,n)
                 dem1 = eQP(a) - eQP(i) + eh_Om(n) 
                 dem2 = eQP(j) - eQP(b) - eh_Om(n)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 Ec_eh(5) = Ec_eh(5) + num * (reg1/dem1) * (reg2/dem2)

                 ! 2h2p(d) * 2p2h
                 num  = 0.5d0 * (ERI(j,i,a,b) - ERI(j,i,b,a)) * eh_rho(a,i,n) * eh_rho(b,j,nS+n)
                 dem1 = eQP(a) - eQP(i) + eh_Om(n)
                 dem2 = eQP(i) + eQP(j) - eQP(a) - eQP(b)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 Ec_eh(6) = Ec_eh(6) + num * (reg1/dem1) * (reg2/dem2)

             end do ! n
           end do ! b
        end do ! a
     end do ! j
  end do ! i
   ! !$OMP END DO
   ! !$OMP END PARALLEL

   call wall_time(end_t)
   t = end_t - start_t

   write(*,'(1X,A50,1X,F9.3,A8)') 'Wall time for computing eh correlation energy =',t,' seconds'
   write(*,*) 

!-----------------------------!
!  pp part of the self-energy !
!-----------------------------!

  call wall_time(start_t)

  Ec_pp(:) = 0d0

  ! !$OMP PARALLEL DEFAULT(NONE)    &
  ! !$OMP PRIVATE(p,i,j,k,l,a,b,n,num,dem1,dem2,reg1,reg2) &
  ! !$OMP SHARED(nC,nO,nOrb,nR,nOO,nVV,eta,ERI,eQP,ee_rho,ee_Om,hh_rho,hh_Om,Ec_pp)
  ! !$OMP DO COLLAPSE(2)
  do i=nC+1,nO
     do j=nC+1,nO
        do n=1,nVV

           do k=nC+1,nO
              do l=nC+1,nO

                 ! 2h2p(d) * 2p2h(d)
                 num  = - 0.5d0 * ERI(k,l,i,j) * ee_rho(i,j,n) * ee_rho(k,l,n)
                 dem1 = ee_Om(n) - eQP(i) - eQP(j)
                 dem2 = eQP(k) + eQP(l) - ee_Om(n)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 Ec_pp(1) = Ec_pp(1) + num * (reg1/dem1) * (reg2/dem2)

              end do ! l
           end do ! k
           
           do a=nO+1,nOrb-nR
              do b=nO+1,nOrb-nR

                 ! 2h2p(d) * 2p2h
                 num  = - 0.5d0 * ERI(a,b,i,j) * ee_rho(i,j,n) * ee_rho(a,b,n)
                 dem1 = ee_Om(n) - eQP(i) - eQP(j)
                 dem2 = eQP(i) + eQP(j) - eQP(a) - eQP(b)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 Ec_pp(2) = Ec_pp(2) + num * (reg1/dem1) * (reg2/dem2)
                 
              end do ! a
           end do ! a
           
        end do ! n
        do n=1,nOO
           
           do a=nO+1,nOrb-nR
              do b=nO+1,nOrb-nR

                 ! 2h2p(d) * 2p2h
                 num  = - 0.5d0 * ERI(a,b,i,j) * hh_rho(i,j,n) * hh_rho(a,b,n)
                 dem1 = eQP(a) + eQP(b) - hh_Om(n)
                 dem2 = eQP(i) + eQP(j) - eQP(a) - eQP(b)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 Ec_pp(3) = Ec_pp(3) + num * (reg1/dem1) * (reg2/dem2)
                 
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
                 
                 ! 2p2h(d) * 2p2h(d)
                 num  = 0.5d0 * ERI(c,d,a,b) * hh_rho(a,b,n) * hh_rho(c,d,n)
                 dem1 = hh_Om(n) - eQP(a) - eQP(b)
                 dem2 = hh_Om(n) - eQP(c) - eQP(d)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 Ec_pp(4) = Ec_pp(4) + num * (reg1/dem1) * (reg2/dem2)

              end do ! d
           end do ! c
           
           do i=nC+1,nO
              do j=nC+1,nO

                 ! 2p2h(d) * 2p2h
                 num  = 0.5d0 * ERI(i,j,a,b) * hh_rho(a,b,n) * hh_rho(i,j,n)
                 dem1 = hh_Om(n) - eQP(a) - eQP(b)
                 dem2 = eQP(i) + eQP(j) - eQP(a) - eQP(b)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 Ec_pp(5) = Ec_pp(5) + num * (reg1/dem1) * (reg2/dem2)
                 
              end do ! j
           end do ! i
           
        end do ! n
        do n=1,nVV

           do i=nC+1,nO
              do j=nC+1,nO

                 ! 2p2h * 2p2h(d)
                 num  = 0.5d0 * ERI(i,j,a,b) * ee_rho(a,b,n) * ee_rho(i,j,n)
                 dem1 = eQP(i) + eQP(j) - eQP(a) - eQP(b)
                 dem2 = eQP(i) + eQP(j) - ee_Om(n)
                 reg1 = (1d0 - exp(- 2d0 * eta * dem1 * dem1))
                 reg2 = (1d0 - exp(- 2d0 * eta * dem2 * dem2))

                 Ec_pp(6) = Ec_pp(6) + num * (reg1/dem1) * (reg2/dem2)
                 
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

!--------------------------!
! Total correlation energy !
!--------------------------!

  EcGM = Ec_2 + sum(Ec_eh) + sum(Ec_pp)

!----------------------------------!
! Correlation energy decomposition !
!----------------------------------!

  ! eh decomposition

  write(*,*) '-------------------------------------------------'
  write(*,*) '| Ec^eh components (in au)                      |'
  write(*,*) '-------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A20,1X,A1,1X,A16,1X,A1,1X)') '|','#','|','','|','Ec^eh','|'
  write(*,*) '-------------------------------------------------'
  write(*,'(1X,A1,1X,I3,1X,A1,1X,A20,1X,A1,1X,F16.10,1X,A1,1X)') '|',1,'|','2h2p(d)*2h2p','|',Ec_eh(1),'|'
  write(*,'(1X,A1,1X,I3,1X,A1,1X,A20,1X,A1,1X,F16.10,1X,A1,1X)') '|',2,'|','2h2p(d)*2h2p(d)','|',Ec_eh(2),'|'
  write(*,'(1X,A1,1X,I3,1X,A1,1X,A20,1X,A1,1X,F16.10,1X,A1,1X)') '|',3,'|','2h2p(d)*2h2p','|',Ec_eh(3),'|'
  write(*,'(1X,A1,1X,I3,1X,A1,1X,A20,1X,A1,1X,F16.10,1X,A1,1X)') '|',4,'|','2p2h(d)*2p2h','|',Ec_eh(4),'|'
  write(*,'(1X,A1,1X,I3,1X,A1,1X,A20,1X,A1,1X,F16.10,1X,A1,1X)') '|',5,'|','2h2p(d)*2p2h(d)','|',Ec_eh(5),'|'
  write(*,'(1X,A1,1X,I3,1X,A1,1X,A20,1X,A1,1X,F16.10,1X,A1,1X)') '|',6,'|','2h2p(d)*2p2h','|',Ec_eh(6),'|'
  write(*,*) '-------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A20,1X,A1,1X,F16.10,1X,A1,1X)') '|','Tot','|','','|',sum(Ec_eh),'|'
  write(*,*) '-------------------------------------------------'
  write(*,*) 
  
  ! pp decomposition

  write(*,*) '-------------------------------------------------'
  write(*,*) '| Ec^pp components (in au)                      |'
  write(*,*) '-------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A20,1X,A1,1X,A16,1X,A1,1X)') '|','#','|','','|','Ec^pp','|'
  write(*,*) '-------------------------------------------------'
  write(*,'(1X,A1,1X,I3,1X,A1,1X,A20,1X,A1,1X,F16.10,1X,A1,1X)') '|',1,'|','2h2p(d)*2p2h(d)','|',Ec_pp(1),'|'
  write(*,'(1X,A1,1X,I3,1X,A1,1X,A20,1X,A1,1X,F16.10,1X,A1,1X)') '|',2,'|','2h2p(d)*2p2h','|',Ec_pp(2),'|'
  write(*,'(1X,A1,1X,I3,1X,A1,1X,A20,1X,A1,1X,F16.10,1X,A1,1X)') '|',3,'|','2h2p(d)*2p2h','|',Ec_pp(3),'|'
  write(*,'(1X,A1,1X,I3,1X,A1,1X,A20,1X,A1,1X,F16.10,1X,A1,1X)') '|',4,'|','2p2h(d)*2p2h(d)','|',Ec_pp(4),'|'
  write(*,'(1X,A1,1X,I3,1X,A1,1X,A20,1X,A1,1X,F16.10,1X,A1,1X)') '|',5,'|','2p2h(d)*2p2h','|',Ec_pp(5),'|'
  write(*,'(1X,A1,1X,I3,1X,A1,1X,A20,1X,A1,1X,F16.10,1X,A1,1X)') '|',6,'|','2p2h*2p2h(d)','|',Ec_pp(6),'|'
  write(*,*) '-------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A20,1X,A1,1X,F16.10,1X,A1,1X)') '|','Tot','|','','|',sum(Ec_pp),'|'
  write(*,*) '-------------------------------------------------'
  write(*,*) 
  
  write(*,*) '-----------------------------------------------------'
  write(*,*) '| Parquet correlation energy decomposition (in au)  |'
  write(*,*) '-----------------------------------------------------'
  write(*,'(1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X)') & 
          '|','Ec^(2)','|','Ec^eh','|','Ec^pp','|','Total','|'
  write(*,*) '-----------------------------------------------------'
  write(*,'(1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X)') & 
          '|',Ec_2,'|',sum(Ec_eh),'|',sum(Ec_pp),'|',EcGM,'|'
  write(*,*) '-----------------------------------------------------'
  write(*,*)
  
end subroutine 
