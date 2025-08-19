subroutine R_Parquet_self_energy_eta(eta,nOrb,nC,nO,nV,nR,nS,nOOs,nVVs,nOOt,nVVt,eQP,ERI, &
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
  double precision              :: dem,dem1,dem2
  double precision              :: num,num1,num2
  double precision              :: start_t,end_t,t
  
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
                 
                 num = ERI(p,a,j,i) * ( 2d0*ERI(j,i,q,a) - ERI(j,i,a,q) )
                 dem = eQP(p) + eQP(a) - eQP(i) - eQP(j) 
                 
                 SigC(p,q) = SigC(p,q) + num*dem/(dem**2 + eta**2)
                 
              end do ! a
           end do ! j
        end do ! i
        ! 2p1h sum
        do i=nC+1,nO
           do a=nO+1,nOrb-nR
              do b=nO+1,nOrb-nR
                 
                 num = ERI(p,i,b,a) * ( 2d0*ERI(b,a,q,i) - ERI(b,a,i,q) )
                 dem = eQP(p) + eQP(i) - eQP(a) - eQP(b)
                 
                 SigC(p,q) = SigC(p,q) + num*dem/(dem**2 + eta**2)
                 
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
  !$OMP PARALLEL DEFAULT(NONE)    &
  !$OMP PRIVATE(p,q,i,a,j,b,n,num,dem1,dem2) &
  !$OMP SHARED(nC,nO,nOrb,nR,nS,eta,ERI,eQP,eh_sing_rho,eh_sing_Om,SigC,Z)
  !$OMP DO
  do q=nC+1,nOrb-nR
     do p=nC+1,nOrb-nR
     
        do i=nC+1,nO
           do a=nO+1,nOrb-nR
              
              do n=1,nS
                 !3h2p
                 do j=nC+1,nO

                    num  = (0.5d0*ERI(p,a,i,j) - ERI(p,a,j,i))* &
                         eh_sing_rho(i,a,n) * eh_sing_rho(q,j,n)
                    dem1 = eQP(p) - eQP(j) + eh_sing_Om(n)
                    dem2 = eQP(p) - eQP(i) - eQP(j) + eQP(a)
                    
                    SigC(p,q) = SigC(p,q) + num*dem1/(dem1**2 + eta**2)*dem2/(dem2**2 + eta**2)
                    
                    num  = (0.5d0*ERI(p,i,a,j) - ERI(p,i,j,a)) * &
                         eh_sing_rho(a,i,n) * eh_sing_rho(q,j,n)
                    dem1 = eQP(a) - eQP(i) + eh_sing_Om(n) 
                    dem2 = eQP(p) - eQP(j) + eh_sing_Om(n)
                    
                    SigC(p,q) = SigC(p,q) + num*dem1/(dem1**2 + eta**2)*dem2/(dem2**2 + eta**2)
                    
                    num  = (0.5d0*ERI(p,a,i,j) - ERI(p,a,j,i))* &
                         eh_sing_rho(a,i,n) * eh_sing_rho(j,q,n)
                    dem1 = eQP(a) - eQP(i) + eh_sing_Om(n) 
                    dem2 = eQP(p) - eQP(i) - eQP(j) + eQP(a)
                    
                    SigC(p,q) = SigC(p,q) + num*dem1/(dem1**2 + eta**2)*dem2/(dem2**2 + eta**2)
                    
                 end do ! j
                 !3p2h
                 do b=nO+1,nOrb-nR

                    num  = -(0.5d0*ERI(p,i,a,b) - ERI(p,i,b,a)) * &
                         eh_sing_rho(i,a,n) * eh_sing_rho(b,q,n)
                    dem1 = eQP(p) - eQP(b) - eh_sing_Om(n)
                    dem2 = eQP(p) + eQP(i) - eQP(a) - eQP(b)
                    
                    SigC(p,q) = SigC(p,q) + num*dem1/(dem1**2 + eta**2)*dem2/(dem2**2 + eta**2)
                    
                    num  = (0.5d0*ERI(p,a,i,b) - ERI(p,a,b,i)) * &
                         eh_sing_rho(a,i,n) * eh_sing_rho(b,q,n)
                    dem1 = eQP(a) - eQP(i) + eh_sing_Om(n) 
                    dem2 = eQP(p) - eQP(b) - eh_sing_Om(n)
                    
                    SigC(p,q) = SigC(p,q) + num*dem1/(dem1**2 + eta**2)*dem2/(dem2**2 + eta**2)
                    
                    num  = (0.5d0*ERI(p,i,a,b) - ERI(p,i,b,a)) * &
                         eh_sing_rho(a,i,n) * eh_sing_rho(q,b,n)
                    dem1 = eQP(a) - eQP(i) + eh_sing_Om(n) 
                    dem2 = eQP(p) + eQP(i) - eQP(a) - eQP(b)
                    
                    SigC(p,q) = SigC(p,q) + num*dem1/(dem1**2 + eta**2)*dem2/(dem2**2 + eta**2)
                    
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

  write(*,'(1X,A50,1X,F9.3,A8)') 'Wall time for singlet eh self-energy =',t,' seconds'
  write(*,*)

!-------------------------------------!
!  triplet eh part of the self-energy !
!-------------------------------------!
  
   call wall_time(start_t)
   !$OMP PARALLEL DEFAULT(NONE)    &
   !$OMP PRIVATE(p,q,i,a,j,b,n,num,dem1,dem2) &
   !$OMP SHARED(nC,nO,nOrb,nR,nS,eta,ERI,eQP,eh_trip_rho,eh_trip_Om,SigC,Z)
   !$OMP DO
   do q=nC+1,nOrb-nR
      do p=nC+1,nOrb-nR
    
         do i=nC+1,nO
            do a=nO+1,nOrb-nR
          
               do n=1,nS
                  !3h2p
                  do j=nC+1,nO

                     num  = (1.5d0*ERI(p,a,i,j) - 0d0*ERI(p,a,j,i))* &
                          eh_trip_rho(i,a,n) * eh_trip_rho(q,j,n)
                     dem1 = eQP(p) - eQP(j) + eh_trip_Om(n)
                     dem2 = eQP(p) - eQP(i) - eQP(j) + eQP(a)
                  
                     SigC(p,q) = SigC(p,q) + num*dem1/(dem1**2 + eta**2)*dem2/(dem2**2 + eta**2)
                                 
                     num  = (1.5d0*ERI(p,i,a,j) - 0d0*ERI(p,i,j,a)) * &
                          eh_trip_rho(a,i,n) * eh_trip_rho(q,j,n)
                     dem1 = eQP(a) - eQP(i) + eh_trip_Om(n) 
                     dem2 = eQP(p) - eQP(j) + eh_trip_Om(n)
                     
                     SigC(p,q) = SigC(p,q) + num*dem1/(dem1**2 + eta**2)*dem2/(dem2**2 + eta**2)
                     
                     num  = (1.5d0*ERI(p,a,i,j) - 0d0*ERI(p,a,j,i))* &
                          eh_trip_rho(a,i,n) * eh_trip_rho(j,q,n)
                     dem1 = eQP(a) - eQP(i) + eh_trip_Om(n) 
                     dem2 = eQP(p) - eQP(i) - eQP(j) + eQP(a)
                     
                     SigC(p,q) = SigC(p,q) + num*dem1/(dem1**2 + eta**2)*dem2/(dem2**2 + eta**2)
                     
                  end do ! j
                  !3p2h
                  do b=nO+1,nOrb-nR

                     num  = -(1.5d0*ERI(p,i,a,b) - 0d0*ERI(p,i,b,a)) * &
                          eh_trip_rho(i,a,n) * eh_trip_rho(b,q,n)
                     dem1 = eQP(p) - eQP(b) - eh_trip_Om(n)
                     dem2 = eQP(p) + eQP(i) - eQP(a) - eQP(b)
                     
                     SigC(p,q) = SigC(p,q) + num*dem1/(dem1**2 + eta**2)*dem2/(dem2**2 + eta**2)
                     
                     num  = (1.5d0*ERI(p,a,i,b) - 0d0*ERI(p,a,b,i)) * &
                          eh_trip_rho(a,i,n) * eh_trip_rho(b,q,n)
                     dem1 = eQP(a) - eQP(i) + eh_trip_Om(n) 
                     dem2 = eQP(p) - eQP(b) - eh_trip_Om(n)
                     
                     SigC(p,q) = SigC(p,q) + num*dem1/(dem1**2 + eta**2)*dem2/(dem2**2 + eta**2)
                     
                     num  = (1.5d0*ERI(p,i,a,b) - 0d0*ERI(p,i,b,a)) * &
                          eh_trip_rho(a,i,n) * eh_trip_rho(q,b,n)
                     dem1 = eQP(a) - eQP(i) + eh_trip_Om(n) 
                     dem2 = eQP(p) + eQP(i) - eQP(a) - eQP(b)
                     
                     SigC(p,q) = SigC(p,q) + num*dem1/(dem1**2 + eta**2)*dem2/(dem2**2 + eta**2)
                     
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

   write(*,'(1X,A50,1X,F9.3,A8)') 'Wall time for triplet eh self-energy =',t,' seconds'
   write(*,*)
  
!-------------------------------------!
!  singlet pp part of the self-energy !
!-------------------------------------!

  call wall_time(start_t)
  !$OMP PARALLEL DEFAULT(NONE)    &
  !$OMP PRIVATE(p,q,i,j,k,c,n,num,dem1,dem2) &
  !$OMP SHARED(nC,nO,nOrb,nR,nOOs,nVVs,eta,ERI,eQP,ee_sing_rho,ee_sing_Om,hh_sing_rho,hh_sing_Om,SigC,Z)
  !$OMP DO
  do q=nC+1,nOrb-nR
     do p=nC+1,nOrb-nR
     
        do i=nC+1,nO
           do j=nC+1,nO
              do n=1,nVVs
                 ! 4h1p
                 do k=nC+1,nO

                    num  = - 0.5d0 * ERI(p,k,i,j) * ee_sing_rho(i,j,n) * ee_sing_rho(q,k,n)
                    dem1 = ee_sing_Om(n) - eQP(i) - eQP(j)
                    dem2 = eQP(p) + eQP(k) - ee_sing_Om(n)
                    
                    SigC(p,q) = SigC(p,q) + num*dem1/(dem1**2 + eta**2)*dem2/(dem2**2 + eta**2)
                    
                 end do ! k
                 ! 3h2p
                 do c=nO+1,nOrb-nR

                    num  = - 0.5d0*ERI(p,c,i,j) * ee_sing_rho(i,j,n) * ee_sing_rho(q,c,n)
                    dem1 = ee_sing_Om(n) - eQP(i) - eQP(j)
                    dem2 = eQP(p) + eQP(c) - eQP(i) - eQP(j)
                    
                    SigC(p,q) = SigC(p,q) + num*dem1/(dem1**2 + eta**2)*dem2/(dem2**2 + eta**2)
                    
                 end do ! a
              end do ! n
              do n=1,nOOs
                 ! 3h2p
                 do c=nO+1,nOrb-nR

                    num  = - 0.5d0*ERI(p,c,i,j) * hh_sing_rho(i,j,n) * hh_sing_rho(q,c,n)
                    dem1 = eQP(p) + eQP(c) - hh_sing_Om(n)
                    dem2 = eQP(p) + eQP(c) - eQP(i) - eQP(j)
                    
                    SigC(p,q) = SigC(p,q) + num*dem1/(dem1**2 + eta**2)*dem2/(dem2**2 + eta**2)
                    
                 end do ! c
              end do ! n
           end do ! j
        end do ! i
        
     end do ! p
  end do ! q
  !$OMP END DO
  !$OMP END PARALLEL
  !$OMP PARALLEL DEFAULT(NONE)    &
  !$OMP PRIVATE(p,q,k,a,b,c,n,num,dem1,dem2) &
  !$OMP SHARED(nC,nO,nOrb,nR,nOOs,nVVs,eta,ERI,eQP,ee_sing_rho,ee_sing_Om,hh_sing_rho,hh_sing_Om,SigC,Z)
  !$OMP DO
  do q=nC+1,nOrb-nR
     do p=nC+1,nOrb-nR
        do a=nO+1,nOrb-nR
           do b=nO+1,nOrb-nR
              do n=1,nOOs
                 ! 4p1h
                 do c=nO+1,nOrb-nR

                    num  = 0.5d0*ERI(p,c,a,b) * hh_sing_rho(a,b,n) * hh_sing_rho(q,c,n)
                    dem1 = hh_sing_Om(n) - eQP(a) - eQP(b)
                    dem2 = eQP(p) + eQP(c) - hh_sing_Om(n)

                    SigC(p,q) = SigC(p,q) + num*dem1/(dem1**2 + eta**2)*dem2/(dem2**2 + eta**2)
                    
                 end do ! c
                 ! 3p2h
                 do k=nC+1,nO

                    num  = 0.5d0*ERI(p,k,a,b) * hh_sing_rho(a,b,n) * hh_sing_rho(q,k,n)
                    dem1 = hh_sing_Om(n) - eQP(a) - eQP(b)
                    dem2 = eQP(p) + eQP(k) - eQP(a) - eQP(b)
                    
                    SigC(p,q) = SigC(p,q) + num*dem1/(dem1**2 + eta**2)*dem2/(dem2**2 + eta**2)
                    
                 end do ! k
              end do ! n
              do n=1,nVVs
                 ! 3p2h
                 do k=nC+1,nO

                    num  = 0.5d0*ERI(p,k,a,b) * ee_sing_rho(a,b,n) * ee_sing_rho(q,k,n)
                    dem1 = eQP(p) + eQP(k) - eQP(a) - eQP(b)
                    dem2 = eQP(p) + eQP(k) - ee_sing_Om(n)
                    
                    SigC(p,q) = SigC(p,q) + num*dem1/(dem1**2 + eta**2)*dem2/(dem2**2 + eta**2)
                    
                 end do ! k
              end do ! n
           end do ! b
        end do ! a
        
     end do ! p
  end do ! q
  !$OMP END DO
  !$OMP END PARALLEL
  call wall_time(end_t)
  t = end_t - start_t

  write(*,'(1X,A50,1X,F9.3,A8)') 'Wall time for singlet pp self-energy =',t,' seconds'
  write(*,*)

!-------------------------------------!
!  triplet pp part of the self-energy !
!-------------------------------------!

  call wall_time(start_t)
  !$OMP PARALLEL DEFAULT(NONE)    &
  !$OMP PRIVATE(p,q,i,j,k,c,n,num,dem1,dem2) &
  !$OMP SHARED(nC,nO,nOrb,nR,nOOt,nVVt,eta,ERI,eQP,ee_trip_rho,ee_trip_Om,hh_trip_rho,hh_trip_Om,SigC,Z)
  !$OMP DO
  do q=nC+1,nOrb-nR
     do p=nC+1,nOrb-nR
    
        do i=nC+1,nO
           do j=nC+1,nO
              do n=1,nVVt
                 ! 4h1p
                 do k=nC+1,nO

                    num  = - 1.5d0 * ERI(p,k,i,j) * ee_trip_rho(i,j,n) * ee_trip_rho(q,k,n)
                    dem1 = ee_trip_Om(n) - eQP(i) - eQP(j)
                    dem2 = eQP(p) + eQP(k) - ee_trip_Om(n)
                    
                    SigC(p,q) = SigC(p,q) + num*dem1/(dem1**2 + eta**2)*dem2/(dem2**2 + eta**2)
                    
                 end do ! k
                 ! 3h2p
                 do c=nO+1,nOrb-nR

                    num  = - 1.5d0 * ERI(p,c,i,j) * ee_trip_rho(i,j,n) * ee_trip_rho(q,c,n)
                    dem1 = ee_trip_Om(n) - eQP(i) - eQP(j)
                    dem2 = eQP(p) + eQP(c) - eQP(i) - eQP(j)
                    
                    SigC(p,q) = SigC(p,q) + num*dem1/(dem1**2 + eta**2)*dem2/(dem2**2 + eta**2)
                    
                 end do ! c
              end do ! n
              do n=1,nOOt
                 ! 3h2p
                 do c=nO+1,nOrb-nR

                    num  = - 1.5d0 * ERI(p,c,i,j) * hh_trip_rho(i,j,n) * hh_trip_rho(q,c,n)
                    dem1 = eQP(p) + eQP(c) - hh_trip_Om(n)
                    dem2 = eQP(p) + eQP(c) - eQP(i) - eQP(j)

                    SigC(p,q) = SigC(p,q) + num*dem1/(dem1**2 + eta**2)*dem2/(dem2**2 + eta**2)
                 
                 end do ! c
              end do ! n
           end do ! j
        end do ! i
        
     end do ! p
  end do ! q
  !$OMP END DO
  !$OMP END PARALLEL
  !$OMP PARALLEL DEFAULT(NONE)    &
  !$OMP PRIVATE(p,q,k,a,b,c,n,num,dem1,dem2) &
  !$OMP SHARED(nC,nO,nOrb,nR,nOOt,nVVt,eta,ERI,eQP,ee_trip_rho,ee_trip_Om,hh_trip_rho,hh_trip_Om,SigC,Z)
  !$OMP DO
  do q=nC+1,nOrb-nR
     do p=nC+1,nOrb-nR
        do a=nO+1,nOrb-nR
           do b=nO+1,nOrb-nR
              do n=1,nOOt
                 ! 4p1h
                 do c=nO+1,nOrb-nR

                    num  = 1.5d0 * ERI(p,c,a,b) * hh_trip_rho(a,b,n) * hh_trip_rho(q,c,n)
                    dem1 = hh_trip_Om(n) - eQP(a) - eQP(b)
                    dem2 = eQP(p) + eQP(c) - hh_trip_Om(n)

                    SigC(p,q) = SigC(p,q) + num*dem1/(dem1**2 + eta**2)*dem2/(dem2**2 + eta**2)

                 end do ! c
                 ! 3p2h
                 do k=nC+1,nO

                    num  = 1.5d0 * ERI(p,k,a,b) * hh_trip_rho(a,b,n) * hh_trip_rho(q,k,n)
                    dem1 = hh_trip_Om(n) - eQP(a) - eQP(b)
                    dem2 = eQP(p) + eQP(k) - eQP(a) - eQP(b)

                    SigC(p,q) = SigC(p,q) + num*dem1/(dem1**2 + eta**2)*dem2/(dem2**2 + eta**2)

                 end do ! k
              end do ! n
              do n=1,nVVt
                 ! 3p2h
                 do k=nC+1,nO

                    num  = 1.5d0 * ERI(p,k,a,b) * ee_trip_rho(a,b,n) * ee_trip_rho(q,k,n)
                    dem1 = eQP(p) + eQP(k) - eQP(a) - eQP(b)
                    dem2 = eQP(p) + eQP(k) - ee_trip_Om(n)
                    
                    SigC(p,q) = SigC(p,q) + num*dem1/(dem1**2 + eta**2)*dem2/(dem2**2 + eta**2)
                 
                 end do ! k
              end do ! n
           end do ! b
        end do ! a
        
     end do ! p
  end do ! q
  !$OMP END DO
  !$OMP END PARALLEL 
  call wall_time(end_t)
  t = end_t - start_t

  write(*,'(1X,A50,1X,F9.3,A8)') 'Wall time for triplet pp self-energy =',t,' seconds'
  write(*,*)
 
!-----------------------------!
!   Renormalization factor    !
!-----------------------------!

  Z(:) = 1d0

end subroutine 
