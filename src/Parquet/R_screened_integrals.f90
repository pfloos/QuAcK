subroutine R_eh_singlet_screened_integral(nOrb,nC,nO,nR,nS,ERI,eh_sing_Phi,eh_trip_Phi,pp_sing_Phi,pp_trip_Phi,XpY,XmY,rho)

! Compute excitation densities
  implicit none

! Input variables
  integer,intent(in)            :: nOrb,nC,nO,nR,nS
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: eh_sing_Phi(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: eh_trip_Phi(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: pp_sing_Phi(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: pp_trip_Phi(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: XpY(nS,nS),XmY(nS,nS)

! Local variables
  integer                       :: ia,jb,p,q,j,b
  double precision              :: X,Y

! Output variables
  double precision,intent(out)  :: rho(nOrb,nOrb,nS)
  
  rho(:,:,:) = 0d0   
!  !$OMP PARALLEL &
!  !$OMP SHARED(nC,nOrb,nR,nO,nS,rho,ERI,XpY,eh_sing_Gam) &
!  !$OMP PRIVATE(q,p,jb,ia) &
!  !$OMP DEFAULT(NONE)
!  !$OMP DO
  do q=nC+1,nOrb-nR
     do p=nC+1,nOrb-nR
              
        do ia=1,nS
        
           jb = 0
           do j=nC+1,nO
              do b=nO+1,nOrb-nR
                 jb = jb + 1

                 X = 0.5d0*(XpY(ia,jb) + XmY(ia,jb))
                 Y = 0.5d0*(XpY(ia,jb) - XmY(ia,jb))
                 
                 rho(p,q,ia) = rho(p,q,ia) &
                             + (2d0*ERI(p,b,q,j) - ERI(p,b,j,q)) * X       & 
                             + (- 0.5d0*eh_sing_Phi(p,b,j,q) - 1.5d0*eh_trip_Phi(p,b,j,q) &
                             + 0.5d0*pp_sing_Phi(p,b,q,j) + 1.5d0*pp_trip_Phi(p,b,q,j)) * X &
                             + (2d0*ERI(p,j,q,b) - ERI(p,j,b,q)) * Y              & 
                             + (- 0.5d0*eh_sing_Phi(p,j,b,q) - 1.5d0*eh_trip_Phi(p,j,b,q) &
                             + 0.5d0*pp_sing_Phi(p,j,q,b) + 1.5d0*pp_trip_Phi(p,j,q,b)) * Y
                 
              end do

           end do
        end do

     end do
  end do
!  !$OMP END DO
!  !$OMP END PARALLEL
  
end subroutine 

subroutine R_eh_triplet_screened_integral(nOrb,nC,nO,nR,nS,ERI,eh_sing_Phi,eh_trip_Phi,pp_sing_Phi,pp_trip_Phi,XpY,XmY,rho)

! Compute excitation densities
  implicit none

! Input variables
  integer,intent(in)            :: nOrb,nC,nO,nR,nS
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: eh_sing_Phi(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: eh_trip_Phi(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: pp_sing_Phi(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: pp_trip_Phi(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: XpY(nS,nS),XmY(nS,nS)

! Local variables
  integer                       :: ia,jb,p,q,j,b
  double precision              :: X,Y
  
! Output variables
  double precision,intent(out)  :: rho(nOrb,nOrb,nS)
  
  rho(:,:,:) = 0d0   
!  !$OMP PARALLEL &
!  !$OMP SHARED(nC,nOrb,nR,nO,nS,rho,ERI,XpY,eh_trip_Gam) &
!  !$OMP PRIVATE(q,p,jb,ia) &
!  !$OMP DEFAULT(NONE)
!  !$OMP DO
  do q=nC+1,nOrb-nR
     do p=nC+1,nOrb-nR
              
        do ia=1,nS
        
           jb = 0
           do j=nC+1,nO
              do b=nO+1,nOrb-nR
                 jb = jb + 1

                 X = 0.5d0*(XpY(ia,jb) + XmY(ia,jb))
                 Y = 0.5d0*(XpY(ia,jb) - XmY(ia,jb))
                 
                 rho(p,q,ia) = rho(p,q,ia) &
                             - ERI(p,b,j,q) * X       & 
                             + (- 0.5d0*eh_sing_Phi(p,b,j,q) + 0.5d0*eh_trip_Phi(p,b,j,q) &
                             - 0.5d0*pp_sing_Phi(p,b,q,j) + 0.5d0*pp_trip_Phi(p,b,q,j)) * X &
                             - ERI(p,j,b,q) * Y              & 
                             + (- 0.5d0*eh_sing_Phi(p,j,b,q) + 0.5d0*eh_trip_Phi(p,j,b,q) &
                             - 0.5d0*pp_sing_Phi(p,j,q,b) + 0.5d0*pp_trip_Phi(p,j,q,b)) * Y
                 
                 
              end do
           end do
           
        end do
        
     end do
  end do
!  !$OMP END DO
!  !$OMP END PARALLEL
  
end subroutine 


subroutine R_pp_singlet_screened_integral(nOrb,nC,nO,nR,nOO,nVV,ERI,eh_sing_Phi,eh_trip_Phi,X1,Y1,rho1,X2,Y2,rho2)

! Compute excitation densities in the singlet pp channel

  implicit none

! Input variables


  integer,intent(in)            :: nOrb,nC,nO,nR
  integer,intent(in)            :: nOO,nVV 
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: eh_sing_Phi(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: eh_trip_Phi(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: X1(nVV,nVV)
  double precision,intent(in)   :: Y1(nOO,nVV)
  double precision,intent(in)   :: X2(nVV,nOO)
  double precision,intent(in)   :: Y2(nOO,nOO)

! Local variables

  integer                       :: i,j,k,l
  integer                       :: a,b,c,d
  integer                       :: p,q
  integer                       :: ab,cd,ij,kl
  double precision,external     :: Kronecker_delta

! Output variables

  double precision,intent(out)  :: rho1(nOrb,nOrb,nVV)
  double precision,intent(out)  :: rho2(nOrb,nOrb,nOO)

  integer                       :: dim_1, dim_2

! Initialization

  rho1(:,:,:) = 0d0
  rho2(:,:,:) = 0d0

  ! !$OMP PARALLEL DEFAULT(NONE)                                &
  ! !$OMP PRIVATE(p, q, a, b, ab, c, d, cd, i, j, ij, k, l, kl) &
  ! !$OMP SHARED(nC, nOrb, nR, nO, rho1, rho2, ERI, eh_sing_Phi, eh_trip_Phi, X1, Y1, X2, Y2)
  ! !$OMP DO COLLAPSE(2)
  do q=nC+1,nOrb-nR
     do p=nC+1,nOrb-nR

        ab=0
        do a = nO+1, nOrb-nR
           do b = a, nOrb-nR
              ab = ab + 1
              
              cd = 0
              do c = nO+1, nOrb-nR
                 do d = c, nOrb-nR
                    cd = cd + 1
                    
                    rho1(p,q,ab) = rho1(p,q,ab) + (ERI(p,q,c,d) + ERI(p,q,d,c))*X1(cd,ab)/sqrt(1d0 + Kronecker_delta(c,d)) &
                                 + (0.5d0*eh_sing_Phi(p,q,c,d) - 1.5d0*eh_trip_Phi(p,q,c,d) &
                                 + 0.5d0*eh_sing_Phi(p,q,d,c) - 1.5d0*eh_trip_Phi(p,q,d,c))&
                                 *X1(cd,ab)/sqrt(1d0 + Kronecker_delta(c,d))
                    
                 end do ! d
              end do ! c
   
              kl = 0
              do k = nC+1, nO
                 do l = k, nO
                    
                    kl = kl + 1
  
                    rho1(p,q,ab) = rho1(p,q,ab) + (ERI(p,q,k,l) + ERI(p,q,l,k))*Y1(kl,ab)/sqrt(1d0 + Kronecker_delta(k,l))  &
                                 + (0.5d0*eh_sing_Phi(p,q,k,l) - 1.5d0*eh_trip_Phi(p,q,k,l) &
                                 + 0.5d0*eh_sing_Phi(p,q,l,k) - 1.5d0*eh_trip_Phi(p,q,l,k))&
                                 *Y1(kl,ab)/sqrt(1d0 + Kronecker_delta(k,l))
                 end do ! l
              end do ! k 
           end do ! b
        end do ! a
   
        ij = 0
        do i = nC+1, nO
           do j = i, nO
              ij = ij + 1
            
              cd = 0
              do c = nO+1, nOrb-nR
                 do d = c, nOrb-nR
                    cd = cd + 1
                   
                    rho2(p,q,ij) = rho2(p,q,ij) + (ERI(p,q,c,d) + ERI(p,q,d,c))*X2(cd,ij)/sqrt(1d0 + Kronecker_delta(c,d)) &
                                 + (0.5d0*eh_sing_Phi(p,q,c,d) - 1.5d0*eh_trip_Phi(p,q,c,d) &
                                 + 0.5d0*eh_sing_Phi(p,q,d,c) - 1.5d0*eh_trip_Phi(p,q,d,c))&
                                 *X2(cd,ij)/sqrt(1d0 + Kronecker_delta(c,d))
                 end do ! d
              end do ! c
             
              kl = 0
              do k = nC+1, nO
                 do l = k, nO
                    kl = kl + 1
                   
                    rho2(p,q,ij) = rho2(p,q,ij) + (ERI(p,q,k,l) + ERI(p,q,l,k))*Y2(kl,ij)/sqrt(1d0 + Kronecker_delta(k,l)) &
                                 + (0.5d0*eh_sing_Phi(p,q,k,l) - 1.5d0*eh_trip_Phi(p,q,k,l) &
                                 + 0.5d0*eh_sing_Phi(p,q,l,k) - 1.5d0*eh_trip_Phi(p,q,l,k))&
                                 *Y2(kl,ij)/sqrt(1d0 + Kronecker_delta(k,l))
                 end do ! l
              end do ! k
              
           end do ! j
        end do ! i   
        
     end do
  end do
  ! !$OMP END DO
  ! !$OMP END PARALLEL

end subroutine 




subroutine R_pp_triplet_screened_integral(nOrb,nC,nO,nR,nOO,nVV,ERI,eh_sing_Phi,eh_trip_Phi,X1,Y1,rho1,X2,Y2,rho2)

! Compute excitation densities in the triplet pp channel
  implicit none

! Input variables
  integer,intent(in)            :: nOrb,nC,nO,nR
  integer,intent(in)            :: nOO,nVV 
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: eh_sing_Phi(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: eh_trip_Phi(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: X1(nVV,nVV)
  double precision,intent(in)   :: Y1(nOO,nVV)
  double precision,intent(in)   :: X2(nVV,nOO)
  double precision,intent(in)   :: Y2(nOO,nOO)

! Local variables
  integer                       :: i,j,k,l
  integer                       :: a,b,c,d
  integer                       :: p,q
  integer                       :: ab,cd,ij,kl
  double precision,external     :: Kronecker_delta

! Output variables
  double precision,intent(out)  :: rho1(nOrb,nOrb,nVV)
  double precision,intent(out)  :: rho2(nOrb,nOrb,nOO)
  integer                       :: dim_1, dim_2

! Initialization
  rho1(:,:,:) = 0d0
  rho2(:,:,:) = 0d0

  dim_1 = (nOrb - nO) * (nOrb - nO - 1) / 2
  dim_2 = nO * (nO - 1) / 2

  ! !$OMP PARALLEL DEFAULT(NONE)                                &
  ! !$OMP PRIVATE(p, q, a, b, ab, c, d, cd, i, j, ij, k, l, kl) &
  ! !$OMP SHARED(nC, nOrb, nR, nO, rho1, rho2, ERI, eh_sing_Phi, eh_trip_Phi, X1, Y1, X2, Y2)
  ! !$OMP DO COLLAPSE(2)
  do q = nC+1, nOrb-nR
     do p = nC+1, nOrb-nR

        ab = 0
        do a = nO+1, nOrb-nR
           do b = a+1, nOrb-nR
              ab = ab + 1
              
              cd = 0
              do c = nO+1, nOrb-nR
                 do d = c+1, nOrb-nR
                    cd = cd + 1
                    
                    rho1(p,q,ab) = rho1(p,q,ab) + (ERI(p,q,c,d) - ERI(p,q,d,c))*X1(cd,ab) &
                                 + (0.5d0*eh_sing_Phi(p,q,c,d) + 0.5d0*eh_trip_Phi(p,q,c,d) &
                                 - 0.5d0*eh_sing_Phi(p,q,d,c) - 0.5d0*eh_trip_Phi(p,q,d,c) )*X1(cd,ab)
                    
                 end do ! d
              end do ! c
   
              kl = 0
              do k = nC+1, nO
                 do l = k+1, nO
                    
                    kl = kl + 1
  
                    rho1(p,q,ab) = rho1(p,q,ab) + (ERI(p,q,k,l) - ERI(p,q,l,k))*Y1(kl,ab) &
                                 + (0.5d0*eh_sing_Phi(p,q,k,l) + 0.5d0*eh_trip_Phi(p,q,k,l) &
                                 - 0.5d0*eh_sing_Phi(p,q,l,k) - 0.5d0*eh_trip_Phi(p,q,l,k) )*Y1(kl,ab)
                 end do ! l
              end do ! k 
           end do ! b
        end do ! a
   
        ij = 0
        do i = nC+1, nO
           do j = i+1, nO
              ij = ij + 1
            
              cd = 0
              do c = nO+1, nOrb-nR
                 do d = c+1, nOrb-nR
                    cd = cd + 1
                   
                    rho2(p,q,ij) = rho2(p,q,ij) + (ERI(p,q,c,d) - ERI(p,q,d,c))*X2(cd,ij) &
                                 + (0.5d0*eh_sing_Phi(p,q,c,d) + 0.5d0*eh_trip_Phi(p,q,c,d) &
                                 - 0.5d0*eh_sing_Phi(p,q,d,c) - 0.5d0*eh_trip_Phi(p,q,d,c) )*X2(cd,ij)
                 end do ! d
              end do ! c
             
              kl = 0
              do k = nC+1, nO
                 do l = k+1, nO
                    kl = kl + 1
                   
                    rho2(p,q,ij) = rho2(p,q,ij) + (ERI(p,q,k,l) - ERI(p,q,l,k))*Y2(kl,ij)                            &
                                 + (0.5d0*eh_sing_Phi(p,q,k,l) + 0.5d0*eh_trip_Phi(p,q,k,l) &
                                 - 0.5d0*eh_sing_Phi(p,q,l,k) - 0.5d0*eh_trip_Phi(p,q,l,k) )*Y2(kl,ij)
                 end do ! l
              end do ! k
              
           end do ! j
        end do ! i
        
     end do ! p
  end do ! q
  ! !$OMP END DO
  ! !$OMP END PARALLEL

end subroutine 
