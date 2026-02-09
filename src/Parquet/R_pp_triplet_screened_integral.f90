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
  double precision,allocatable  :: Gamma_pp_cd(:,:,:),Gamma_pp_kl(:,:,:)

! Output variables
  double precision,intent(out)  :: rho1(nOrb,nOrb,nVV)
  double precision,intent(out)  :: rho2(nOrb,nOrb,nOO)

! Initialization
  rho1(:,:,:) = 0d0
  rho2(:,:,:) = 0d0

  allocate(Gamma_pp_cd(nOrb,nOrb,nVV),Gamma_pp_kl(nOrb,nOrb,nOO))
  
  Gamma_pp_cd(:,:,:) = 0d0
  Gamma_pp_kl(:,:,:) = 0d0

  !$OMP PARALLEL DEFAULT(NONE)                                &
  !$OMP PRIVATE(p, q, c, d, cd, k, l, kl) &
  !$OMP SHARED(nC, nOrb, nR, nO, Gamma_pp_cd, Gamma_pp_kl, ERI, eh_sing_Phi, eh_trip_Phi)
  !$OMP DO COLLAPSE(2)
  do q = nC+1, nOrb-nR
     do p = nC+1, nOrb-nR
        
        cd = 0
        do c = nO+1, nOrb-nR
           do d = c+1, nOrb-nR
              cd = cd + 1
              
              Gamma_pp_cd(p,q,cd) = Gamma_pp_cd(p,q,cd) + ( ERI(p,q,c,d) - ERI(p,q,d,c)      &
                                  + 0.5d0*eh_sing_Phi(p,q,c,d) + 0.5d0*eh_trip_Phi(p,q,c,d)  &
                                  - 0.5d0*eh_sing_Phi(p,q,d,c) - 0.5d0*eh_trip_Phi(p,q,d,c) )
              
           end do ! d
        end do ! c
        
        kl = 0
        do k = nC+1, nO
           do l = k+1, nO
              kl = kl + 1
                    
              Gamma_pp_kl(p,q,kl) = Gamma_pp_kl(p,q,kl) + ( ERI(p,q,k,l) - ERI(p,q,l,k)      &
                                  + 0.5d0*eh_sing_Phi(p,q,k,l) + 0.5d0*eh_trip_Phi(p,q,k,l)  &
                                  - 0.5d0*eh_sing_Phi(p,q,l,k) - 0.5d0*eh_trip_Phi(p,q,l,k) )
              
           end do ! l
        end do ! k 

     end do ! p
  end do ! q
  !$OMP END DO
  !$OMP END PARALLEL

  ! rho1 = Gamma_pp_cd * X1(cd,ab) + Gamma_pp_kl * Y1(kl,ab)

  call dgemm('N','N', nOrb*nOrb, nVV, nVV, 1d0, Gamma_pp_cd(1,1,1), nOrb * nOrb, X1(1,1), nVV, 0d0, rho1(1,1,1), nOrb*nOrb)

  call dgemm('N','N', nOrb*nOrb, nVV, nOO, 1d0, Gamma_pp_kl(1,1,1), nOrb * nOrb, Y1(1,1), nOO, 1d0, rho1(1,1,1), nOrb*nOrb)
  
  ! rho2 = Gamma_pp_cd * X2(cd,ij) + Gamma_pp_kl * Y1(kl,ij)

  call dgemm('N','N', nOrb*nOrb, nOO, nVV, 1d0, Gamma_pp_cd(1,1,1), nOrb * nOrb, X2(1,1), nVV, 0d0, rho2(1,1,1), nOrb*nOrb)

  call dgemm('N','N', nOrb*nOrb, nOO, nOO, 1d0, Gamma_pp_kl(1,1,1), nOrb * nOrb, Y2(1,1), nOO, 1d0, rho2(1,1,1), nOrb*nOrb)
  
  deallocate(Gamma_pp_cd,Gamma_pp_kl)

! This is the old code that uses openmp loops instead of dgemm to do the contraction
  
  ! !$OMP PARALLEL DEFAULT(NONE)                                &
  ! !$OMP PRIVATE(p, q, a, b, ab, c, d, cd, i, j, ij, k, l, kl) &
  ! !$OMP SHARED(nC, nOrb, nR, nO, rho1, rho2, ERI, eh_sing_Phi, eh_trip_Phi, X1, Y1, X2, Y2)
  ! !$OMP DO COLLAPSE(2)
  ! do q = nC+1, nOrb-nR
  !    do p = nC+1, nOrb-nR

  !       ab = 0
  !       do a = nO+1, nOrb-nR
  !          do b = a+1, nOrb-nR
  !             ab = ab + 1
              
  !             cd = 0
  !             do c = nO+1, nOrb-nR
  !                do d = c+1, nOrb-nR
  !                   cd = cd + 1
                    
  !                   rho1(p,q,ab) = rho1(p,q,ab) + ( ERI(p,q,c,d) - ERI(p,q,d,c)            &
  !                                + 0.5d0*eh_sing_Phi(p,q,c,d) + 0.5d0*eh_trip_Phi(p,q,c,d) &
  !                                - 0.5d0*eh_sing_Phi(p,q,d,c) - 0.5d0*eh_trip_Phi(p,q,d,c) )*X1(cd,ab)
                    
  !                end do ! d
  !             end do ! c
   
  !             kl = 0
  !             do k = nC+1, nO
  !                do l = k+1, nO
                    
  !                   kl = kl + 1
  
  !                   rho1(p,q,ab) = rho1(p,q,ab) + ( ERI(p,q,k,l) - ERI(p,q,l,k)            &
  !                                + 0.5d0*eh_sing_Phi(p,q,k,l) + 0.5d0*eh_trip_Phi(p,q,k,l) &
  !                                - 0.5d0*eh_sing_Phi(p,q,l,k) - 0.5d0*eh_trip_Phi(p,q,l,k) )*Y1(kl,ab)
  !                end do ! l
  !             end do ! k 
  !          end do ! b
  !       end do ! a
   
  !       ij = 0
  !       do i = nC+1, nO
  !          do j = i+1, nO
  !             ij = ij + 1
            
  !             cd = 0
  !             do c = nO+1, nOrb-nR
  !                do d = c+1, nOrb-nR
  !                   cd = cd + 1
                   
  !                   rho2(p,q,ij) = rho2(p,q,ij) + ( ERI(p,q,c,d) - ERI(p,q,d,c)            &
  !                                + 0.5d0*eh_sing_Phi(p,q,c,d) + 0.5d0*eh_trip_Phi(p,q,c,d) &
  !                                - 0.5d0*eh_sing_Phi(p,q,d,c) - 0.5d0*eh_trip_Phi(p,q,d,c) )*X2(cd,ij)
  !                end do ! d
  !             end do ! c
             
  !             kl = 0
  !             do k = nC+1, nO
  !                do l = k+1, nO
  !                   kl = kl + 1
                   
  !                   rho2(p,q,ij) = rho2(p,q,ij) + ( ERI(p,q,k,l) - ERI(p,q,l,k)            &
  !                                + 0.5d0*eh_sing_Phi(p,q,k,l) + 0.5d0*eh_trip_Phi(p,q,k,l) &
  !                                - 0.5d0*eh_sing_Phi(p,q,l,k) - 0.5d0*eh_trip_Phi(p,q,l,k) )*Y2(kl,ij)
  !                end do ! l
  !             end do ! k
              
  !          end do ! j
  !       end do ! i
        
  !    end do ! p
  ! end do ! q
  ! !$OMP END DO
  ! !$OMP END PARALLEL

end subroutine 
