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
  double precision,allocatable  :: X(:,:),Y(:,:)
  double precision,allocatable  :: Gamma_eh_bj(:,:,:),Gamma_eh_jb(:,:,:)
  
! Output variables
  double precision,intent(out)  :: rho(nOrb,nOrb,nS)
  
  rho(:,:,:) = 0d0

  allocate(X(nS,nS),Y(nS,nS))

  X(:,:) = 0.5d0*(XpY(:,:) + XmY(:,:))
  Y(:,:) = 0.5d0*(XpY(:,:) - XmY(:,:))
  
  allocate(Gamma_eh_bj(nOrb,nOrb,nS),Gamma_eh_jb(nOrb,nOrb,nS))
  
  Gamma_eh_bj(:,:,:) = 0d0
  Gamma_eh_jb(:,:,:) = 0d0

  !$OMP PARALLEL &
  !$OMP SHARED(nC,nOrb,nR,nO,ERI,eh_sing_Phi,eh_trip_Phi,pp_sing_Phi,pp_trip_Phi,Gamma_eh_jb,Gamma_eh_bj) &
  !$OMP PRIVATE(q,p,jb) &
  !$OMP DEFAULT(NONE)
  !$OMP DO COLLAPSE(2)
  do q=nC+1,nOrb-nR
     do p=nC+1,nOrb-nR
        
           jb = 0
           do j=nC+1,nO
              do b=nO+1,nOrb-nR
                 jb = jb + 1

                 
                 Gamma_eh_bj(p,q,jb) = Gamma_eh_bj(p,q,jb)                              &
                             + ( - ERI(p,b,j,q)                                         & 
                             - 0.5d0*eh_sing_Phi(p,b,j,q) + 0.5d0*eh_trip_Phi(p,b,j,q)  &
                             - 0.5d0*pp_sing_Phi(p,b,q,j) + 0.5d0*pp_trip_Phi(p,b,q,j) )
                 Gamma_eh_jb(p,q,jb) = Gamma_eh_jb(p,q,jb)                              &
                             + ( - ERI(p,j,b,q)                                         & 
                             - 0.5d0*eh_sing_Phi(p,j,b,q) + 0.5d0*eh_trip_Phi(p,j,b,q)  &
                             - 0.5d0*pp_sing_Phi(p,j,q,b) + 0.5d0*pp_trip_Phi(p,j,q,b) )
                 
           end do
        end do

     end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  ! rho = Gamma_eh_bj * X^T(ia,bj) + Gamma_eh_jb * Y^T(ia,jb)

  call dgemm('N','T', nOrb*nOrb, nS, nS, 1d0, Gamma_eh_bj(1,1,1), nOrb * nOrb, X(1,1), nS, 0d0, rho(1,1,1), nOrb*nOrb)

  call dgemm('N','T', nOrb*nOrb, nS, nS, 1d0, Gamma_eh_jb(1,1,1), nOrb * nOrb, Y(1,1), nS, 1d0, rho(1,1,1), nOrb*nOrb)
  
! This is the old code that uses openmp loops instead of dgemm to do the contraction
  
  ! !$OMP PARALLEL &
  ! !$OMP SHARED(nC,nOrb,nR,nO,nS,rho,ERI,XpY,XmY,eh_sing_Phi,eh_trip_Phi,pp_sing_Phi,pp_trip_Phi) &
  ! !$OMP PRIVATE(q,p,jb,ia,X,Y) &
  ! !$OMP DEFAULT(NONE)
  ! !$OMP DO COLLAPSE(2)
  ! do q=nC+1,nOrb-nR
  !    do p=nC+1,nOrb-nR
              
  !       do ia=1,nS
        
  !          jb = 0
  !          do j=nC+1,nO
  !             do b=nO+1,nOrb-nR
  !                jb = jb + 1

  !                X = 0.5d0*(XpY(ia,jb) + XmY(ia,jb))
  !                Y = 0.5d0*(XpY(ia,jb) - XmY(ia,jb))
                 
  !                rho(p,q,ia) = rho(p,q,ia)                                                   &
  !                            + ( - ERI(p,b,j,q)                                              & 
  !                            - 0.5d0*eh_sing_Phi(p,b,j,q) + 0.5d0*eh_trip_Phi(p,b,j,q)       &
  !                            - 0.5d0*pp_sing_Phi(p,b,q,j) + 0.5d0*pp_trip_Phi(p,b,q,j) ) * X &
  !                            + ( - ERI(p,j,b,q)                                              & 
  !                            - 0.5d0*eh_sing_Phi(p,j,b,q) + 0.5d0*eh_trip_Phi(p,j,b,q)       &
  !                            - 0.5d0*pp_sing_Phi(p,j,q,b) + 0.5d0*pp_trip_Phi(p,j,q,b) ) * Y
                 
                 
  !             end do
  !          end do
           
  !       end do
        
  !    end do
  ! end do
  ! !$OMP END DO
  ! !$OMP END PARALLEL
  
end subroutine 
