subroutine R_pp_singlet_Phi(eta,nOrb,nC,nR,nOO,nVV,ee_sing_Om,ee_sing_rho,hh_sing_Om,hh_sing_rho,pp_sing_Phi)


! Compute irreducible vertex in the triplet pp channel
  implicit none

! Input variables
  double precision,intent(in)   :: eta
  integer,intent(in)            :: nOrb,nC,nR,nOO,nVV
  double precision,intent(in)   :: ee_sing_Om(nVV)
  double precision,intent(in)   :: ee_sing_rho(nOrb,nOrb,nVV)
  double precision,intent(in)   :: hh_sing_Om(nOO)
  double precision,intent(in)   :: hh_sing_rho(nOrb,nOrb,nOO)

! Local variables
  integer                       :: p,q,r,s
  integer                       :: n
  double precision,allocatable  :: tmp_ee(:,:,:),tmp_hh(:,:,:)

! Output variables
  double precision,intent(out)   :: pp_sing_Phi(nOrb,nOrb,nOrb,nOrb)

! Initialization
  pp_sing_Phi(:,:,:,:) = 0d0

  allocate(tmp_ee(nOrb,nOrb,nVV),tmp_hh(nOrb,nOrb,nOO))
  
  tmp_ee(:,:,:) = 0d0
  tmp_hh(:,:,:) = 0d0

  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(p, q, n) &
  !$OMP SHARED(eta, nC, nOrb, nR, nVV, nOO, tmp_ee, tmp_hh, ee_sing_rho, ee_sing_Om, hh_sing_rho, hh_sing_Om)
  !$OMP DO COLLAPSE(2)
  do q = nC+1, nOrb-nR
     do p = nC+1, nOrb-nR
        
        do n=1,nVV
           tmp_ee(p,q,n) = - (ee_sing_rho(p,q,n)/ee_sing_Om(n)) * (1d0 - exp(- 2d0 * eta * ee_sing_Om(n) * ee_sing_Om(n)))
        end do
        
        do n=1,nOO
           tmp_hh(p,q,n) = + (hh_sing_rho(p,q,n)/hh_sing_Om(n)) * (1d0 - exp(- 2d0 * eta * hh_sing_Om(n) * hh_sing_Om(n)))
        end do
        
     enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  
  call dgemm('N','T', nOrb*nOrb,  nOrb*nOrb, nVV, 1d0, tmp_ee(1,1,1),  nOrb*nOrb, ee_sing_rho(1,1,1),  nOrb*nOrb, 0d0, pp_sing_Phi(1,1,1,1),  nOrb*nOrb)

  call dgemm('N','T', nOrb*nOrb,  nOrb*nOrb, nOO, 1d0, tmp_hh(1,1,1),  nOrb*nOrb, hh_sing_rho(1,1,1),  nOrb*nOrb, 1d0, pp_sing_Phi(1,1,1,1),  nOrb*nOrb)
  
  deallocate(tmp_ee,tmp_hh)
  
  ! !$OMP PARALLEL DEFAULT(NONE) &
  ! !$OMP PRIVATE(p, q, r, s, n) &
  ! !$OMP SHARED(eta, nC, nOrb, nR, nVV, nOO, pp_sing_Phi, ee_sing_rho, ee_sing_Om, hh_sing_rho, hh_sing_Om)
  ! !$OMP DO COLLAPSE(2)
  ! do s = nC+1, nOrb-nR
  !    do r = nC+1, nOrb-nR
  !       do q = nC+1, nOrb-nR
  !          do p = nC+1, nOrb-nR

  !             do n=1,nVV
  !                pp_sing_Phi(p,q,r,s) = pp_sing_Phi(p,q,r,s) &
  !                     - (ee_sing_rho(p,q,n)*ee_sing_rho(r,s,n)/ee_sing_Om(n)) * (1d0 - exp(- 2d0 * eta * ee_sing_Om(n) * ee_sing_Om(n))) 
  !             end do

  !             do n=1,nOO
  !                pp_sing_Phi(p,q,r,s) = pp_sing_Phi(p,q,r,s) &
  !                     + (hh_sing_rho(p,q,n)*hh_sing_rho(r,s,n)/hh_sing_Om(n)) * (1d0 - exp(- 2d0 * eta * hh_sing_Om(n) * hh_sing_Om(n))) 
  !             end do
              
  !          enddo
  !       enddo
  !    enddo
  ! enddo
  ! !$OMP END DO
  ! !$OMP END PARALLEL

end subroutine 
