subroutine R_pp_triplet_Phi(eta,nOrb,nC,nR,nOO,nVV,ee_trip_Om,ee_trip_rho,hh_trip_Om,hh_trip_rho,omega,pp_trip_Phi)


! Compute irreducible vertex in the triplet pp channel
  implicit none

! Input variables
  double precision,intent(in)   :: eta
  integer,intent(in)            :: nOrb,nC,nR,nOO,nVV
  double precision,intent(in)   :: ee_trip_Om(nVV)
  double precision,intent(in)   :: ee_trip_rho(nOrb,nOrb,nVV)
  double precision,intent(in)   :: hh_trip_Om(nOO)
  double precision,intent(in)   :: hh_trip_rho(nOrb,nOrb,nOO)
  double precision,intent(in)   :: omega

! Local variables
  integer                       :: p,q,r,s
  integer                       :: n
  double precision              :: dem
  double precision,allocatable  :: tmp_ee(:,:,:),tmp_hh(:,:,:)

! Output variables
  double precision,intent(out)   :: pp_trip_Phi(nOrb,nOrb,nOrb,nOrb)

! Initialization
  pp_trip_Phi(:,:,:,:) = 0d0

  allocate(tmp_ee(nOrb,nOrb,nVV),tmp_hh(nOrb,nOrb,nOO))
  
  tmp_ee(:,:,:) = 0d0
  tmp_hh(:,:,:) = 0d0

  !$OMP PARALLEL DEFAULT(NONE) &
  !$OMP PRIVATE(p, q, n, dem) &
  !$OMP SHARED(eta, nC, nOrb, nR, nVV, nOO, tmp_ee, tmp_hh, ee_trip_rho, ee_trip_Om, hh_trip_rho, hh_trip_Om, omega)
  !$OMP DO COLLAPSE(2)
  do q = nC+1, nOrb-nR
     do p = nC+1, nOrb-nR
        
        do n=1,nVV
           dem = omega - ee_trip_Om(n)
           tmp_ee(p,q,n) = + (ee_trip_rho(p,q,n)/dem) * (1d0 - exp(- 2d0 * eta * dem * dem))
        end do
        
        do n=1,nOO
           dem = omega - hh_trip_Om(n)
           tmp_hh(p,q,n) = - (hh_trip_rho(p,q,n)/dem) * (1d0 - exp(- 2d0 * eta * dem * dem))
        end do
        
     enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  call dgemm('N','T', nOrb*nOrb,  nOrb*nOrb, nVV, 1d0, tmp_ee(1,1,1),  nOrb*nOrb, ee_trip_rho(1,1,1),  nOrb*nOrb, 0d0, pp_trip_Phi(1,1,1,1),  nOrb*nOrb)

  call dgemm('N','T', nOrb*nOrb,  nOrb*nOrb, nOO, 1d0, tmp_hh(1,1,1),  nOrb*nOrb, hh_trip_rho(1,1,1),  nOrb*nOrb, 1d0, pp_trip_Phi(1,1,1,1),  nOrb*nOrb)
  
  deallocate(tmp_ee,tmp_hh)

  ! This is the old code that uses openmp loops instead of dgemm to do the contraction
  
  ! !$OMP PARALLEL DEFAULT(NONE) &
  ! !$OMP PRIVATE(p, q, r, s, n) &
  ! !$OMP SHARED(eta, nC, nOrb, nR, nVV, nOO, pp_trip_Phi, ee_trip_rho, ee_trip_Om, hh_trip_rho, hh_trip_Om)
  ! !$OMP DO COLLAPSE(2)
  ! do s = nC+1, nOrb-nR
  !    do r = nC+1, nOrb-nR
  !       do q = nC+1, nOrb-nR
  !          do p = nC+1, nOrb-nR
              
  !             do n=1,nVV
  !                pp_trip_Phi(p,q,r,s) = pp_trip_Phi(p,q,r,s) &
  !                     - (ee_trip_rho(p,q,n)*ee_trip_rho(r,s,n)/ee_trip_Om(n)) * (1d0 - exp(- 2d0 * eta * ee_trip_Om(n) * ee_trip_Om(n)))
  !             end do

  !             do n=1,nOO
  !                pp_trip_Phi(p,q,r,s) = pp_trip_Phi(p,q,r,s) &
  !                     + (hh_trip_rho(p,q,n)*hh_trip_rho(r,s,n)/hh_trip_Om(n)) * (1d0 - exp(- 2d0 * eta * hh_trip_Om(n) * hh_trip_Om(n)))
  !             end do
              
  !          enddo
  !       enddo
  !    enddo
  ! enddo
  ! !$OMP END DO
  ! !$OMP END PARALLEL

end subroutine 
