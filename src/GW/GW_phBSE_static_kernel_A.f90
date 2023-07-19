subroutine GW_phBSE_static_kernel_A(eta,nBas,nC,nO,nV,nR,nS,lambda,ERI,Om,rho,KA)

! Compute the BSE static kernel for the resonant block

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas,nC,nO,nV,nR,nS
  double precision,intent(in)   :: eta
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: Om(nS)
  double precision,intent(in)   :: rho(nBas,nBas,nS)

! Local variables

  double precision              :: chi
  double precision              :: eps
  integer                       :: i,j,a,b,ia,jb,kc

! Output variables

  double precision,intent(out)  :: KA(nS,nS)

! Compute static kernel

  ia = 0
  do i=nC+1,nO
    do a=nO+1,nBas-nR
      ia = ia + 1
      jb = 0
      do j=nC+1,nO
        do b=nO+1,nBas-nR
          jb = jb + 1

          chi = 0d0
          do kc=1,nS
            eps = Om(kc)**2 + eta**2
            chi = chi + rho(i,j,kc)*rho(a,b,kc)*Om(kc)/eps
          enddo

          KA(ia,jb) = 4d0*lambda*chi

        enddo
      enddo
    enddo
  enddo

end subroutine 
