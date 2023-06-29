subroutine XBSE_static_kernel_KA(eta,nBas,nC,nO,nV,nR,nS,lambda,ERI,Om,rho,WA,eW,eGW)

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
  double precision,intent(in)   :: eW(nBas)
  double precision,intent(in)   :: eGW(nBas)

! Local variables

  double precision              :: chi
  double precision              :: eps
  double precision              :: num,den,ei,ea
  integer                       :: i,j,k,a,b,c,ia,jb,m
  double precision,external     :: Kronecker_delta

! Output variables

  double precision,intent(out)  :: WA(nS,nS)

! Initialize 

  WA(:,:) = 0d0

  ia = 0
  do i=nC+1,nO
    do a=nO+1,nBas-nR
      ia = ia + 1
      jb = 0
      do j=nC+1,nO
        do b=nO+1,nBas-nR
          jb = jb + 1

          ! virtual quasiparticle term

          ea = 0d0
          do m=1,nS
            do k=nC+1,nO
              num = 1d0*rho(a,k,m)*rho(b,k,m)
              den = eW(a) - eW(k) + Om(m)

              ea = ea + Kronecker_delta(i,j)*num*den/(den**2+eta**2) 

              den = eW(b) - eW(k) + Om(m)

              ea = ea + Kronecker_delta(i,j)*num*den/(den**2+eta**2) 

            end do
            do c=nO+1,nBas-nR
              num = 1d0*rho(a,c,m)*rho(b,c,m)
              den = eW(a) - eW(c) - Om(m)

              ea = ea + Kronecker_delta(i,j)*num*den/(den**2+eta**2) 

              den = eW(b) - eW(c) - Om(m)

              ea = ea + Kronecker_delta(i,j)*num*den/(den**2+eta**2) 

            end do
          end do

          ! occupied quasiparticle term

          ei = 0d0
          do m=1,nS
            do k=nC+1,nO
              num = 1d0*rho(i,k,m)*rho(j,k,m)
              den = eW(i) - eW(k) + Om(m)

              ei = ei + Kronecker_delta(a,b)*num*den/(den**2+eta**2) 

              den = eW(j) - eW(k) + Om(m)

              ei = ei + Kronecker_delta(a,b)*num*den/(den**2+eta**2) 

            end do
            do c=nO+1,nBas-nR
              num = 1d0*rho(i,c,m)*rho(j,c,m)
              den = eW(i) - eW(c) - Om(m)

              ei = ei + Kronecker_delta(a,b)*num*den/(den**2+eta**2) 

              den = eW(j) - eW(c) - Om(m)

              ei = ei + Kronecker_delta(a,b)*num*den/(den**2+eta**2) 

            end do
          end do

          ! kernel term

          chi = 0d0
          do m=1,nS
            eps = Om(m)**2 + eta**2
            chi = chi + rho(i,j,m)*rho(a,b,m)*Om(m)/eps!&
!                     - rho(i,b,m)*rho(a,j,m)*Om(m)/eps
          enddo

!         WA(ia,jb) = WA(ia,jb) + lambda*ERI(i,b,j,a) - 4d0*lambda*chi &
!                   - (eGW(a) - eW(a))*Kronecker_delta(i,j)*Kronecker_delta(a,b) &
!                   + (eGW(i) - eW(i))*Kronecker_delta(i,j)*Kronecker_delta(a,b)
          WA(ia,jb) = WA(ia,jb) - ea + ei + lambda*ERI(i,b,j,a) - 4d0*lambda*chi

        enddo
      enddo
    enddo
  enddo

end subroutine 
