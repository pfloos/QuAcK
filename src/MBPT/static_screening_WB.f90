subroutine static_screening_WB(eta,nBas,nC,nO,nV,nR,nS,lambda,ERI,Omega,rho,WB)

! Compute the static screening W for the coupling block

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas,nC,nO,nV,nR,nS
  double precision,intent(in)   :: eta
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: Omega(nS)
  double precision,intent(in)   :: rho(nBas,nBas,nS)

! Local variables

  double precision              :: chi
  double precision              :: eps
  integer                       :: i,j,a,b,ia,jb,kc

! Output variables

  double precision,intent(out)  :: WB(nS,nS)

! Initialize 

  WB(:,:) = 0d0

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
            eps = Omega(kc)**2 + eta**2
            chi = chi + rho(i,b,kc)*rho(a,j,kc)*Omega(kc)/eps
          enddo

          WB(ia,jb) = WB(ia,jb) + lambda*ERI(i,j,b,a) - 4d0*lambda*chi

        enddo
      enddo
    enddo
  enddo

end subroutine static_screening_WB
