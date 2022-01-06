subroutine Bethe_Salpeter_B_matrix(eta,nBas,nC,nO,nV,nR,nS,lambda,ERI,OmRPA,rho_RPA,B_lr)

! Compute the extra term for Bethe-Salpeter equation for linear response 

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas,nC,nO,nV,nR,nS
  double precision,intent(in)   :: eta
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: OmRPA(nS)
  double precision,intent(in)   :: rho_RPA(nBas,nBas,nS)
  
! Local variables

  double precision              :: chi
  double precision              :: eps
  integer                       :: i,j,a,b,ia,jb,kc

! Output variables

  double precision,intent(out)  :: B_lr(nS,nS)

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
            eps = OmRPA(kc)**2 + eta**2
!           chi = chi + lambda*rho_RPA(i,b,kc)*rho_RPA(a,j,kc)*OmRPA(kc)/eps
            chi = chi + rho_RPA(i,b,kc)*rho_RPA(a,j,kc)*OmRPA(kc)/eps
          enddo

          B_lr(ia,jb) = B_lr(ia,jb) - lambda*ERI(i,j,b,a) + 4d0*lambda*chi

        enddo
      enddo
    enddo
  enddo

end subroutine Bethe_Salpeter_B_matrix
