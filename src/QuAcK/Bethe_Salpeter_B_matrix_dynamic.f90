subroutine Bethe_Salpeter_B_matrix_dynamic(eta,nBas,nC,nO,nV,nR,nS,lambda,OmRPA,OmBSE,rho,B_lr)

! Compute the dynamic part of the Bethe-Salpeter equation matrices

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas,nC,nO,nV,nR,nS
  double precision,intent(in)   :: eta
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: OmRPA(nS)
  double precision,intent(in)   :: OmBSE(nS)
  double precision,intent(in)   :: rho(nBas,nBas,nS)
  
! Local variables

  double precision              :: chi
  double precision              :: eps
  integer                       :: i,j,a,b,ia,jb,kc

! Output variables

  double precision,intent(out)  :: B_lr(nS,nS)

  B_lr(:,:) = 0d0

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
            chi = chi + rho(i,b,kc)*rho(a,j,kc)*OmRPA(kc)/eps

          enddo

          B_lr(ia,jb) = B_lr(ia,jb) - 4d0*lambda*chi

          chi = 0d0
          do kc=1,nS

            eps = (OmBSE(kc) - OmRPA(kc))**2 + eta**2
            chi = chi + rho(i,b,kc)*rho(a,j,kc)*(OmBSE(kc) - OmRPA(kc))/eps

            eps = (OmBSE(kc) + OmRPA(kc))**2 + eta**2
            chi = chi - rho(i,b,kc)*rho(a,j,kc)*(OmBSE(kc) + OmRPA(kc))/eps


          enddo

          B_lr(ia,jb) = B_lr(ia,jb) - 2d0*lambda*chi

        enddo
      enddo
    enddo
  enddo

end subroutine Bethe_Salpeter_B_matrix_dynamic
