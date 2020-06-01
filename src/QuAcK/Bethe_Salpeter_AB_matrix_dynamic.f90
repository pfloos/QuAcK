subroutine Bethe_Salpeter_AB_matrix_dynamic(eta,nBas,nC,nO,nV,nR,nS,lambda,eGW,OmRPA,OmBSE,rho,A_dyn,B_dyn)

! Compute the dynamic part of the Bethe-Salpeter equation matrices

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas,nC,nO,nV,nR,nS
  double precision,intent(in)   :: eta
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: eGW(nBas)
  double precision,intent(in)   :: OmRPA(nS)
  double precision,intent(in)   :: OmBSE
  double precision,intent(in)   :: rho(nBas,nBas,nS)
  
! Local variables

  integer                       :: maxS
  double precision              :: chi_A,chi_B,eps,eps_A,eps_B
  integer                       :: i,j,a,b,ia,jb,kc

! Output variables

  double precision,intent(out)  :: A_dyn(nS,nS)
  double precision,intent(out)  :: B_dyn(nS,nS)

! Initialization

  A_dyn(:,:) = 0d0
  B_dyn(:,:) = 0d0

! Number of poles taken into account 

  maxS = nS

! Build dynamic A matrix

  ia = 0
  do i=nC+1,nO
    do a=nO+1,nBas-nR
      ia = ia + 1
      jb = 0
      do j=nC+1,nO
        do b=nO+1,nBas-nR
          jb = jb + 1
 
          chi_A = 0d0
          chi_B = 0d0

          do kc=1,maxS

            eps = OmRPA(kc)**2 + eta**2
            chi_A = chi_A + rho(i,j,kc)*rho(a,b,kc)*OmRPA(kc)/eps
            chi_B = chi_B + rho(i,b,kc)*rho(a,j,kc)*OmRPA(kc)/eps

          enddo

          A_dyn(ia,jb) = A_dyn(ia,jb) - 4d0*lambda*chi_A
          B_dyn(ia,jb) = B_dyn(ia,jb) - 4d0*lambda*chi_B

          chi_A = 0d0
          chi_B = 0d0

          do kc=1,maxS

            eps_A = (+ OmBSE - OmRPA(kc) - (eGW(a) - eGW(j)))**2 + eta**2
            chi_A = chi_A + rho(i,j,kc)*rho(a,b,kc)*(+ OmBSE - OmRPA(kc) - (eGW(a) - eGW(j)))/eps_A

            eps_A = (+ OmBSE - OmRPA(kc) - (eGW(b) - eGW(i)))**2 + eta**2
            chi_A = chi_A + rho(i,j,kc)*rho(a,b,kc)*(+ OmBSE - OmRPA(kc) - (eGW(b) - eGW(i)))/eps_A

            eps_B = (+ OmBSE - OmRPA(kc) - (eGW(a) - eGW(b)))**2 + eta**2
            chi_B = chi_B + rho(i,b,kc)*rho(a,j,kc)*(+ OmBSE - OmRPA(kc) - (eGW(a) - eGW(b)))/eps_B

            eps_B = (+ OmBSE - OmRPA(kc) - (eGW(j) - eGW(i)))**2 + eta**2
            chi_B = chi_B + rho(i,b,kc)*rho(a,j,kc)*(+ OmBSE - OmRPA(kc) - (eGW(j) - eGW(i)))/eps_B

          enddo

          A_dyn(ia,jb) = A_dyn(ia,jb) - 2d0*lambda*chi_A

          B_dyn(ia,jb) = B_dyn(ia,jb) - 2d0*lambda*chi_B

        enddo
      enddo
    enddo
  enddo

end subroutine Bethe_Salpeter_AB_matrix_dynamic
