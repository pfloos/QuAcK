subroutine Sangalli_A_matrix_dynamic(eta,nBas,nC,nO,nV,nR,nS,lambda,eGW,OmRPA,OmBSE,rho,A_dyn)

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
  double precision              :: chi
  double precision              :: eps
  integer                       :: i,j,a,b,ia,jb,kc

! Output variables

  double precision,intent(out)  :: A_dyn(nS,nS)

! Initialization

  A_dyn(:,:) = 0d0

! Number of poles taken into account 

  maxS = nS

! Build dynamic A matrix

  do 


  (ERI(i,k,j,c)*KroneckerDelta(a,b) + ERI(k,a,c,b)*KroneckerDelta(i,j))*(X)

  ia = 0
  do i=nC+1,nO
    do a=nO+1,nBas-nR
      ia = ia + 1
      jb = 0
      do j=nC+1,nO
        do b=nO+1,nBas-nR
          jb = jb + 1
 
          chi = 0d0
          do kc=1,maxS

            chi = chi + rho(i,j,kc)*rho(a,b,kc)*OmRPA(kc)/(OmRPA(kc)**2 + eta**2)

          enddo

          A_dyn(ia,jb) = A_dyn(ia,jb) - 4d0*lambda*chi

          chi = 0d0
          do kc=1,maxS
            do ld=1,maxS

            eps = OmBSE - (OmRPA(kc) + OmRPA(ld))
            chi = chi + cRPA(ia,kc,ld)*cRPA(jb,kc,ld)*eps/(eps**2 + (2d0*eta)**2)

          enddo

          A_dyn(ia,jb) = A_dyn(ia,jb) - 2d0*lambda*chi

        enddo
      enddo
    enddo
  enddo

end subroutine Sangalli_A_matrix_dynamic
