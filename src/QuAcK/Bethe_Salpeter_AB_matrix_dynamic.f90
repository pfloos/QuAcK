subroutine Bethe_Salpeter_AB_matrix_dynamic(eta,nBas,nC,nO,nV,nR,nS,lambda,eGW,OmRPA,OmBSE,rho,Ap,Am,Bp,Bm)

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
  double precision              :: chi_A,chi_B,eps
  double precision              :: chi_Ap,chi_Am,chi_Bp,chi_Bm
  double precision              :: eps_Ap,eps_Am,eps_Bp,eps_Bm
  integer                       :: i,j,a,b,ia,jb,kc

! Output variables

  double precision,intent(out)  :: Ap(nS,nS)
  double precision,intent(out)  :: Am(nS,nS)
  double precision,intent(out)  :: Bp(nS,nS)
  double precision,intent(out)  :: Bm(nS,nS)

! Initialization

  Ap(:,:) = 0d0
  Am(:,:) = 0d0
  Bp(:,:) = 0d0
  Bm(:,:) = 0d0

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

          Ap(ia,jb) = Ap(ia,jb) - 4d0*lambda*chi_A
          Am(ia,jb) = Am(ia,jb) - 4d0*lambda*chi_A
          Bp(ia,jb) = Bp(ia,jb) - 4d0*lambda*chi_B
          Bm(ia,jb) = Bm(ia,jb) - 4d0*lambda*chi_B

          chi_Ap = 0d0
          chi_Am = 0d0
          chi_Bp = 0d0
          chi_Bm = 0d0

          do kc=1,maxS

            eps_Ap = (+ OmBSE - OmRPA(kc) - (eGW(a) - eGW(j)))**2 + eta**2
            eps_Am = (+ OmBSE - OmRPA(kc) - (eGW(a) - eGW(j)))**2 + eta**2
            chi_Ap = chi_Ap + rho(i,j,kc)*rho(a,b,kc)*(+ OmBSE - OmRPA(kc) - (eGW(a) - eGW(j)))/eps_Ap
            chi_Am = chi_Am + rho(i,j,kc)*rho(a,b,kc)*(+ OmBSE - OmRPA(kc) - (eGW(a) - eGW(j)))/eps_Am

            eps_Ap = (+ OmBSE - OmRPA(kc) - (eGW(b) - eGW(i)))**2 + eta**2
            eps_Am = (+ OmBSE - OmRPA(kc) - (eGW(b) - eGW(i)))**2 + eta**2
            chi_Ap = chi_Ap + rho(i,j,kc)*rho(a,b,kc)*(+ OmBSE - OmRPA(kc) - (eGW(b) - eGW(i)))/eps_Ap
            chi_Am = chi_Am + rho(i,j,kc)*rho(a,b,kc)*(+ OmBSE - OmRPA(kc) - (eGW(b) - eGW(i)))/eps_Am

            eps_Bp = (+ OmBSE - OmRPA(kc) - (eGW(a) - eGW(b)))**2 + eta**2
            eps_Bm = (+ OmBSE - OmRPA(kc) - (eGW(a) - eGW(b)))**2 + eta**2
            chi_Bp = chi_Bp + rho(i,b,kc)*rho(a,j,kc)*(+ OmBSE - OmRPA(kc) - (eGW(a) - eGW(b)))/eps_Bp
            chi_Bm = chi_Bm + rho(i,b,kc)*rho(a,j,kc)*(+ OmBSE - OmRPA(kc) - (eGW(a) - eGW(b)))/eps_Bm

            eps_Bp = (+ OmBSE - OmRPA(kc) - (eGW(j) - eGW(i)))**2 + eta**2
            eps_Bm = (+ OmBSE - OmRPA(kc) - (eGW(j) - eGW(i)))**2 + eta**2
            chi_Bp = chi_Bp + rho(i,b,kc)*rho(a,j,kc)*(+ OmBSE - OmRPA(kc) - (eGW(j) - eGW(i)))/eps_Bp
            chi_Bm = chi_Bm + rho(i,b,kc)*rho(a,j,kc)*(+ OmBSE - OmRPA(kc) - (eGW(j) - eGW(i)))/eps_Bm

          enddo

          Ap(ia,jb) = Ap(ia,jb) - 2d0*lambda*chi_Ap
          Am(ia,jb) = Am(ia,jb) - 2d0*lambda*chi_Am

          Bp(ia,jb) = Bp(ia,jb) - 2d0*lambda*chi_Bp
          Bm(ia,jb) = Bm(ia,jb) - 2d0*lambda*chi_Bm

        enddo
      enddo
    enddo
  enddo

end subroutine Bethe_Salpeter_AB_matrix_dynamic
