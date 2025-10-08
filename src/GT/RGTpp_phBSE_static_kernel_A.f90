subroutine RGTpp_phBSE_static_kernel_A(ispin,eta,nBas,nC,nO,nV,nR,nS,nOOs,nVVs,nOOt,nVVt,lambda,Om1s,rho1s,Om2s,rho2s,Om1t,rho1t,Om2t,rho2t,KA)

! Compute the OOVV block of the static T-matrix

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: ispin
  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  integer,intent(in)            :: nOOs,nVVs
  integer,intent(in)            :: nOOt,nVVt
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: Om1s(nVVs)
  double precision,intent(in)   :: rho1s(nBas,nBas,nVVs)
  double precision,intent(in)   :: Om2s(nOOs)
  double precision,intent(in)   :: rho2s(nBas,nBas,nOOs)
  double precision,intent(in)   :: Om1t(nVVt)
  double precision,intent(in)   :: rho1t(nBas,nBas,nVVt)
  double precision,intent(in)   :: Om2t(nOOt)
  double precision,intent(in)   :: rho2t(nBas,nBas,nOOt)

! Local variables

  double precision              :: chi
  double precision              :: eps
  integer                       :: i,j,a,b,ia,jb,kl,cd,c,d

! Output variables

  double precision,intent(out)  :: KA(nS,nS)

  KA(:,:) = 0d0


! Build A matrix for single manifold

  if(ispin == 1) then 
  
  jb = 0
  !$omp parallel do default(private) shared(KA,Om1s,rho1s,Om2s,rho2s,Om1t,rho1t,Om2t,rho2t,nO,nC,nR,nBas,nVVs,nOOs,nVVt,nOOt,chi,eps,eta,lambda)
  do j=nC+1,nO
    do b=nO+1,nBas-nR
      jb = (b-nO) + (j-1)*(nBas-nO) 

      ia = 0
      do i=nC+1,nO
        do a=nO+1,nBas-nR
          ia = (a-nO) + (i-1)*(nBas-nO) 

          chi = 0d0

          do cd=1,nVVs
            eps = - Om1s(cd)
            chi = chi + 0.5d0*rho1s(i,b,cd)*rho1s(a,j,cd)*eps/(eps**2 + eta**2)

          end do

          do kl=1,nOOs
            eps = + Om2s(kl)
            chi = chi + 0.5d0*rho2s(i,b,kl)*rho2s(a,j,kl)*eps/(eps**2 + eta**2)
          end do

          do cd=1,nVVt
            eps = - Om1t(cd)
            chi = chi + 1.5d0*rho1t(i,b,cd)*rho1t(a,j,cd)*eps/(eps**2 + eta**2)

          end do

          do kl=1,nOOt
            eps = + Om2t(kl)
            chi = chi + 1.5d0*rho2t(i,b,kl)*rho2t(a,j,kl)*eps/(eps**2 + eta**2)
          end do
          
          KA(ia,jb) = lambda*chi

        end do
      end do
    end do
  end do

  !$omp end parallel do

  end if

! Build A matrix for triplet manifold
  
  if(ispin == 2) then

     jb = 0
     !$omp parallel do default(private) shared(KA,Om1s,rho1s,Om2s,rho2s,Om1t,rho1t,Om2t,rho2t,nO,nC,nR,nBas,nVVs,nOOs,nVVt,nOOt,chi,eps,eta,lambda)
     do j=nC+1,nO
        do b=nO+1,nBas-nR
           jb = (b-nO) + (j-1)*(nBas-nO) 
           
           ia = 0
           do i=nC+1,nO
              do a=nO+1,nBas-nR
                 ia = (a-nO) + (i-1)*(nBas-nO) 
                 
                 chi = 0d0
                 
                 do cd=1,nVVs
                    eps = - Om1s(cd)
                    chi = chi - 0.5d0*rho1s(i,b,cd)*rho1s(a,j,cd)*eps/(eps**2 + eta**2)
                    
                 end do
                 
                 do kl=1,nOOs
                    eps = + Om2s(kl)
                    chi = chi - 0.5d0*rho2s(i,b,kl)*rho2s(a,j,kl)*eps/(eps**2 + eta**2)
                 end do
                 
                 do cd=1,nVVt
                    eps = - Om1t(cd)
                    chi = chi + 0.5d0*rho1t(i,b,cd)*rho1t(a,j,cd)*eps/(eps**2 + eta**2)
                 end do

                 do kl=1,nOOt
                    eps = + Om2t(kl)
                    chi = chi + 0.5d0*rho2t(i,b,kl)*rho2t(a,j,kl)*eps/(eps**2 + eta**2)
                 end do
                 
                 KA(ia,jb) = lambda*chi
                 
              end do
           end do
        end do
     end do
     
     !$omp end parallel do
     
  end if
  
end subroutine 
