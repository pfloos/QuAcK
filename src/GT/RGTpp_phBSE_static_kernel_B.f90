subroutine RGTpp_phBSE_static_kernel_B(ispin,eta,nBas,nC,nO,nV,nR,nS,nOOs,nVVs,nOOt,nVVt,lambda,Om1s,rho1s,Om2s,rho2s,Om1t,rho1t,Om2t,rho2t,KB)

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

  double precision,intent(out)  :: KB(nS,nS)

  KB(:,:) = 0d0

! Build B matrix for single manifold

  if(ispin == 1) then 
  
  jb = 0
  !$omp parallel do default(private) shared(KB,Om1s,rho1s,Om2s,rho2s,Om1t,rho1t,Om2t,rho2t,nO,nC,nR,nBas,nVVs,nOOs,nVVt,nOOt,chi,eps,eta,lambda)
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
            chi = chi + 0.5d0*rho1s(i,j,cd)*rho1s(a,b,cd)*eps/(eps**2 + eta**2)

          end do

          do kl=1,nOOs
            eps = + Om2s(kl)
            chi = chi + 0.5d0*rho2s(i,j,kl)*rho2s(a,b,kl)*eps/(eps**2 + eta**2)
          end do

          do cd=1,nVVt
            eps = - Om1t(cd)
            chi = chi + 1.5d0*rho1t(i,j,cd)*rho1t(a,b,cd)*eps/(eps**2 + eta**2)

          end do

          do kl=1,nOOt
            eps = + Om2t(kl)
            chi = chi + 1.5d0*rho2t(i,j,kl)*rho2t(a,b,kl)*eps/(eps**2 + eta**2)
          end do
          
          KB(ia,jb) = lambda*chi

        end do
      end do
    end do
  end do

  !$omp end parallel do

  end if

! Build B matrix for triplet manifold
  
  if(ispin == 2) then

     jb = 0
     !$omp parallel do default(private) shared(KB,Om1s,rho1s,Om2s,rho2s,Om1t,rho1t,Om2t,rho2t,nO,nC,nR,nBas,nVVs,nOOs,nVVt,nOOt,chi,eps,eta,lambda)
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
                    chi = chi - 0.5d0*rho1s(i,j,cd)*rho1s(a,b,cd)*eps/(eps**2 + eta**2)
                    
                 end do
                 
                 do kl=1,nOOs
                    eps = + Om2s(kl)
                    chi = chi - 0.5d0*rho2s(i,j,kl)*rho2s(a,b,kl)*eps/(eps**2 + eta**2)
                 end do
                 
                 do cd=1,nVVt
                    eps = - Om1t(cd)
                    chi = chi + 0.5d0*rho1t(i,j,cd)*rho1t(a,b,cd)*eps/(eps**2 + eta**2)
                 end do

                 do kl=1,nOOt
                    eps = + Om2t(kl)
                    chi = chi + 0.5d0*rho2t(i,j,kl)*rho2t(a,b,kl)*eps/(eps**2 + eta**2)
                 end do
                 
                 KB(ia,jb) = lambda*chi
                 
              end do
           end do
        end do
     end do
     
     !$omp end parallel do
     
  end if
  
end subroutine 
