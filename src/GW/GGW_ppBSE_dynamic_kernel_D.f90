subroutine GGW_ppBSE_dynamic_kernel_D(eta,nBas,nC,nO,nV,nR,nS,nOO,lambda,eGW,Om,rho,OmBSE,KD_dyn,ZD_dyn)

! Compute the dynamic part of the Bethe-Salpeter equation matrices

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  integer,intent(in)            :: nOO
  double precision,intent(in)   :: eta
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: eGW(nBas)
  double precision,intent(in)   :: Om(nS)
  double precision,intent(in)   :: rho(nBas,nBas,nS)
  double precision,intent(in)   :: OmBSE
  
! Local variables

  double precision,external     :: Kronecker_delta
  double precision              :: dem,num
  integer                       :: m
  integer                       :: i,j,k,l
  integer                       :: ij,kl

! Output variables

  double precision,intent(out)  :: KD_dyn(nOO,nOO)
  double precision,intent(out)  :: ZD_dyn(nOO,nOO)

! Initialization

  KD_dyn(:,:) = 0d0
  ZD_dyn(:,:) = 0d0

! Build dynamic D matrix

  ij = 0
  do i=nC+1,nO
    do j=i+1,nO
      ij = ij + 1
 
      kl = 0
      do k=nC+1,nO
        do l=k+1,nO
          kl = kl + 1
  
          do m=1,nS
             
               num = (rho(i,k,m)*rho(j,l,m) - rho(j,k,m)*rho(i,l,m))/2
               dem = - OmBSE - Om(m) + eGW(j) + eGW(l)
               KD_dyn(ij,kl) = KD_dyn(ij,kl) + num*dem/(dem**2 + eta**2)
               ZD_dyn(ij,kl) = ZD_dyn(ij,kl) + num*(dem**2 - eta**2)/(dem**2 + eta**2)**2

               dem = - OmBSE - Om(m) + eGW(i) + eGW(k)
               KD_dyn(ij,kl) = KD_dyn(ij,kl) + num*dem/(dem**2 + eta**2)
               ZD_dyn(ij,kl) = ZD_dyn(ij,kl) + num*(dem**2 - eta**2)/(dem**2 + eta**2)**2

               dem = - OmBSE - Om(m) + eGW(i) + eGW(l)
               KD_dyn(ij,kl) = KD_dyn(ij,kl) + num*dem/(dem**2 + eta**2)
               ZD_dyn(ij,kl) = ZD_dyn(ij,kl) + num*(dem**2 - eta**2)/(dem**2 + eta**2)**2
       
               dem = - OmBSE - Om(m) + eGW(j) + eGW(k)
               KD_dyn(ij,kl) = KD_dyn(ij,kl) + num*dem/(dem**2 + eta**2)
               ZD_dyn(ij,kl) = ZD_dyn(ij,kl) + num*(dem**2 - eta**2)/(dem**2 + eta**2)**2

          end do            

          KD_dyn(ij,kl) = 0.5d0*KD_dyn(ij,kl)
          ZD_dyn(ij,kl) = 0.5d0*ZD_dyn(ij,kl)
          
        end do
      end do

    end do
  end do

end subroutine 
