subroutine GGTpp_ppBSE_static_kernel_D(eta,nOrb,nC,nO,nV,nR,nOO,nVV,lambda,ERI,eGF,Om1,rho1,Om2,rho2,T,KD_sta)

! Compute the OOOO block of the T-matrix static pp kernel

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nOO
  integer,intent(in)            :: nVV
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: eGF(nOrb)
  double precision,intent(in)   :: Om1(nVV)
  double precision,intent(in)   :: rho1(nOrb,nOrb,nVV)
  double precision,intent(in)   :: Om2(nOO)
  double precision,intent(in)   :: rho2(nOrb,nOrb,nOO)
  double precision,intent(in)   :: T(nOrb,nOrb,nOrb,nOrb)

! Local variables

  double precision              :: dem,num
  integer                       :: p,q,r,s,e,m
  integer                       :: i,j,k,l
  integer                       :: ij,kl,cd

! Output variables

  double precision,intent(out)  :: KD_sta(nOO,nOO)
  
! Initialization
  KD_sta(:,:) = 0d0

! Computing the kernel
! This is the same code as for the GF(2) kernel with elements T instead of ERI
  ij = 0
  do i=nC+1,nO
    do j=i+1,nO
      ij = ij + 1

      kl = 0
      do k=nC+1,nO
        do l=k+1,nO
          kl = kl + 1

          do m=nC+1,nO
            do e=nO+1,nOrb-nR
   
              dem = eGF(m) - eGF(e)
              num =       T(i,m,k,e) * T(e,j,m,l) + T(i,e,k,m) * T(m,j,e,l)
              num = num - T(j,m,k,e) * T(e,i,m,l) - T(j,e,k,m) * T(m,i,e,l)
                                                                               
              KD_sta(ij,kl) = KD_sta(ij,kl) + num*dem/(dem**2 + eta**2)
            end do
          end do
        
        end do
      end do
    end do
  end do
  
end subroutine 
