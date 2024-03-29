subroutine ppULR_D(ispin,nBas,nC,nO,nV,nR,nHaa,nHab,nHbb,nHt,lambda,e,ERI_aaaa,ERI_aabb,ERI_bbbb,Dpp)

! Compute linear response

  implicit none
  include 'parameters.h'

! Input variables
 
  integer,intent(in)            :: ispin
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC(nspin)
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nV(nspin)
  integer,intent(in)            :: nR(nspin)
  integer,intent(in)            :: nHaa
  integer,intent(in)            :: nHab
  integer,intent(in)            :: nHbb
  integer,intent(in)            :: nHt
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: e(nBas,nspin)
  double precision,intent(in)   :: ERI_aaaa(nBas,nBas,nBas,nBas) 
  double precision,intent(in)   :: ERI_aabb(nBas,nBas,nBas,nBas) 
  double precision,intent(in)   :: ERI_bbbb(nBas,nBas,nBas,nBas) 
  
! Local variables
  
  double precision              :: eF 
  double precision,external     :: Kronecker_delta

  integer                       :: i,j,k,l,ij,kl

! Output variables

  double precision,intent(out)  :: Dpp(nHt,nHt)


  eF = 0d0 

  if(ispin == 3) then 

    ! abab block

    ij = 0
    do i=nC(1)+1,nO(1)
      do j=nC(2)+1,nO(2)
        ij = ij + 1
        kl = 0
        do k=nC(1)+1,nO(1)
          do l=nC(2)+1,nO(2)
            kl = kl + 1
            Dpp(ij,kl) = -(e(i,1) + e(j,2))*Kronecker_delta(i,k)&
                          *Kronecker_delta(j,l) +lambda*ERI_aabb(i,j,k,l)
          end  do
        end  do
      end  do
    end  do

  end if

  if(ispin == 4) then 

    ! aaaa block

    ij = 0
    do i=nC(1)+1,nO(1)
      do j=i+1,nO(1)
        ij = ij + 1
        kl = 0
        do k=nC(1)+1,nO(1)
          do l=k+1,nO(1)
            kl = kl + 1
 
            Dpp(ij,kl) = -(e(i,1) + e(j,1) - eF)*Kronecker_delta(i,k)&
                          *Kronecker_delta(j,l) + lambda*(ERI_aaaa(i,j,k,l) &
                          - ERI_aaaa(i,j,l,k))

          end  do
        end  do
      end  do
    end  do
  end if

  if (ispin == 7) then 

    ! bbbb block

    ij = 0
    do i=nC(2)+1,nO(2)
      do j=i+1,nO(2)
        ij = ij + 1
        kl = 0
        do k=nC(2)+1,nO(2)
          do l=k+1,nO(2)
            kl = kl + 1
 
           Dpp(ij,kl) = -(e(i,2) + e(j,2) - eF)*Kronecker_delta(i,k) &
                         *Kronecker_delta(j,l) + lambda*(ERI_bbbb(i,j,k,l) &
                         - ERI_bbbb(i,j,l,k))

          end  do
        end  do
      end  do
    end  do

  end if

end subroutine 
