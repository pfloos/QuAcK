subroutine unrestricted_linear_response_B_pp(ispin,nBas,nC,nO,nV,nR,nPaa,nPab,nPbb,nPt,nHaa, &  
 nHab,nHbb,nHt,lambda,ERI_aaaa,ERI_aabb,ERI_bbbb,B_pp)

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
  integer,intent(in)            :: nPaa
  integer,intent(in)            :: nPab
  integer,intent(in)            :: nPbb
  integer,intent(in)            :: nPt 
  integer,intent(in)            :: nHaa
  integer,intent(in)            :: nHab
  integer,intent(in)            :: nHbb
  integer,intent(in)            :: nHt
  double precision,intent(in)   :: lambda 
  double precision,intent(in)   :: ERI_aaaa(nBas,nBas,nBas,nBas) 
  double precision,intent(in)   :: ERI_aabb(nBas,nBas,nBas,nBas) 
  double precision,intent(in)   :: ERI_bbbb(nBas,nBas,nBas,nBas) 
  
! Local variables
  
  double precision              :: eF 
  double precision,external     :: Kronecker_delta

  integer                       :: i,j,a,b,c,d,ij,ab

! Output variables

  double precision,intent(out)  :: B_pp(nPt,nHt)


  eF = 0d0 

!-----------------------------------------------
! Build B matrix for spin-conserved transitions
!-----------------------------------------------

  if(ispin == 1) then

    ! abab block

    ab = 0
    do a=nO(1)+1,nBas-nR(1)
      do b=nO(2)+1,nBas-nR(2)
        ab = ab + 1
        ij = 0
        do i=nC(1)+1,nO(1)
          do j=nC(2)+1,nO(2)
            ij = ij + 1

            B_pp(ab,ij) = lambda*ERI_aabb(a,b,i,j)

          end  do
        end  do
      end  do
    end  do

  end if

!-----------------------------------------------
! Build B matrix for spin-flip transitions
!-----------------------------------------------

  if(ispin == 2) then  

    ! aaaa block

    ab = 0
    do a=nO(1)+1,nBas-nR(1)
      do b=a+1,nBas-nR(1)
        ab = ab + 1
        ij = 0
        do i=nC(1)+1,nO(1)
          do j=i+1,nO(1)
            ij = ij + 1
 
            B_pp(ab,ij) = lambda*(ERI_aaaa(a,b,i,j) - ERI_aaaa(a,b,j,i))
            
          end  do
        end  do
      end  do
    end  do
 end if 

  if (ispin == 3) then

    ! bbbb block

    ab = 0
    do a=nO(2)+1,nBas-nR(2)
      do b=a+1,nBas-nR(2)
        ab = ab + 1
        ij = 0
        do i=nC(2)+1,nO(2)
          do j=i+1,nO(2)
            ij = ij + 1
 
            B_pp(ab,ij) = lambda*(ERI_bbbb(a,b,i,j) - ERI_bbbb(a,b,j,i)) 

          end  do
        end  do
      end  do
    end  do

  end if

end subroutine unrestricted_linear_response_B_pp
