subroutine ppULR_C(ispin,nBas,nC,nO,nV,nR,nPaa,nPab,nPbb,nPt,lambda,e,ERI_aaaa,ERI_aabb,ERI_bbbb,Cpp)

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
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: e(nBas,nspin)
  double precision,intent(in)   :: ERI_aaaa(nBas,nBas,nBas,nBas) 
  double precision,intent(in)   :: ERI_aabb(nBas,nBas,nBas,nBas) 
  double precision,intent(in)   :: ERI_bbbb(nBas,nBas,nBas,nBas) 
  
! Local variables
  
  double precision              :: eF 
  double precision,external     :: Kronecker_delta

  integer                       :: i,j,a,b,c,d,ab,cd

! Output variables

  double precision,intent(out)  :: Cpp(nPt,nPt)


  eF = 0d0 

  if(ispin == 3) then 

    ! abab block

    ab = 0
    do a=nO(1)+1,nBas-nR(1)
      do b=nO(2)+1,nBas-nR(2)
        ab = ab + 1
        cd = 0
        do c=nO(1)+1,nBas-nR(1)
          do d=nO(2)+1,nBas-nR(2)
            cd = cd + 1
            Cpp(ab,cd) = (e(a,1) + e(b,2))*Kronecker_delta(a,c) &
                          *Kronecker_delta(b,d) + lambda*ERI_aabb(a,b,c,d)
          end  do
        end  do
      end  do
    end  do

  end if

  if(ispin == 4) then  

    ! aaaa block

    ab = 0
    do a=nO(1)+1,nBas-nR(1)
      do b=a+1,nBas-nR(1)
        ab = ab + 1
        cd = 0
        do c=nO(1)+1,nBas-nR(1)
          do d=c+1,nBas-nR(1)
            cd = cd + 1
 
            Cpp(ab,cd) = (e(a,1) + e(b,1) - eF)*Kronecker_delta(a,c)&
                          *Kronecker_delta(b,d) + lambda*(ERI_aaaa(a,b,c,d) &
                          - ERI_aaaa(a,b,d,c))

          end  do
        end  do
      end  do
    end  do
  end if 


  if (ispin == 7) then

    ! bbbb block

    ab = 0
    do a=nO(2)+1,nBas-nR(2)
      do b=a+1,nBas-nR(2)
        ab = ab + 1
        cd = 0
        do c=nO(2)+1,nBas-nR(2)
          do d=c+1,nBas-nR(2)
            cd = cd + 1

           Cpp(ab,cd) = (e(a,2) + e(b,2) - eF)*Kronecker_delta(a,c) &
                         *Kronecker_delta(b,d) + lambda*(ERI_bbbb(a,b,c,d) &
                         - ERI_bbbb(a,b,d,c))

          end  do
        end  do
      end  do
    end  do

  end if

end subroutine 
