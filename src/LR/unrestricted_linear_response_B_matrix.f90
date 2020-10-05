subroutine unrestricted_linear_response_B_matrix(ispin,dRPA,nBas,nC,nO,nV,nR,nSa,nSb,nSt,lambda, & 
                                                 ERI_aaaa,ERI_aabb,ERI_bbbb,B_lr)

! Compute linear response

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dRPA
  integer,intent(in)            :: ispin
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC(nspin)
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nV(nspin)
  integer,intent(in)            :: nR(nspin)
  integer,intent(in)            :: nSa
  integer,intent(in)            :: nSb
  integer,intent(in)            :: nSt
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: ERI_aaaa(nBas,nBas,nBas,nBas) 
  double precision,intent(in)   :: ERI_aabb(nBas,nBas,nBas,nBas) 
  double precision,intent(in)   :: ERI_bbbb(nBas,nBas,nBas,nBas) 
  
! Local variables

  double precision              :: delta_dRPA
  double precision,external     :: Kronecker_delta

  integer                       :: i,j,a,b,ia,jb

! Output variables

  double precision,intent(out)  :: B_lr(nSt,nSt)

! Direct RPA

  delta_dRPA = 0d0
  if(dRPA) delta_dRPA = 1d0

!-----------------------------------------------
! Build B matrix for spin-conserving transitions
!-----------------------------------------------

  if(ispin == 1) then 

    ! aaaa block

    ia = 0
    do i=nC(1)+1,nO(1)
      do a=nO(1)+1,nBas-nR(1)
        ia = ia + 1
        jb = 0
        do j=nC(1)+1,nO(1)
          do b=nO(1)+1,nBas-nR(1)
            jb = jb + 1
 
            B_lr(ia,jb) = lambda*ERI_aaaa(i,j,a,b) - (1d0 - delta_dRPA)*lambda*ERI_aaaa(i,j,b,a)

          end  do
        end  do
      end  do
    end  do

    ! aabb block

    ia = 0
    do i=nC(1)+1,nO(1)
      do a=nO(1)+1,nBas-nR(1)
        ia = ia + 1
        jb = 0
        do j=nC(2)+1,nO(2)
          do b=nO(2)+1,nBas-nR(2)
            jb = jb + 1
 
            B_lr(ia,nSa+jb) = lambda*ERI_aabb(i,j,a,b) 

          end  do
        end  do
      end  do
    end  do

    ! bbaa block

    ia = 0
    do i=nC(2)+1,nO(2)
      do a=nO(2)+1,nBas-nR(2)
        ia = ia + 1
        jb = 0
        do j=nC(1)+1,nO(1)
          do b=nO(1)+1,nBas-nR(1)
            jb = jb + 1
 
            B_lr(nSa+ia,jb) = lambda*ERI_aabb(j,i,b,a)

          end  do
        end  do
      end  do
    end  do

    ! bbbb block

    ia = 0
    do i=nC(2)+1,nO(2)
      do a=nO(2)+1,nBas-nR(2)
        ia = ia + 1
        jb = 0
        do j=nC(2)+1,nO(2)
          do b=nO(2)+1,nBas-nR(2)
            jb = jb + 1
 
            B_lr(nSa+ia,nSa+jb) = lambda*ERI_bbbb(i,j,a,b) - (1d0 - delta_dRPA)*lambda*ERI_bbbb(i,j,b,a)

          end  do
        end  do
      end  do
    end  do

  end if

!-----------------------------------------------
! Build B matrix for spin-flip transitions
!-----------------------------------------------

  if(ispin == 2) then

    B_lr(:,:) = 0d0

    ! abba block

    ia = 0
    do i=nC(1)+1,nO(1)
      do a=nO(2)+1,nBas-nR(2)
        ia = ia + 1
        jb = 0
        do j=nC(2)+1,nO(2)
          do b=nO(1)+1,nBas-nR(1)
            jb = jb + 1
 
            B_lr(ia,nSa+jb) = - (1d0 - delta_dRPA)*lambda*ERI_aabb(i,j,b,a)

          end  do
        end  do
      end  do
    end  do

    ! baab block

    ia = 0
    do i=nC(2)+1,nO(2)
      do a=nO(1)+1,nBas-nR(1)
        ia = ia + 1
        jb = 0
        do j=nC(1)+1,nO(1)
          do b=nO(2)+1,nBas-nR(2)
            jb = jb + 1
 
            B_lr(nSa+ia,jb) = - (1d0 - delta_dRPA)*lambda*ERI_aabb(j,i,a,b)

          end  do
        end  do
      end  do
    end  do
  
  end if


end subroutine unrestricted_linear_response_B_matrix
