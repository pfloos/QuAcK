subroutine CVS_phULR_A(ispin,dRPA,nBas,nC,nO,nV,nR,nSa,nSb,nSt,nCVS,occupations,virtuals,lambda,eHF,ERI_aaaa,ERI_aabb,ERI_bbbb,Aph)

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
  integer,intent(in)            :: occupations(maxval(nO),nspin)
  integer,intent(in)            :: virtuals(nBas-minval(nO),nspin)
  integer,intent(in)            :: nSa
  integer,intent(in)            :: nSb
  integer,intent(in)            :: nSt
  integer,intent(in)            :: nCVS(nspin)
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: eHF(nBas,nspin)
  double precision,intent(in)   :: ERI_aaaa(nBas,nBas,nBas,nBas) 
  double precision,intent(in)   :: ERI_aabb(nBas,nBas,nBas,nBas) 
  double precision,intent(in)   :: ERI_bbbb(nBas,nBas,nBas,nBas) 
  
! Local variables

  double precision              :: delta_dRPA
  double precision,external     :: Kronecker_delta

  integer                       :: i,j,a,b,ia,jb

! Output variables

  double precision,intent(out)  :: Aph(nSt,nSt)

! Direct RPA

  delta_dRPA = 0d0
  if(dRPA) delta_dRPA = 1d0

! Initialization

  Aph(:,:) = 0d0

!----------------------------------------------
! Build A matrix for spin-conserved transitions
!----------------------------------------------

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
 
            Aph(ia,jb) = (eHF(a,1) - eHF(i,1))*Kronecker_delta(i,j)*Kronecker_delta(a,b) &
                       + lambda*ERI_aaaa(i,b,a,j) - (1d0 - delta_dRPA)*lambda*ERI_aaaa(i,b,j,a)

          end do
        end do
      end do
    end do

    ! aabb block

    ia = 0
    do i=nC(1)+1,nO(1)
      do a=nO(1)+1,nBas-nR(1)
        ia = ia + 1
        jb = 0
        do j=nC(2)+1,nO(2)
          do b=nO(2)+1,nBas-nR(2)
            jb = jb + 1
 
            Aph(ia,nSa+jb) = lambda*ERI_aabb(i,b,a,j) 

          end do
        end do
      end do
    end do

    ! bbaa block

    ia = 0
    do i=nC(2)+1,nO(2)
      do a=nO(2)+1,nBas-nR(2)
        ia = ia + 1
        jb = 0
        do j=nC(1)+1,nO(1)
          do b=nO(1)+1,nBas-nR(1)
            jb = jb + 1
 
            Aph(nSa+ia,jb) = lambda*ERI_aabb(b,i,j,a)

          end do
        end do
      end do
    end do

    ! bbbb block

    ia = 0
    do i=nC(2)+1,nO(2)
      do a=nO(2)+1,nBas-nR(2)
        ia = ia + 1
        jb = 0
        do j=nC(2)+1,nO(2)
          do b=nO(2)+1,nBas-nR(2)
            jb = jb + 1
 
            Aph(nSa+ia,nSa+jb) = (eHF(a,2) - eHF(i,2))*Kronecker_delta(i,j)*Kronecker_delta(a,b) &
                               + lambda*ERI_bbbb(i,b,a,j) - (1d0 - delta_dRPA)*lambda*ERI_bbbb(i,b,j,a)

          end do
        end do
      end do
    end do

  end if

!-----------------------------------------
! Build A matrix for spin-flip transitions
!-----------------------------------------

  if(ispin == 2) then 

    ! abab block

    ia = 0
    do i=nC(1)+1,nO(1)
      do a=nO(2)+1,nBas-nR(2)
        ia = ia + 1
        jb = 0
        do j=nC(1)+1,nO(1)
          do b=nO(2)+1,nBas-nR(2)
            jb = jb + 1
            Aph(ia,jb) = (eHF(a,2) - eHF(i,1))*Kronecker_delta(i,j)*Kronecker_delta(a,b) &
                       - (1d0 - delta_dRPA)*lambda*ERI_aabb(i,b,j,a)
           end do
         end do
       end do
     end do

    ! baba block

    ia = 0
    do i=nC(2)+1,nO(2)
      do a=nO(1)+1,nBas-nR(1)
        ia = ia + 1
        jb = 0
        do j=nC(2)+1,nO(2)
          do b=nO(1)+1,nBas-nR(1)
            jb = jb + 1
 
            Aph(nSa+ia,nSa+jb) = (eHF(a,1) - eHF(i,2))*Kronecker_delta(i,j)*Kronecker_delta(a,b) &
                               - (1d0 - delta_dRPA)*lambda*ERI_aabb(b,i,a,j)

           end do
         end do
       end do
     end do

  end if

end subroutine 
