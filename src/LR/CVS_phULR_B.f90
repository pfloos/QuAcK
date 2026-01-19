subroutine CVS_phULR_B(ispin,dRPA,nBas,nC,nO,nV,nR,nSa,nSb,nSt,nCVS,occupations,virtuals,lambda,ERI_aaaa,ERI_aabb,ERI_bbbb,Bph)

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
  integer,intent(in)            :: nCVS(nspin)
  integer,intent(in)            :: nSa
  integer,intent(in)            :: nSb
  integer,intent(in)            :: nSt
  integer,intent(in)            :: occupations(maxval(nO),nspin)
  integer,intent(in)            :: virtuals(nBas- minval(nO),nspin)
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: ERI_aaaa(nBas,nBas,nBas,nBas) 
  double precision,intent(in)   :: ERI_aabb(nBas,nBas,nBas,nBas) 
  double precision,intent(in)   :: ERI_bbbb(nBas,nBas,nBas,nBas) 
  
! Local variables

  double precision              :: delta_dRPA
  double precision,external     :: Kronecker_delta

  integer                       :: i,j,a,b,ia,jb

! Output variables

  double precision,intent(out)  :: Bph(nSt,nSt)

! Direct RPA

  delta_dRPA = 0d0
  if(dRPA) delta_dRPA = 1d0

! Initialization

  Bph(:,:) = 0d0

!-----------------------------------------------
! Build B matrix for spin-conserving transitions
!-----------------------------------------------

  if(ispin == 1) then 

    ! aaaa block

    ia = 0
    do i=1,nO(1)
      do a=nCVS(1)+1,nBas-nO(1)
        ia = ia + 1
        jb = 0
        do j=1,nO(1)
          do b=nCVS(1)+1,nBas-nO(1)
            jb = jb + 1
 
            Bph(ia,jb) = lambda*ERI_aaaa(occupations(i,1),occupations(j,1),virtuals(a,1),virtuals(b,1))&
                       - (1d0 - delta_dRPA)*lambda*ERI_aaaa(occupations(i,1),occupations(j,1),virtuals(b,1),virtuals(a,1))

           end do
         end do
       end do
     end do

    ! aabb block

    ia = 0
    do i=1,nO(1)
      do a=nCVS(1)+1,nBas-nO(1)
        ia = ia + 1
        jb = 0
        do j=1,nO(2)
          do b=nCVS(2)+1,nBas-nO(2)
            jb = jb + 1
 
            Bph(ia,nSa+jb) = lambda*ERI_aabb(occupations(i,1),occupations(j,2),virtuals(a,1),virtuals(b,2)) 

          end do
        end do
      end do
    end do

    ! bbaa block

    ia = 0
    do i=1,nO(2)
      do a=nCVS(2)+1,nBas-nO(2)
        ia = ia + 1
        jb = 0
        do j=1,nO(1)
          do b=nCVS(1)+1,nBas-nO(1)
            jb = jb + 1
 
            Bph(nSa+ia,jb) = lambda*ERI_aabb(occupations(j,1),occupations(i,2),virtuals(b,1),virtuals(a,2))

          end do
        end do
      end do
    end do

    ! bbbb block

    ia = 0
    do i=1,nO(2)
      do a=nCVS(2)+1,nBas-nO(2)
        ia = ia + 1
        jb = 0
        do j=1,nO(2)
          do b=nCVS(2)+1,nBas-nO(2)
            jb = jb + 1
 
            Bph(nSa+ia,nSa+jb) = lambda*ERI_bbbb(occupations(i,2),occupations(j,2),virtuals(a,2),virtuals(b,2))&
                               - (1d0 - delta_dRPA)*lambda*ERI_bbbb(occupations(i,2),occupations(j,2),virtuals(b,2),virtuals(a,2))

         end do
       end do
     end do
   end do

  end if

!-----------------------------------------------
! Build B matrix for spin-flip transitions
!-----------------------------------------------

  if(ispin == 2) then

    ! abba block

    ia = 0
    do i=1,nO(1)
      do a=nCVS(2)+1,nBas-nO(2)
        ia = ia + 1
        jb = 0
        do j=1,nO(2)
          do b=1,nBas-nO(1)
            jb = jb + 1
 
            Bph(ia,nSa+jb) = - (1d0 - delta_dRPA)*lambda*ERI_aabb(occupations(i,1),occupations(j,2),virtuals(b,1),virtuals(a,2))

          end do
        end do
      end do
    end do

    ! baab block

    ia = 0
    do i=1,nO(2)
      do a=nCVS(1)+1,nBas-nO(1)
        ia = ia + 1
        jb = 0
        do j=1,nO(1)
          do b=nCVS(2)+1,nBas-nO(2)
            jb = jb + 1
 
            Bph(nSa+ia,jb) = - (1d0 - delta_dRPA)*lambda*ERI_aabb(occupations(j,1),occupations(i,2),virtuals(a,1),virtuals(b,2))

          end do
        end do
      end do
    end do
  
  end if


end subroutine 
