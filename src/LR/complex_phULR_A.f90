subroutine complex_phULR_A(ispin,dRPA,nBas,nC,nO,nV,nR,nSa,nSb,nSt,nCVS,nFC,occupations,virtuals,lambda,eHF,ERI_aaaa,ERI_aabb,ERI_bbbb,Aph)

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
  integer,intent(in)            :: nCVS(nspin),nFC(nspin)
  integer,intent(in)            :: occupations(maxval(nO-nFC),nspin)
  integer,intent(in)            :: virtuals(nBas-minval(nO),nspin)
  integer,intent(in)            :: nSa
  integer,intent(in)            :: nSb
  integer,intent(in)            :: nSt
  double precision,intent(in)   :: lambda
  complex*16,intent(in)         :: eHF(nBas,nspin)
  complex*16,intent(in)         :: ERI_aaaa(nBas,nBas,nBas,nBas) 
  complex*16,intent(in)         :: ERI_aabb(nBas,nBas,nBas,nBas) 
  complex*16,intent(in)         :: ERI_bbbb(nBas,nBas,nBas,nBas) 
  
! Local variables

  double precision              :: delta_dRPA
  double precision,external     :: Kronecker_delta

  integer                       :: i,j,a,b,ia,jb

! Output variables

  complex*16,intent(out)        :: Aph(nSt,nSt)

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
    do i=1,nO(1) - nFC(1)
      do a=nCVS(1)+1,nBas - nO(1)
        ia = ia + 1
        jb = 0
        do j=1,nO(1) - nFC(1)
          do b=nCVS(1)+1,nBas - nO(1)
            jb = jb + 1
            
            Aph(ia,jb) = (eHF(virtuals(a,1),1) - eHF(occupations(i,1),1))*Kronecker_delta(i,j)*Kronecker_delta(a,b) &
                       + lambda*ERI_aaaa(occupations(i,1),virtuals(b,1),virtuals(a,1),occupations(j,1)) &
                       - (1d0 - delta_dRPA)*lambda*ERI_aaaa(occupations(i,1),virtuals(b,1),occupations(j,1),virtuals(a,1))

          end do
        end do
      end do
    end do

    ! aabb block

    ia = 0
    do i=1,nO(1) - nFC(1)
      do a=nCVS(1)+1,nBas - nO(1)
        ia = ia + 1
        jb = 0
        do j=1,nO(2) - nFC(2)
          do b=nCVS(2)+1,nBas - nO(2)
            jb = jb + 1
            Aph(ia,nSa+jb) = lambda*ERI_aabb(occupations(i,1),virtuals(b,2),virtuals(a,1),occupations(j,2)) 

          end do
        end do
      end do
    end do

    ! bbaa block

    ia = 0
    do i=1,nO(2) - nFC(2)
      do a=nCVS(2)+1,nBas-nO(2)
        ia = ia + 1
        jb = 0
        do j=1,nO(1) - nFC(1)
          do b=nCVS(1)+1,nBas-nO(1)
            jb = jb + 1
 
            Aph(nSa+ia,jb) = lambda*ERI_aabb(virtuals(b,1),occupations(i,2),occupations(j,1),virtuals(a,2))

          end do
        end do
      end do
    end do

    ! bbbb block

    ia = 0
    do i=1,nO(2) - nFC(2)
      do a=nCVS(2)+1,nBas-nO(2)
        ia = ia + 1
        jb = 0
        do j=1,nO(2) - nFC(2)
          do b=nCVS(2)+1,nBas-nO(2)
            jb = jb + 1
 
            Aph(nSa+ia,nSa+jb) = (eHF(virtuals(a,2),2) - eHF(occupations(i,2),2))*Kronecker_delta(i,j)*Kronecker_delta(a,b) &
                               + lambda*ERI_bbbb(occupations(i,2),virtuals(b,2),virtuals(a,2),occupations(j,2))&
                               - (1d0 - delta_dRPA)*lambda*ERI_bbbb(occupations(i,2),virtuals(b,2),occupations(j,2),virtuals(a,2))

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
    do i=1,nO(1) -nFC(1)
      do a=nCVS(2)+1,nBas-nO(2)
        ia = ia + 1
        jb = 0
        do j=1,nO(1)-nFC(1)
          do b=nCVS(2)+1,nBas-nO(2)
            jb = jb + 1
            Aph(ia,jb) = (eHF(virtuals(a,2),2) - eHF(occupations(i,1),1))*Kronecker_delta(i,j)*Kronecker_delta(a,b) &
                       - (1d0 - delta_dRPA)*lambda*ERI_aabb(occupations(i,1),virtuals(b,2),occupations(j,1),virtuals(a,2))
           end do
         end do
       end do
     end do

    ! baba block

    ia = 0
    do i=1,nO(2) - nFC(2)
      do a=nCVS(1)+1,nBas-nO(1)
        ia = ia + 1
        jb = 0
        do j=1,nO(2) - nFC(2)
          do b=nCVS(1)+1,nBas-nO(1)
            jb = jb + 1
 
            Aph(nSa+ia,nSa+jb) = (eHF(virtuals(a,1),1) - eHF(occupations(i,2),2))*Kronecker_delta(i,j)*Kronecker_delta(a,b) &
                               - (1d0 - delta_dRPA)*lambda*ERI_aabb(virtuals(b,1),occupations(i,2),virtuals(a,1),occupations(j,2))

           end do
         end do
       end do
     end do

  end if

end subroutine 
