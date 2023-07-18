subroutine MOM_overlap(nBas,nO,S,cG,c,ON)

! Compute overlap between old and new MO coefficients

  implicit none

! Input variables

  integer,intent(in)            :: nBas,nO
  double precision,intent(in)   :: S(nBas,nBas),cG(nBas,nBas),c(nBas,nBas)

! Local variables

  integer                       :: i,j,ploc
  double precision,allocatable  :: Ov(:,:),pOv(:)

! Output variables

  double precision,intent(inout):: ON(nBas)

  allocate(Ov(nBas,nBas),pOv(nBas))

  Ov = matmul(transpose(cG),matmul(S,c))

  pOv(:) = 0d0

  do i=1,nBas
    do j=1,nBas
      pOv(j) = pOv(j) + Ov(i,j)**2
    enddo
  enddo

  pOv(:) = sqrt(pOV(:))

! print*,'--- MOM overlap ---'
! call matout(nBas,1,pOv)

  ON(:) = 0d0

  do i=1,nO
    ploc = maxloc(pOv,nBas)
    ON(ploc) = 1d0
    pOv(ploc) = 0d0
 enddo


! print*,'--- Occupation numbers ---'
! call matout(nBas,1,ON)

  

end subroutine 
