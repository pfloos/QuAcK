subroutine AOtoMO_integral_transform(bra1,bra2,ket1,ket2,nBas,c,ERI_AO_basis,ERI_MO_basis)

! AO to MO transformation of two-electron integrals via the semi-direct O(N^5) algorithm
! bra and ket are the spin of (bra1 bra2|ket1 ket2)

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: bra1,bra2
  integer,intent(in)            :: ket1,ket2
  integer,intent(in)            :: nBas
  double precision,intent(in)   :: ERI_AO_basis(nBas,nBas,nBas,nBas),c(nBas,nBas,nspin)

! Local variables

  double precision,allocatable  :: scr(:,:,:,:)
  integer                       :: mu,nu,la,si,i,j,k,l

! Output variables

  double precision,intent(out)  :: ERI_MO_basis(nBas,nBas,nBas,nBas)

! Memory allocation

  allocate(scr(nBas,nBas,nBas,nBas))

! Four-index transform via semi-direct O(N^5) algorithm

  scr(:,:,:,:) = 0d0

  do l=1,nBas
    do si=1,nBas
      do la=1,nBas
        do nu=1,nBas
          do mu=1,nBas
            scr(mu,nu,la,l) = scr(mu,nu,la,l) + ERI_AO_basis(mu,nu,la,si)*c(si,l,ket2)
          enddo
        enddo
      enddo
    enddo
  enddo

  ERI_MO_basis(:,:,:,:) = 0d0

  do l=1,nBas
    do la=1,nBas
      do nu=1,nBas
        do i=1,nBas
          do mu=1,nBas
            ERI_MO_basis(i,nu,la,l) = ERI_MO_basis(i,nu,la,l) + c(mu,i,bra1)*scr(mu,nu,la,l)
          enddo
        enddo
      enddo
    enddo
  enddo

  scr(:,:,:,:) = 0d0

  do l=1,nBas 
    do k=1,nBas
      do la=1,nBas
        do nu=1,nBas
          do i=1,nBas
            scr(i,nu,k,l) = scr(i,nu,k,l) + ERI_MO_basis(i,nu,la,l)*c(la,k,bra2)
          enddo
        enddo
      enddo
    enddo
  enddo

  ERI_MO_basis(:,:,:,:) = 0d0

  do l=1,nBas
    do k=1,nBas
      do j=1,nBas
        do i=1,nBas
          do nu=1,nBas
            ERI_MO_basis(i,j,k,l) = ERI_MO_basis(i,j,k,l) + c(nu,j,ket1)*scr(i,nu,k,l)
          enddo
!         print*,i,k,j,l,ERI_MO_basis(i,j,k,l)
        enddo
      enddo
    enddo
  enddo


end subroutine AOtoMO_integral_transform
