subroutine AOtoMO_ERI_GHF(nBas,nBas2,c1,c2,ERI_AO_basis,ERI_MO_basis)

! AO to MO transformation of two-electron integrals via the semi-direct O(N^5) algorithm
! bra and ket are the spin of (bra1 bra2|ket1 ket2)

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nBas2
  double precision,intent(in)   :: ERI_AO_basis(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: c1(nBas,nBas2)
  double precision,intent(in)   :: c2(nBas,nBas2)

! Local variables

  double precision,allocatable  :: scr(:,:,:,:)
  integer                       :: mu,nu,la,si,i,j,k,l

! Output variables

  double precision,intent(out)  :: ERI_MO_basis(nBas2,nBas2,nBas2,nBas2)

! Memory allocation

  allocate(scr(nBas2,nBas2,nBas2,nBas2))

! Four-index transform via semi-direct O(N^5) algorithm

  scr(:,:,:,:) = 0d0

  do l=1,nBas2
    do si=1,nBas
      do la=1,nBas
        do nu=1,nBas
          do mu=1,nBas
            scr(mu,nu,la,l) = scr(mu,nu,la,l) + ERI_AO_basis(mu,nu,la,si)*c2(si,l)
          enddo
        enddo
      enddo
    enddo
  enddo

  ERI_MO_basis(:,:,:,:) = 0d0

  do l=1,nBas2
    do la=1,nBas
      do nu=1,nBas
        do i=1,nBas2
          do mu=1,nBas
            ERI_MO_basis(i,nu,la,l) = ERI_MO_basis(i,nu,la,l) + c1(mu,i)*scr(mu,nu,la,l)
          enddo
        enddo
      enddo
    enddo
  enddo

  scr(:,:,:,:) = 0d0

  do l=1,nBas2
    do k=1,nBas2
      do la=1,nBas
        do nu=1,nBas
          do i=1,nBas2
            scr(i,nu,k,l) = scr(i,nu,k,l) + ERI_MO_basis(i,nu,la,l)*c1(la,k)
          enddo
        enddo
      enddo
    enddo
  enddo

  ERI_MO_basis(:,:,:,:) = 0d0

  do l=1,nBas2
    do k=1,nBas2
      do j=1,nBas2
        do i=1,nBas2
          do nu=1,nBas
            ERI_MO_basis(i,j,k,l) = ERI_MO_basis(i,j,k,l) + c2(nu,j)*scr(i,nu,k,l)
          enddo
        enddo
      enddo
    enddo
  enddo

end subroutine 
