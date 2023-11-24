subroutine AOtoMO_ERI_GHF(nBas,nBas2,c1,c2,ERI_AO,ERI_MO)

! AO to MO transformation of two-electron integrals via the semi-direct O(N^5) algorithm
! bra and ket are the spin of (bra1 bra2|ket1 ket2)

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nBas2
  double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: c1(nBas,nBas2)
  double precision,intent(in)   :: c2(nBas,nBas2)

! Local variables

  double precision,allocatable  :: scr(:,:,:,:)
  integer                       :: mu,nu,la,si
  integer                       :: i,j,k,l

! Output variables

  double precision,intent(out)  :: ERI_MO(nBas2,nBas2,nBas2,nBas2)

! Memory allocation

  allocate(scr(nBas2,nBas2,nBas2,nBas2))

  call dgemm('T','N',nBas**3,nBas2,nBas,1d0,ERI_AO,nBas,c1(1,1),nBas,0d0,scr,nBas**3)

  call dgemm('T','N',nBas**2*nBas2,nBas2,nBas,1d0,scr,nBas,c2(1,1),nBas,0d0,ERI_MO,nBas**2*nBas2)

  call dgemm('T','N',nBas*nBas2**2,nBas2,nBas,1d0,ERI_MO,nBas,c1(1,1),nBas,0d0,scr,nBas*nBas2**2)

  call dgemm('T','N',nBas2**3,nBas2,nBas,1d0,scr,nBas,c2(1,1),nBas,0d0,ERI_MO,nBas2**3)

! Four-index transform via semi-direct O(N^5) algorithm

! scr(:,:,:,:) = 0d0

! do l=1,nBas2
!   do si=1,nBas
!     do la=1,nBas
!       do nu=1,nBas
!         do mu=1,nBas
!           scr(mu,nu,la,l) = scr(mu,nu,la,l) + ERI_AO(mu,nu,la,si)*c2(si,l)
!         enddo
!       enddo
!     enddo
!   enddo
! enddo

! ERI_MO(:,:,:,:) = 0d0

! do l=1,nBas2
!   do la=1,nBas
!     do nu=1,nBas
!       do i=1,nBas2
!         do mu=1,nBas
!           ERI_MO(i,nu,la,l) = ERI_MO(i,nu,la,l) + c1(mu,i)*scr(mu,nu,la,l)
!         enddo
!       enddo
!     enddo
!   enddo
! enddo


! scr(:,:,:,:) = 0d0

! do l=1,nBas2
!   do k=1,nBas2
!     do la=1,nBas
!       do nu=1,nBas
!         do i=1,nBas2
!           scr(i,nu,k,l) = scr(i,nu,k,l) + ERI_MO(i,nu,la,l)*c1(la,k)
!         enddo
!       enddo
!     enddo
!   enddo
! enddo

! ERI_MO(:,:,:,:) = 0d0

! do l=1,nBas2
!   do k=1,nBas2
!     do j=1,nBas2
!       do i=1,nBas2
!         do nu=1,nBas
!           ERI_MO(i,j,k,l) = ERI_MO(i,j,k,l) + c2(nu,j)*scr(i,nu,k,l)
!         enddo
!       enddo
!     enddo
!   enddo
! enddo

end subroutine 
