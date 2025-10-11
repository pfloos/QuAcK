
subroutine Gitau2Giw(nBas,ntimes,tweight,tcoord,weval,Gitau,Giw)

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: ntimes

  double precision,intent(in)   :: weval
  double precision,intent(in)   :: tweight(ntimes)
  double precision,intent(in)   :: tcoord(ntimes)

  complex*16,intent(in)         :: Gitau(ntimes,nBas,nBas)

! Local variables

  integer                       :: itau

! Output variables

  complex*16,intent(out)        :: Giw(nBas,nBas)
  
!------------------------------------------------------------------------
! Build G(i w) as the Fourier transform of  G(i tau)
!------------------------------------------------------------------------

  Giw=czero
  
  do itau=1,ntimes
   Giw(:,:) = Giw(:,:) - im*tweight(itau)*Gitau(2*itau-1,:,:)*Exp( im*tcoord(itau)*weval) &  ! G(i tau)
                       - im*tweight(itau)*Gitau(2*itau  ,:,:)*Exp(-im*tcoord(itau)*weval)    ! G(-i tau)
  enddo

end subroutine

