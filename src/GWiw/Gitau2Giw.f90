
subroutine Gitau2Giw(nBas,ntimes,ntimes_twice,tweight,tcoord,weval,Gitau,Giw1,Giw2)

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: ntimes
  integer,intent(in)            :: ntimes_twice

  double precision,intent(in)   :: weval(2)
  double precision,intent(in)   :: tweight(ntimes)
  double precision,intent(in)   :: tcoord(ntimes)

  complex*16,intent(in)         :: Gitau(ntimes_twice,nBas,nBas)

! Local variables

  integer                       :: itau

! Output variables

  complex*16,intent(out)        :: Giw1(nBas,nBas)
  complex*16,intent(out)        :: Giw2(nBas,nBas)
  
!------------------------------------------------------------------------
! Build G(i w) as the Fourier transform of  G(i tau)
!------------------------------------------------------------------------

  Giw1=czero; Giw2=czero;
  
  do itau=1,ntimes

   Giw1(:,:) = Giw1(:,:) - im*tweight(itau)*Gitau(2*itau-1,:,:)*Exp( im*tcoord(itau)*weval(1)) &  ! G(i tau)
                         - im*tweight(itau)*Gitau(2*itau  ,:,:)*Exp(-im*tcoord(itau)*weval(1))    ! G(-i tau)

   Giw2(:,:) = Giw2(:,:) - im*tweight(itau)*Gitau(2*itau-1,:,:)*Exp( im*tcoord(itau)*weval(2)) &  ! G(i tau)
                         - im*tweight(itau)*Gitau(2*itau  ,:,:)*Exp(-im*tcoord(itau)*weval(2))    ! G(-i tau)
  enddo

  Giw1=conjg(Giw1)
  Giw2=conjg(Giw2)

end subroutine

