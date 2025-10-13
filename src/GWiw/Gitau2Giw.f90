
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

  double precision              :: fact ! [to cancell large oscillations due to iw]

! Output variables

  complex*16,intent(out)        :: Giw1(nBas,nBas)
  complex*16,intent(out)        :: Giw2(nBas,nBas)
  
!------------------------------------------------------------------------
! Build G(i w) as the Fourier transform of  G(i tau)
!------------------------------------------------------------------------

  Giw1=czero; Giw2=czero;
  
  do itau=1,ntimes

   fact=1d0
   if(weval(1)>2d3) fact=0d0
   Giw1(:,:) = Giw1(:,:) - im*fact*tweight(itau)*Gitau(2*itau-1,:,:)*Exp( im*tcoord(itau)*weval(1)) &  ! G(i tau)
                         - im*fact*tweight(itau)*Gitau(2*itau  ,:,:)*Exp(-im*tcoord(itau)*weval(1))    ! G(-i tau)

   fact=1d0
   if(weval(2)>2d3) fact=0d0
   Giw2(:,:) = Giw2(:,:) - im*fact*tweight(itau)*Gitau(2*itau-1,:,:)*Exp( im*tcoord(itau)*weval(2)) &  ! G(i tau)
                         - im*fact*tweight(itau)*Gitau(2*itau  ,:,:)*Exp(-im*tcoord(itau)*weval(2))    ! G(-i tau)
  enddo

  Giw1=conjg(Giw1)
  Giw2=conjg(Giw2)

end subroutine

