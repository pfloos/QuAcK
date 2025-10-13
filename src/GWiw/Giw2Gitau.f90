
subroutine Giw2Gitau(nBas,nfreqs,wweight,wcoord,tcoord,Giw,Gitau_plus,Gitau_minus)

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nfreqs

  double precision,intent(in)   :: tcoord
  double precision,intent(in)   :: wweight(nfreqs)
  double precision,intent(in)   :: wcoord(nfreqs)

  complex*16,intent(in)         :: Giw(nfreqs,nBas,nBas)

! Local variables

  integer                       :: ifreq
  integer                       :: ibas
  double precision              :: teval(2)

  double precision              :: fact ! [to cancell large oscillations due to itau or iw]

! Output variables

  complex*16,intent(out)        :: Gitau_plus(nBas,nBas)
  complex*16,intent(out)        :: Gitau_minus(nBas,nBas)
  
!------------------------------------------------------------------------
! Build G(i w) as the Fourier transform of  G(i tau)
!------------------------------------------------------------------------

  Gitau_plus=czero; Gitau_minus=czero;
  teval(1)=  tcoord  
  teval(2)= -tcoord  
 
  do ifreq=1,nfreqs
    
   fact=1d0
!   if(wcoord(ifreq)>1d2 .or. tcoord>1d2) fact=0d0

   Gitau_plus(:,:)  = Gitau_plus(:,:)                                                              &
                    + im*fact*wweight(ifreq)*      Giw(ifreq,:,:) *Exp( im*wcoord(ifreq)*teval(1)) &  ! G(i w)
                    + im*fact*wweight(ifreq)*conjg(Giw(ifreq,:,:))*Exp(-im*wcoord(ifreq)*teval(1))    ! G(-i w)

   Gitau_minus(:,:) = Gitau_minus(:,:)                                                             &
                    + im*fact*wweight(ifreq)*      Giw(ifreq,:,:) *Exp( im*wcoord(ifreq)*teval(2)) &  ! G(i w)
                    + im*fact*wweight(ifreq)*conjg(Giw(ifreq,:,:))*Exp(-im*wcoord(ifreq)*teval(2))    ! G(-i tau)
  enddo

  Gitau_plus=Gitau_plus/(2d0*pi)
  Gitau_minus=Gitau_minus/(2d0*pi)

end subroutine

