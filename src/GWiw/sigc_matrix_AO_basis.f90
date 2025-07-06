subroutine sigc_AO_basis(nBas,nOrb,nOrb_twice,c,U_QP,eqsGWB_state,vMAT,nfreqs,ntimes,wcoord,wweight, &
                        Sigc_he,Sigc_hh)

! Compute Sigma_c matrix in the AO basis

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nOrb_twice
  integer,intent(in)            :: nfreqs
  integer,intent(in)            :: ntimes

  double precision,intent(in)   :: wcoord(nfreqs),wweight(nfreqs)
  double precision,intent(in)   :: U_QP(nOrb_twice,nOrb_twice)
  double precision,intent(in)   :: eqsGWB_state(nOrb_twice)
  double precision,intent(in)   :: vMAT(nOrb*nOrb,nOrb*nOrb)
  double precision,intent(in)   :: c(nBas,nOrb)

! Local variables

! Output variables

  double precision,intent(out)  :: Sigc_he(nBas,nBas)
  double precision,intent(out)  :: Sigc_hh(nBas,nBas)

!


 Sigc_he=0d0
 Sigc_hh=0d0

end subroutine 

! ---

