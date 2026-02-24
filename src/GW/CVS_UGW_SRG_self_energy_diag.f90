subroutine CVS_UGW_SRG_self_energy_diag(flow,nBas,nC,nO,nV,nR,nSt,nCVS,nFC,occupations,virtuals,e,Om,rho,EcGM,Sig,Z)

! Compute diagonal of the correlation part of the self-energy

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: flow
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC(nspin)
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nV(nspin)
  integer,intent(in)            :: nR(nspin)
  integer,intent(in)            :: nSt
  double precision,intent(in)   :: e(nBas,nspin)
  double precision,intent(in)   :: Om(nSt)
  double precision,intent(in)   :: rho(nBas,nBas,nSt,nspin)
  
  integer,intent(in)            :: nCVS(nspin),nFC(nspin)
  integer,intent(in)            :: occupations(maxval(nO-nFC),nspin)
  integer,intent(in)            :: virtuals(nBas-minval(nO),nspin)

! Local variables

  integer                       :: i,a,p,m
  double precision              :: num,eps
  double precision              :: s
  double precision              :: Dpim,Dpam,Diam

! Output variables

  double precision,intent(out)  :: Sig(nBas,nspin)
  double precision,intent(out)  :: Z(nBas,nspin)
  double precision              :: EcGM(nspin)

! Initialize 

  Sig(:,:) = 0d0
  Z(:,:)   = 0d0
  EcGM(:)  = 0d0
  s = flow
  
  print *, "SRG not implemented for CVS/MOM yet"

end subroutine 
