subroutine sigc_AO_basis_RHF(nBas,nOrb,nO,eta,shift,c,eqsGW_state,S,vMAT,nfreqs,ntimes, &
                             wcoord,wweight,Sigc_ao)

! Compute Sigma_c matrix in the AO basis

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nO
  integer,intent(in)            :: nfreqs
  integer,intent(in)            :: ntimes

  double precision,intent(in)   :: eta
  double precision,intent(in)   :: shift
  double precision,intent(in)   :: wcoord(nfreqs),wweight(nfreqs)
  double precision,intent(in)   :: eqsGW_state(nOrb)
  double precision,intent(in)   :: vMAT(nOrb*nOrb,nOrb*nOrb)
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: c(nBas,nOrb)

! Local variables

  integer                       :: iorb,jorb
  integer                       :: iE
  integer                       :: nE_eval_global

  double precision              :: chem_pot
  double precision,allocatable  :: Sigc_mo_tmp(:,:,:)
  double precision,allocatable  :: Sigc_mo(:,:)

  complex *16,allocatable       :: Sigc_mo_cpx(:,:,:)
  complex *16,allocatable       :: E_eval_global_cpx(:)

! Output variables

  double precision,intent(out)  :: Sigc_ao(nBas,nBas)

! Initialize variables
    
  chem_pot = 0.5d0*(eqsGW_state(nO)+eqsGW_state(nO+1))

! Se energies using cluster method or just using two shifts
  if(.false.) then
   allocate(E_eval_global_cpx(1))
   call set_Eeval_cluster(nOrb,nOrb,1,shift,eqsGW_state,nE_eval_global,&
                          E_eval_global_cpx,chem_pot)
   deallocate(E_eval_global_cpx)
   allocate(E_eval_global_cpx(nE_eval_global))
   call set_Eeval_cluster(nOrb,nOrb,nE_eval_global,shift,eqsGW_state,&
                           nE_eval_global,E_eval_global_cpx,chem_pot)
  else
   nE_eval_global=2*nOrb
   allocate(E_eval_global_cpx(nE_eval_global))
   do iorb=1,nOrb
    E_eval_global_cpx(2*iorb-1)=eqsGW_state(iorb)-shift-chem_pot
    E_eval_global_cpx(2*iorb)  =eqsGW_state(iorb)+shift-chem_pot
   enddo
  endif

! Run over energies
  allocate(Sigc_mo_cpx(nE_eval_global,nOrb,nOrb))
  call build_Sigmac_w_RHF(nOrb,nO,nE_eval_global,eta,0,E_eval_global_cpx,eqsGW_state,nfreqs,ntimes,&
                          wweight,wcoord,vMAT,Sigc_mo_cpx)
  deallocate(E_eval_global_cpx)

! Interpolate and transform Sigma from MO to AO basis
! TODO interpolate Sigma with clusters method
  allocate(Sigc_mo_tmp(nOrb,nOrb,nOrb))
  allocate(Sigc_mo(nOrb,nOrb))
  do iorb=1,nOrb
   Sigc_mo_tmp(iorb,:,:)=0.5d0*(Real(Sigc_mo_cpx(2*iorb-1,:,:))+Real(Sigc_mo_cpx(2*iorb,:,:)))
  enddo
  deallocate(Sigc_mo_cpx)
  do iorb=1,nOrb
   Sigc_mo(iorb,:)=Sigc_mo_tmp(iorb,iorb,:)
  enddo
  Sigc_mo = 0.5d0 * (Sigc_mo + transpose(Sigc_mo))
  call MOtoAO(nBas,nOrb,S,c,Sigc_mo,Sigc_ao)
  deallocate(Sigc_mo_tmp,Sigc_mo)

end subroutine 

