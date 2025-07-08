subroutine sigc_AO_basis_HFB(nBas,nOrb,nOrb_twice,eta,shift,c,U_QP,eqsGWB_state,S,vMAT, &
                             nfreqs,ntimes,wcoord,wweight,Sigc_ao_he,Sigc_ao_hh)

! Compute Sigma_c matrix in the AO basis

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nOrb_twice
  integer,intent(in)            :: nfreqs
  integer,intent(in)            :: ntimes

  double precision,intent(in)   :: eta
  double precision,intent(in)   :: shift
  double precision,intent(in)   :: wcoord(nfreqs),wweight(nfreqs)
  double precision,intent(in)   :: U_QP(nOrb_twice,nOrb_twice)
  double precision,intent(in)   :: eqsGWB_state(nOrb_twice)
  double precision,intent(in)   :: vMAT(nOrb*nOrb,nOrb*nOrb)
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: c(nBas,nOrb)

! Local variables

  integer                       :: iE
  integer                       :: iorb,jorb
  integer                       :: nE_eval_global

  double precision,allocatable  :: Sigc_mo_tmp(:,:,:)
  double precision,allocatable  :: Sigc_mo(:,:)

  complex *16,allocatable       :: Sigc_mo_he_cpx(:,:,:)
  complex *16,allocatable       :: Sigc_mo_hh_cpx(:,:,:)
  complex *16,allocatable       :: E_eval_global_cpx(:)

! Output variables

  double precision,intent(out)  :: Sigc_ao_he(nBas,nBas)
  double precision,intent(out)  :: Sigc_ao_hh(nBas,nBas)

! Se energies using cluster method or just using two shifts
  if(.false.) then
   allocate(E_eval_global_cpx(1))
   call set_Eeval_cluster(nOrb,nOrb_twice,1,shift,eqsGWB_state,nE_eval_global,&
                          E_eval_global_cpx,0d0)
   deallocate(E_eval_global_cpx)
   allocate(E_eval_global_cpx(nE_eval_global))
   call set_Eeval_cluster(nOrb,nOrb_twice,nE_eval_global,shift,eqsGWB_state,&
                           nE_eval_global,E_eval_global_cpx,0d0)
  else
   nE_eval_global=2*nOrb
   allocate(E_eval_global_cpx(nE_eval_global))
   do iorb=1,nOrb
    E_eval_global_cpx(2*iorb-1)=eqsGWB_state(iorb)-shift
    E_eval_global_cpx(2*iorb)  =eqsGWB_state(iorb)+shift
   enddo
  endif

  ! Run over unique energies
  allocate(Sigc_mo_he_cpx(nE_eval_global,nOrb,nOrb),Sigc_mo_hh_cpx(nE_eval_global,nOrb,nOrb))
  call build_Sigmac_w_HFB(nOrb,nOrb_twice,nE_eval_global,eta,0,E_eval_global_cpx,eqsGWB_state, &
                          nfreqs,ntimes,wweight,wcoord,vMAT,U_QP,Sigc_mo_he_cpx,Sigc_mo_hh_cpx)
  deallocate(E_eval_global_cpx)

! Interpolate and transform Sigma from MO to AO basis
! TODO interpolate Sigma with clusters method
  allocate(Sigc_mo_tmp(nOrb,nOrb,nOrb))
  allocate(Sigc_mo(nOrb,nOrb))
  ! Sigma_c_he
  do iorb=1,nOrb
   Sigc_mo_tmp(iorb,:,:)=0.5d0*(Real(Sigc_mo_he_cpx(2*iorb-1,:,:))+Real(Sigc_mo_he_cpx(2*iorb,:,:)))
   Sigc_mo_tmp(iorb,:,:)=0.5d0*(Sigc_mo_tmp(iorb,:,:)+Transpose(Sigc_mo_tmp(iorb,:,:)))
  enddo
  deallocate(Sigc_mo_he_cpx)
  do iorb=1,nOrb
   do jorb=iorb,nOrb
    Sigc_mo(iorb,jorb)=0.5d0*(Sigc_mo_tmp(iorb,iorb,jorb)+Sigc_mo_tmp(jorb,iorb,jorb))
    Sigc_mo(jorb,iorb)=Sigc_mo(iorb,jorb)
   enddo
  enddo
  call MOtoAO(nBas,nOrb,S,c,Sigc_mo,Sigc_ao_he)
  ! Sigma_c_hh
  do iorb=1,nOrb
   Sigc_mo_tmp(iorb,:,:)=0.5d0*(Real(Sigc_mo_hh_cpx(2*iorb-1,:,:))+Real(Sigc_mo_hh_cpx(2*iorb,:,:)))
   Sigc_mo_tmp(iorb,:,:)=0.5d0*(Sigc_mo_tmp(iorb,:,:)+Transpose(Sigc_mo_tmp(iorb,:,:)))
  enddo
  deallocate(Sigc_mo_hh_cpx)
  do iorb=1,nOrb
   do jorb=iorb,nOrb
    Sigc_mo(iorb,jorb)=0.5d0*(Sigc_mo_tmp(iorb,iorb,jorb)+Sigc_mo_tmp(jorb,iorb,jorb))
    Sigc_mo(jorb,iorb)=Sigc_mo(iorb,jorb)
   enddo
  enddo
  call MOtoAO(nBas,nOrb,S,c,Sigc_mo,Sigc_ao_hh)
  deallocate(Sigc_mo_tmp,Sigc_mo)


end subroutine 

