subroutine sigc_AO_basis_HFB(nBas,nOrb,nOrb_twice,eta,shift,c,U_QP,eqsGWB_state,vMAT, &
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
  double precision,intent(in)   :: c(nBas,nOrb)

! Local variables

  integer                       :: iE
  integer                       :: nE_eval_global

  complex *16,allocatable       :: Sigc_mo_he(:,:,:)
  complex *16,allocatable       :: Sigc_mo_hh(:,:,:)
  complex *16,allocatable       :: E_eval_global_cpx(:)

! Output variables

  double precision,intent(out)  :: Sigc_ao_he(nBas,nBas)
  double precision,intent(out)  :: Sigc_ao_hh(nBas,nBas)

!
  allocate(E_eval_global_cpx(1))
  call set_Eeval_cluster(nOrb,nOrb_twice,1,shift,eqsGWB_state,nE_eval_global,&
                         E_eval_global_cpx,0d0)
  deallocate(E_eval_global_cpx)
  allocate(E_eval_global_cpx(nE_eval_global))
  call set_Eeval_cluster(nOrb,nOrb_twice,nE_eval_global,shift,eqsGWB_state,&
                          nE_eval_global,E_eval_global_cpx,0d0)

  ! Run over unique energies
  allocate(Sigc_mo_he(nE_eval_global,nOrb,nOrb),Sigc_mo_hh(nE_eval_global,nOrb,nOrb))
  call build_Sigmac_w_HFB(nOrb,nOrb_twice,nE_eval_global,eta,0,E_eval_global_cpx,eqsGWB_state, &
                          nfreqs,ntimes,wweight,wcoord,vMAT,U_QP,Sigc_mo_he,Sigc_mo_hh)


  ! TODO remove this printing
  write(*,*)
  write(*,*) 'w   Sigma_c(1,1)    Sigma_c(2,2) '
  do iE=1,nE_eval_global
   write(*,'(*(f10.5))') Real(E_eval_global_cpx(iE)),Real(Sigc_mo_he(iE,1,1)),Real(Sigc_mo_he(iE,2,2))
  enddo
  write(*,*)


! TODO interpolate Sigma

! TODO transform Sigma from MO to AO basis

! Delete dynamic arrays
 
  deallocate(Sigc_mo_he,Sigc_mo_hh)
  deallocate(E_eval_global_cpx)



! TODO remove these lines
 Sigc_ao_he=0d0
 Sigc_ao_hh=0d0


end subroutine 

