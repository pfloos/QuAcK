subroutine sigc_AO_basis_RHF(nBas,nOrb,nO,eta,shift,c,eqsGW_state,vMAT,nfreqs,ntimes, &
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
  double precision,intent(in)   :: c(nBas,nOrb)

! Local variables

  integer                       :: icluster
  integer                       :: iorb
  integer                       :: iE
  integer                       :: iE_
  integer                       :: nE_eval_global
  integer                       :: ncluster
  integer,allocatable           :: nE_per_cluster(:)

  double precision              :: step_E=2.5d-2
  double precision              :: chem_pot
  double precision,allocatable  :: E_eval_global(:)

  complex *16,allocatable       :: Sigc_mo(:,:,:)
  complex *16,allocatable       :: E_eval_global_cpx(:)

  type                          :: t_cluster
   integer                      :: nE
   integer                      :: nE_eval
   double precision,allocatable :: E_QP(:)
   double precision,allocatable :: E_eval(:)
  end type
  type(t_cluster),allocatable   :: clusters(:)

! Output variables

  double precision,intent(out)  :: Sigc_ao(nBas,nBas)

! Initialize variables
    
  chem_pot = 0.5d0*(eqsGW_state(nO)+eqsGW_state(nO+1))

  allocate(E_eval_global_cpx(1))
  call set_Eeval_cluster(nOrb,nOrb,1,shift,eqsGW_state,nE_eval_global,&
                         E_eval_global_cpx,chem_pot)
  deallocate(E_eval_global_cpx)
  allocate(E_eval_global_cpx(nE_eval_global))
  call set_Eeval_cluster(nOrb,nOrb,nE_eval_global,shift,eqsGW_state,&
                          nE_eval_global,E_eval_global_cpx,chem_pot)

  ! Run over unique energies
  allocate(Sigc_mo(nE_eval_global,nOrb,nOrb))
   call build_Sigmac_w_RHF(nOrb,nO,nE_eval_global,eta,0,E_eval_global_cpx,eqsGW_state,nfreqs,ntimes,&
                           wweight,wcoord,vMAT,Sigc_mo)


  ! TODO remove this printing
  write(*,*)
  write(*,*) 'w(adjusted with chem_pot)   Sigma_c(1,1)    Sigma_c(2,2) '
  do iE=1,nE_eval_global
   write(*,'(*(f10.5))') Real(E_eval_global_cpx(iE)),Real(Sigc_mo(iE,1,1)),Real(Sigc_mo(iE,2,2))
  enddo
  write(*,*)


! TODO interpolate Sigma

! TODO transform Sigma from MO to AO basis

! Delete dynamic arrays
 
  deallocate(Sigc_mo)
  deallocate(E_eval_global_cpx)


! TODO remove these lines
 Sigc_ao=0d0


end subroutine 

