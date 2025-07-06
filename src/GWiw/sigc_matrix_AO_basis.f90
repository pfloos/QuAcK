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

  integer                       :: icluster
  integer                       :: iorb
  integer                       :: iE
  integer                       :: ncluster
  integer                       :: nEnergies
  integer,allocatable           :: nE_per_cluster(:)

  double precision              :: step_E=2.5d-2

  type                          :: t_cluster
   integer                      :: nE
   integer                      :: nE_eval
   double precision,allocatable :: E_QP(:)
  end type
  type(t_cluster),allocatable   :: clusters(:)

! Output variables

  double precision,intent(out)  :: Sigc_he(nBas,nBas)
  double precision,intent(out)  :: Sigc_hh(nBas,nBas)

!

! Set the number of clusters [i.e., energies separated more than (or equal to) 0.1 a.u.]
  ncluster=1
  do iorb=2,nOrb
   if(abs(eqsGWB_state(iorb)-eqsGWB_state(iorb-1))>=0.1d0) ncluster=ncluster+1
  enddo
  allocate(clusters(ncluster))

! For each cluster find the number of QP energies associated
  allocate(nE_per_cluster(ncluster))
  icluster=1; nE_per_cluster(:)=1;
  do iorb=2,nOrb
   if(abs(eqsGWB_state(iorb)-eqsGWB_state(iorb-1))>=0.1d0) then
    icluster=icluster+1
   else
    nE_per_cluster(icluster)=nE_per_cluster(icluster)+1
   endif
  enddo
  do icluster=1,ncluster
   clusters(icluster)%nE=nE_per_cluster(icluster)
  enddo
  deallocate(nE_per_cluster)

! For each cluster store the QP energies associated
  iorb=1
  do icluster=1,ncluster
   allocate(clusters(icluster)%E_QP(clusters(icluster)%nE))
   do iE=1,clusters(icluster)%nE
    clusters(icluster)%E_QP(iE)=eqsGWB_state(iorb)
    iorb=iorb+1
   enddo
  enddo

! For each cluster find the energies (e_EVAL) to be used in Sigma_c [i.e., Sigma_c(e_EVAl)]
  do icluster=1,ncluster
   do iE=1,clusters(icluster)%nE
!    Use nEnergies 
!    clusters(icluster)%E_QP(iE)
   enddo
  enddo

  do icluster=1,ncluster
   write(*,'(i5,*(f10.5))') clusters(icluster)%nE,clusters(icluster)%E_QP(:)
  enddo

! Delete dynamic arrays
 
  do icluster=1,ncluster
   deallocate(clusters(icluster)%E_QP)
  enddo
  deallocate(clusters)




! TODO remove these lines
 Sigc_he=0d0
 Sigc_hh=0d0



end subroutine 

! ---

