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
  integer                       :: iE_
  integer                       :: ncluster
  integer                       :: nEnergies
  integer,allocatable           :: nE_per_cluster(:)

  double precision              :: shift=1d-2
  double precision              :: step_E=2.5d-2

  complex *16,allocatable       :: Sigc_mo_he(:,:,:)
  complex *16,allocatable       :: Sigc_mo_hh(:,:,:)

  type                          :: t_cluster
   integer                      :: nE
   integer                      :: nE_eval
   double precision,allocatable :: E_QP(:)
   complex *16,allocatable      :: E_eval(:)
  end type
  type(t_cluster),allocatable   :: clusters(:)

! Output variables

  double precision,intent(out)  :: Sigc_he(nBas,nBas)
  double precision,intent(out)  :: Sigc_hh(nBas,nBas)

!

  ! TODO
  !allocate(Sigc_mo_he(nOrb,nOrb),Sigc_mo_hh(nOrb,nOrb))

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

! For each cluster find the number of energies (e_EVAL) to be used in Sigma_c [i.e., Sigma_c(e_EVAl)]
! and allocate the energies to be evaluated array
  do icluster=1,ncluster
   clusters(icluster)%nE_eval=2
   do iE=2,clusters(icluster)%nE
     if( abs(clusters(icluster)%E_QP(iE)-clusters(icluster)%E_QP(iE-1)) >= 7.5d-2 ) then
      clusters(icluster)%nE_eval=clusters(icluster)%nE_eval+2
     else if( abs(clusters(icluster)%E_QP(iE)-clusters(icluster)%E_QP(iE-1)) < 7.5d-2 .and. &
              abs(clusters(icluster)%E_QP(iE)-clusters(icluster)%E_QP(iE-1)) >= 5d-2 ) then
      clusters(icluster)%nE_eval=clusters(icluster)%nE_eval+1
     else ! abs(clusters(icluster)%E_QP(iE)-clusters(icluster)%E_QP(iE-1)) < 0.05
      clusters(icluster)%nE_eval=clusters(icluster)%nE_eval
     endif   
   enddo
   clusters(icluster)%nE_eval=clusters(icluster)%nE_eval+2
   clusters(icluster)%nE_eval=clusters(icluster)%nE_eval+2*clusters(icluster)%nE
   allocate(clusters(icluster)%E_eval(clusters(icluster)%nE_eval))
  enddo

! For each cluster set the evaluation energies (e_EVAL) to be used in Sigma_c [i.e., Sigma_c(e_EVAl)]
  do icluster=1,ncluster
   clusters(icluster)%E_eval(:)=czero
   clusters(icluster)%E_eval(1)=clusters(icluster)%E_QP(1)-2d0*step_E
   clusters(icluster)%E_eval(2)=clusters(icluster)%E_QP(1)-step_E
   iE_=3
   do iE=2,clusters(icluster)%nE
     if( abs(clusters(icluster)%E_QP(iE)-clusters(icluster)%E_QP(iE-1)) >= 7.5d-2 ) then
       clusters(icluster)%E_eval(iE_)=clusters(icluster)%E_QP(iE-1)+step_E
       iE_=iE_+1
       clusters(icluster)%E_eval(iE_)=clusters(icluster)%E_QP(iE-1)+2d0*step_E
       iE_=iE_+1
     else if( abs(clusters(icluster)%E_QP(iE)-clusters(icluster)%E_QP(iE-1)) < 7.5d-2 .and. &
              abs(clusters(icluster)%E_QP(iE)-clusters(icluster)%E_QP(iE-1)) >= 5d-2 ) then
       clusters(icluster)%E_eval(iE_)=0.5d0*(clusters(icluster)%E_QP(iE)+clusters(icluster)%E_QP(iE-1))
       iE_=iE_+1
     else ! abs(clusters(icluster)%E_QP(iE)-clusters(icluster)%E_QP(iE-1)) < 0.05
       ! Nth
     endif   
   enddo
   clusters(icluster)%E_eval(iE_)=clusters(icluster)%E_QP(clusters(icluster)%nE)+step_E
   iE_=iE_+1
   clusters(icluster)%E_eval(iE_)=clusters(icluster)%E_QP(clusters(icluster)%nE)+2d0*step_E
   iE_=iE_+1
   ! Add the points where we want the self-energy [Sigma_c(E_QP)] with small +-shifts
   do iE=1,clusters(icluster)%nE
    clusters(icluster)%E_eval(iE_)=clusters(icluster)%E_QP(iE)-shift
    iE_=iE_+1
    clusters(icluster)%E_eval(iE_)=clusters(icluster)%E_QP(iE)+shift
    iE_=iE_+1
   enddo
  enddo

  ! TODO select unique points

  ! TODO run only once over unique points
  do icluster=1,ncluster
  allocate(Sigc_mo_he(clusters(icluster)%nE_eval,nOrb,nOrb),Sigc_mo_hh(clusters(icluster)%nE_eval,nOrb,nOrb))
  call build_Sigmac_w_HFB(nOrb,nOrb_twice,clusters(icluster)%nE_eval,0,clusters(icluster)%E_eval,eqsGWB_state, &
                          nfreqs,ntimes,wweight,wcoord,vMAT,U_QP,Sigc_mo_he,Sigc_mo_hh)
  write(*,*) 'New cluster'
  write(*,*) 'w   Sigma_c(1,1)    Sigma_c(2,2) '
   do iE=1,clusters(icluster)%nE_eval
    write(*,'(*(f10.5))') Real(clusters(icluster)%E_eval(iE)),Real(Sigc_mo_he(iE,1,1)),Real(Sigc_mo_he(iE,2,2))
   enddo
  deallocate(Sigc_mo_he,Sigc_mo_hh)
  enddo

!


! Delete dynamic arrays
 
  do icluster=1,ncluster
   deallocate(clusters(icluster)%E_QP)
   deallocate(clusters(icluster)%E_eval)
  enddo
  deallocate(clusters)




! TODO remove these lines
 Sigc_he=0d0
 Sigc_hh=0d0


end subroutine 

