subroutine scGF2B_AO_itau_iw(nBas,nOrb,nOrb_twice,maxSCF,thresh_in,maxDIIS,dolinGF2,restart_scGF2B,verbose_scGF2B,chem_pot_scG, &
                            no_h_hfb,ENuc,Hc,S,X_in,P_in,Pan_in,cHFB,eQP_state,chem_pot,sigma,nfreqs,wcoord,wweight,   &
                            U_QP,vMAT,ERI_AO)

! Restricted scGF2B

  implicit none
  include 'parameters.h'

! Input variables
 
  logical,intent(in)            :: dolinGF2
  logical,intent(in)            :: no_h_hfb
  logical,intent(in)            :: restart_scGF2B
  logical,intent(in)            :: verbose_scGF2B
  logical,intent(in)            :: chem_pot_scG

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nOrb_twice
  integer,intent(in)            :: maxSCF
  integer,intent(in)            :: maxDIIS

  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: thresh_in
  double precision,intent(in)   :: sigma
  double precision,intent(in)   :: Hc(nBas,nBas)
  double precision,intent(in)   :: P_in(nBas,nBas)
  double precision,intent(in)   :: Pan_in(nBas,nBas)
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: X_in(nBas,nOrb)
  double precision,intent(in)   :: vMAT(nBas*nBas,nBas*nBas)
  double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: U_QP(nOrb_twice,nOrb_twice)

! Local variables

  logical                       :: file_exists
  logical                       :: read_HFB_chkp

  integer                       :: ifreq,itau
  integer                       :: ibas,jbas,kbas,lbas,obas,qbas
  integer                       :: iorb,jorb,korb,lorb
  integer                       :: iter,iter_hfb
  integer                       :: verbose
  integer                       :: ntimes
  integer                       :: nBasSq
  integer                       :: nBas_twice,nBas_twiceSq
  integer                       :: ntimes_twice
  integer                       :: nfreqs_int
  integer                       :: imax_error_gw2gt
  integer                       :: imax_error_st2sw
  integer                       :: idiis_indexR,n_diisR
  integer                       :: idiis_index,n_diis
  integer                       :: nBas2Sqntimes2

  double precision              :: rcond
  double precision              :: rcondR
  double precision              :: alpha_mixing
  double precision              :: N_anom
  double precision              :: eta,diff_Rao
  double precision              :: thrs_N,thrs_Ngrad,thrs_Rao
  double precision              :: nElectrons
  double precision              :: chem_pot_saved
  double precision              :: error_gw2gt
  double precision              :: error_st2sw
  double precision              :: max_error_gw2gt
  double precision              :: sum_error_gw2gt
  double precision              :: max_error_st2sw
  double precision              :: EcGM,Ehfbl,Ecore,Eh,Ex,Epair
  double precision              :: trace_1_rdm
  double precision,external     :: trace_matrix
  double precision              :: start_scGF2Bitauiw,end_scGF2Bitauiw,t_scGF2Bitauiw
  double precision,allocatable  :: Occ(:)
  double precision,allocatable  :: cHFBinv(:,:)
  double precision,allocatable  :: cNO(:,:)
  double precision,allocatable  :: U_mo(:,:)
  double precision,allocatable  :: Mat1(:,:)
  double precision,allocatable  :: Mat2(:,:)
  double precision,allocatable  :: R_ao(:,:)
  double precision,allocatable  :: R_ao_iter(:,:)
  double precision,allocatable  :: R_ao_hfb(:,:)
  double precision,allocatable  :: R_ao_old(:,:)
  double precision,allocatable  :: H_ao_hfb(:,:)
  double precision,allocatable  :: R_mo(:,:)
  double precision,allocatable  :: cHFB_gorkov(:,:)
  double precision,allocatable  :: U_QP_tmp(:,:)
  double precision,allocatable  :: err_currentR(:)
  double precision,allocatable  :: err_diisR(:,:)
  double precision,allocatable  :: R_ao_extrap(:)
  double precision,allocatable  :: R_ao_old_diis(:,:)
  double precision,allocatable  :: wcoord_int(:),wweight_int(:)
  double precision,allocatable  :: tweight(:),tcoord(:)
  double precision,allocatable  :: sint2w_weight(:,:)
  double precision,allocatable  :: cost2w_weight(:,:)
  double precision,allocatable  :: cosw2t_weight(:,:)
  double precision,allocatable  :: sinw2t_weight(:,:)
  double precision,allocatable  :: vMAT_mo(:,:)
  double precision,allocatable  :: ERI_MO(:,:,:,:)

  complex*16                    :: product
  complex*16                    :: weval_cpx
  complex*16                    :: EcGM_itau
  complex*16,allocatable        :: wtest(:)
  complex*16,allocatable        :: Sigma_c_w_ao(:,:,:)
  complex*16,allocatable        :: Sigma_c_c(:,:),Sigma_c_s(:,:)
  complex*16,allocatable        :: Sigma_c_plus(:,:),Sigma_c_minus(:,:)
  complex*16,allocatable        :: G_ao_tmp(:,:)
  complex*16,allocatable        :: Mat_gorkov_tmp(:,:)
  complex*16,allocatable        :: Mat_gorkov_tmp2(:,:)
  complex*16,allocatable        :: G_itau_extrap(:)
  complex*16,allocatable        :: err_current(:)
  complex*16,allocatable        :: G_ao1(:,:)
  complex*16,allocatable        :: G_ao2(:,:)
  complex*16,allocatable        :: G_ao3(:,:)
  complex*16,allocatable        :: G_ao_itau(:,:,:)
  complex*16,allocatable        :: G_ao_itau_old(:,:,:)
  complex*16,allocatable        :: G_ao_itau_hfb(:,:,:)
  complex*16,allocatable        :: err_diis(:,:)
  complex*16,allocatable        :: G_itau_old_diis(:,:)
  complex*16,allocatable        :: G_ao_iw_hfb(:,:,:)
  complex*16,allocatable        :: DeltaG_ao_iw(:,:,:)
  complex*16,allocatable        :: Sigma_c_he(:,:,:)
  complex*16,allocatable        :: Sigma_c_hh(:,:,:)
  complex*16,allocatable        :: Sigma_c_ee(:,:,:)
  complex*16,allocatable        :: Sigma_c_eh(:,:,:)
  complex*16,allocatable        :: Ainter(:,:,:,:)
  complex*16,allocatable        :: Binter(:,:,:,:)
  complex*16,allocatable        :: Cinter(:,:,:,:)

! Output variables
  integer,intent(inout)         :: nfreqs
  double precision,intent(inout):: chem_pot
  double precision,intent(inout):: eQP_state(nOrb_twice)
  double precision,intent(inout):: wcoord(nfreqs)
  double precision,intent(inout):: wweight(nfreqs)
  double precision,intent(inout):: cHFB(nBas,nOrb)
  
 call wall_time(start_scGF2Bitauiw)

 write(*,*)     
 write(*,*)'*****************************************'
 write(*,*)'*   scGF2B ( using it and iw grids )    *'
 write(*,*)'*****************************************'
 write(*,*)

 write(*,*)
 write(*,'(A33,1X,F16.10,A3)') ' Initial chemical potential  = ',chem_pot,' au'
 if(chem_pot_scG) then
  write(*,'(A)') '   Adjusting the chemical potential is activated'
 else
  write(*,'(A)') '   Adjusting the chemical potential is deactivated'
 endif
 write(*,*)

 write(*,'(a,F15.8)') '  emin ',abs(eQP_state(nOrb))
 write(*,'(a,F15.8)') '  emax ',abs(eQP_state(nOrb)+eQP_state(1))

 ! Initialize variables
 ntimes=nfreqs
 nBasSq=nBas*nBas
 ntimes_twice=2*ntimes
 nBas_twice=2*nBas
 nBas_twiceSq=nBas_twice*nBas_twice
 verbose=0
 if(verbose_scGF2B) verbose=1     
 eta=0d0
 thrs_N=1d-10
 thrs_Ngrad=1d-6
 thrs_Rao=thresh_in
 nElectrons=0d0
 do ibas=1,nBas
  do jbas=1,nBas
   nElectrons=nElectrons+P_in(ibas,jbas)*S(ibas,jbas)
  enddo
 enddo
 nElectrons=nint(0.5d0*nElectrons) ! Here we prefer to use 1 spin-channel
 chem_pot_saved=chem_pot
 alpha_mixing=0.6d0
 rcond=0d0
 n_diis=0
 nBas2Sqntimes2=nBas_twiceSq*ntimes_twice

 ! Allocate arrays
 allocate(Occ(nOrb))
 allocate(cNO(nBas,nOrb))
 allocate(cHFBinv(nOrb,nBas))
 allocate(cHFB_gorkov(nBas_twice,nOrb_twice))
 allocate(U_QP_tmp(nOrb_twice,nOrb_twice))
 allocate(U_mo(nOrb,nOrb))
 allocate(R_ao(nBas_twice,nBas_twice))
 allocate(R_ao_iter(nBas_twice,nBas_twice))
 allocate(R_ao_hfb(nBas_twice,nBas_twice))
 allocate(R_ao_old(nBas_twice,nBas_twice))
 allocate(R_mo(nOrb_twice,nOrb_twice))
 allocate(H_ao_hfb(nBas_twice,nBas_twice))
 allocate(G_ao_tmp(nBas,nBas))
 allocate(Mat_gorkov_tmp(nBas_twice,nBas_twice))
 allocate(Mat_gorkov_tmp2(nBas_twice,nBas_twice))
 allocate(G_ao_iw_hfb(nfreqs,nBas_twice,nBas_twice))
 allocate(DeltaG_ao_iw(nfreqs,nBas_twice,nBas_twice))
 allocate(G_ao_itau_hfb(ntimes_twice,nBas_twice,nBas_twice))
 allocate(G_ao_itau(ntimes_twice,nBas_twice,nBas_twice))
 allocate(G_ao_itau_old(ntimes_twice,nBas_twice,nBas_twice))
 allocate(G_ao1(nBas,nBas))
 allocate(G_ao2(nBas,nBas))
 allocate(G_ao3(nBas,nBas))
 allocate(Sigma_c_w_ao(nfreqs,nBas_twice,nBas_twice))
 allocate(Sigma_c_plus(nBas_twice,nBas_twice)) 
 allocate(Sigma_c_minus(nBas_twice,nBas_twice)) 
 allocate(Sigma_c_c(nBas_twice,nBas_twice))      
 allocate(Sigma_c_s(nBas_twice,nBas_twice))
 allocate(tweight(ntimes),tcoord(ntimes))
 allocate(sint2w_weight(nfreqs,ntimes))
 allocate(cost2w_weight(nfreqs,ntimes))
 allocate(cosw2t_weight(ntimes,nfreqs))
 allocate(sinw2t_weight(ntimes,nfreqs))
 allocate(Mat1(nOrb,nOrb),Mat2(nOrb,nOrb))
 allocate(Ainter(nBas,nBas,nBas,nBas))
 allocate(Binter(nBas,nBas,nBas,nBas))
 allocate(Cinter(nBas,nBas,nBas,nBas))
 allocate(err_currentR(1))
 allocate(err_diisR(1,1))
 allocate(R_ao_extrap(1))
 allocate(R_ao_old_diis(1,1))
 allocate(G_itau_extrap(1))
 allocate(err_current(1))
 allocate(err_diis(1,1))
 allocate(G_itau_old_diis(1,1))
 if(maxDIIS>0) then
  deallocate(err_currentR)
  deallocate(R_ao_extrap)
  deallocate(err_diisR)
  deallocate(R_ao_old_diis)
  deallocate(G_itau_extrap)
  deallocate(err_current)
  deallocate(err_diis)
  deallocate(G_itau_old_diis)
  allocate(err_currentR(nBas_twiceSq))
  allocate(R_ao_extrap(nBas_twiceSq))
  allocate(err_diisR(nBas_twiceSq,maxDIIS))
  allocate(R_ao_old_diis(nBas_twiceSq,maxDIIS))
  allocate(G_itau_extrap(nBas2Sqntimes2))
  allocate(err_current(nBas2Sqntimes2))
  allocate(err_diis(nBas2Sqntimes2,maxDIIS))
  allocate(G_itau_old_diis(nBas2Sqntimes2,maxDIIS))
 endif
 err_diis=czero
 G_itau_old_diis=czero

 ! Initialize arrays
 DeltaG_ao_iw(:,:,:)=czero
 Mat1(1:nOrb,1:nOrb)=U_QP(1:nOrb,1:nOrb)
 Mat2(1:nOrb,1:nOrb)=U_QP(nOrb+1:nOrb_twice,1:nOrb)
 R_ao_hfb(1:nBas           ,1:nBas           ) = 0.5d0*P_in(1:nBas,1:nBas)
 R_ao_hfb(nBas+1:nBas_twice,nBas+1:nBas_twice) = matmul(X_in(1:nBas,1:nOrb), transpose(X_in(1:nBas,1:nOrb)))-0.5d0*P_in(1:nBas,1:nBas)
 R_ao_hfb(1:nBas           ,nBas+1:nBas_twice) = Pan_in(1:nBas,1:nBas)
 R_ao_hfb(nBas+1:nBas_twice,1:nBas           ) = Pan_in(1:nBas,1:nBas)
 R_ao=R_ao_hfb
 R_ao_iter=R_ao_hfb
 cHFBinv=matmul(transpose(cHFB),S)
 cHFB_gorkov=0d0
 cHFB_gorkov(1:nBas           ,1:nOrb           ) = cHFB(1:nBas,1:nOrb)
 cHFB_gorkov(nBas+1:nBas_twice,nOrb+1:nOrb_twice) = cHFB(1:nBas,1:nOrb)
 H_ao_hfb=0d0
 H_ao_hfb(1:nBas,1:nBas)=Hc(1:nBas,1:nBas)
 Ehfbl=0d0
 do ibas=1,nBas
  do jbas=1,nBas
   obas=nBas+1+(jbas-1)
   Ehfbl=Ehfbl+2d0*R_ao(ibas,jbas)*Hc(ibas,jbas)
   do kbas=1,nBas
    do lbas=1,nBas
     qbas=nBas+1+(lbas-1)
     H_ao_hfb(ibas,jbas)=H_ao_hfb(ibas,jbas)+2d0*R_ao(kbas,lbas)*vMAT(1+(lbas-1)+(kbas-1)*nBas,1+(jbas-1)+(ibas-1)*nBas) &
                        -R_ao(kbas,lbas)*vMAT(1+(jbas-1)+(kbas-1)*nBas,1+(lbas-1)+(ibas-1)*nBas)
     H_ao_hfb(ibas,obas)=H_ao_hfb(ibas,obas)+sigma*R_ao(kbas,qbas)*vMAT(1+(ibas-1)+(kbas-1)*nBas,1+(jbas-1)+(lbas-1)*nBas)
     Ehfbl=Ehfbl+2d0*R_ao(kbas,lbas)*R_ao(ibas,jbas)*vMAT(1+(lbas-1)+(kbas-1)*nBas,1+(jbas-1)+(ibas-1)*nBas) &
          -R_ao(kbas,lbas)*R_ao(ibas,jbas)*vMAT(1+(jbas-1)+(kbas-1)*nBas,1+(lbas-1)+(ibas-1)*nBas)           &
          +sigma*R_ao(kbas,qbas)*R_ao(ibas,obas)*vMAT(1+(ibas-1)+(kbas-1)*nBas,1+(jbas-1)+(lbas-1)*nBas)
    enddo
   enddo
  enddo
 enddo
 H_ao_hfb(nBas+1:nBas_twice,nBas+1:nBas_twice) = -H_ao_hfb(1:nBas,1:nBas           )
 H_ao_hfb(nBas+1:nBas_twice,1:nBas           ) =  H_ao_hfb(1:nBas,nBas+1:nBas_twice)

 ! Read grids 
 call read_scGX_grids(ntimes,nfreqs,tcoord,tweight,wcoord,wweight,sint2w_weight,cost2w_weight, &
                      cosw2t_weight,sinw2t_weight,verbose)

 ! Build Go(i w)
 do ifreq=1,nfreqs
  weval_cpx=im*wcoord(ifreq)
  ! G_he(iw)
  call G_AO_RHFB_w(nBas,nOrb,nOrb_twice,eta,cHFB,eQP_state,weval_cpx, Mat1, Mat1, Mat2, Mat2,G_ao_tmp)
  G_ao_iw_hfb(ifreq,1:nBas           ,1:nBas           ) = G_ao_tmp(1:nBas,1:nBas)
  ! G_hh(iw)
  call G_AO_RHFB_w(nBas,nOrb,nOrb_twice,eta,cHFB,eQP_state,weval_cpx, Mat1, Mat2,-Mat2, Mat1,G_ao_tmp)
  G_ao_iw_hfb(ifreq,1:nBas           ,nBas+1:nBas_twice) = G_ao_tmp(1:nBas,1:nBas)
  ! G_ee(iw)
  call G_AO_RHFB_w(nBas,nOrb,nOrb_twice,eta,cHFB,eQP_state,weval_cpx, Mat2, Mat1, Mat1,-Mat2,G_ao_tmp)
  G_ao_iw_hfb(ifreq,nBas+1:nBas_twice,1:nBas           ) = G_ao_tmp(1:nBas,1:nBas)
  ! G_eh(iw)
  call G_AO_RHFB_w(nBas,nOrb,nOrb_twice,eta,cHFB,eQP_state,weval_cpx, Mat2, Mat2, Mat1, Mat1,G_ao_tmp)
  G_ao_iw_hfb(ifreq,nBas+1:nBas_twice,nBas+1:nBas_twice) = G_ao_tmp(1:nBas,1:nBas)
 enddo
 ! Build Go(i tau)
 do itau=1,ntimes
  ! tau > 0
  ! G_he(i tau)
  call G_AO_RHFB_itau(nBas,nOrb,nOrb_twice, tcoord(itau),G_ao_tmp,cHFB,eQP_state, Mat1, Mat1, Mat2, Mat2)
  G_ao_itau_hfb(2*itau-1,1:nBas           ,1:nBas           ) = G_ao_tmp(1:nBas,1:nBas)
  ! G_hh(i tau)
  call G_AO_RHFB_itau(nBas,nOrb,nOrb_twice, tcoord(itau),G_ao_tmp,cHFB,eQP_state, Mat1, Mat2,-Mat2, Mat1)
  G_ao_itau_hfb(2*itau-1,1:nBas           ,nBas+1:nBas_twice) = G_ao_tmp(1:nBas,1:nBas)
  ! G_ee(i tau)
  call G_AO_RHFB_itau(nBas,nOrb,nOrb_twice, tcoord(itau),G_ao_tmp,cHFB,eQP_state, Mat2, Mat1, Mat1,-Mat2)
  G_ao_itau_hfb(2*itau-1,nBas+1:nBas_twice,1:nBas           ) = G_ao_tmp(1:nBas,1:nBas)
  ! G_eh(i tau)
  call G_AO_RHFB_itau(nBas,nOrb,nOrb_twice, tcoord(itau),G_ao_tmp,cHFB,eQP_state, Mat2, Mat2, Mat1, Mat1)
  G_ao_itau_hfb(2*itau-1,nBas+1:nBas_twice,nBas+1:nBas_twice) = G_ao_tmp(1:nBas,1:nBas)
  ! tau < 0
  ! G_he(i tau)
  call G_AO_RHFB_itau(nBas,nOrb,nOrb_twice,-tcoord(itau),G_ao_tmp,cHFB,eQP_state, Mat1, Mat1, Mat2, Mat2)
  G_ao_itau_hfb(2*itau  ,1:nBas           ,1:nBas           ) = G_ao_tmp(1:nBas,1:nBas)
  ! G_hh(i tau)
  call G_AO_RHFB_itau(nBas,nOrb,nOrb_twice,-tcoord(itau),G_ao_tmp,cHFB,eQP_state, Mat1, Mat2,-Mat2, Mat1)
  G_ao_itau_hfb(2*itau  ,1:nBas           ,nBas+1:nBas_twice) = G_ao_tmp(1:nBas,1:nBas)
  ! G_ee(i tau)
  call G_AO_RHFB_itau(nBas,nOrb,nOrb_twice,-tcoord(itau),G_ao_tmp,cHFB,eQP_state, Mat2, Mat1, Mat1,-Mat2)
  G_ao_itau_hfb(2*itau  ,nBas+1:nBas_twice,1:nBas           ) = G_ao_tmp(1:nBas,1:nBas)
  ! G_eh(i tau)
  call G_AO_RHFB_itau(nBas,nOrb,nOrb_twice,-tcoord(itau),G_ao_tmp,cHFB,eQP_state, Mat2, Mat2, Mat1, Mat1)
  G_ao_itau_hfb(2*itau  ,nBas+1:nBas_twice,nBas+1:nBas_twice) = G_ao_tmp(1:nBas,1:nBas)
 enddo
 ! Initialize G(i tau)
 G_ao_itau(:,:,:)=G_ao_itau_hfb(:,:,:)
 G_ao_itau_old(:,:,:)=G_ao_itau_hfb(:,:,:)
 ! Check error in the Fourier transformation Go(iw) -> Go(it)
 write(*,*)
 write(*,'(a)') ' Error test for the Go(iw) -> Go(it) transformation'
 write(*,*)
 max_error_gw2gt=-1d0
 sum_error_gw2gt=0d0
 imax_error_gw2gt=1
 do itau=1,ntimes
  ! tau > 0
  Mat_gorkov_tmp=czero
  do ifreq=1,nfreqs
   Mat_gorkov_tmp(:,:) = Mat_gorkov_tmp(:,:) + im*cosw2t_weight(itau,ifreq)*Real(G_ao_iw_hfb(ifreq,:,:))  &
                                             - im*sinw2t_weight(itau,ifreq)*Aimag(G_ao_iw_hfb(ifreq,:,:))
  enddo
  Mat_gorkov_tmp(:,:)=abs(Mat_gorkov_tmp(:,:)-G_ao_itau_hfb(2*itau-1,:,:))
  error_gw2gt=real(sum(Mat_gorkov_tmp(:,:)))
  ! tau < 0
  Mat_gorkov_tmp=czero
  do ifreq=1,nfreqs
   Mat_gorkov_tmp(:,:) = Mat_gorkov_tmp(:,:) + im*cosw2t_weight(itau,ifreq)*Real(G_ao_iw_hfb(ifreq,:,:))  &
                                             + im*sinw2t_weight(itau,ifreq)*Aimag(G_ao_iw_hfb(ifreq,:,:))
  enddo
  Mat_gorkov_tmp(:,:)=abs(Mat_gorkov_tmp(:,:)-G_ao_itau_hfb(2*itau  ,:,:))
  error_gw2gt=error_gw2gt+real(sum(Mat_gorkov_tmp(:,:)))
  sum_error_gw2gt=sum_error_gw2gt+error_gw2gt
  if(error_gw2gt>max_error_gw2gt) then
   imax_error_gw2gt=itau
   max_error_gw2gt=error_gw2gt
  endif
 enddo
 write(*,'(a,*(f20.8))') ' Sum error ',sum_error_gw2gt
 write(*,'(a,f20.8,a,2f20.8,a)') ' Max CAE   ',max_error_gw2gt,' is in the time +/-',0d0,tcoord(imax_error_gw2gt),'i'
 write(*,'(a,*(f20.8))') ' MAE       ',sum_error_gw2gt/(2*ntimes*nBas_twice*nBas_twice)

 ! If required, read restart files
 if(restart_scGF2B) then
  inquire(file='read_HFB_scGWB', exist=file_exists)
  read_HFB_chkp=.false.
  if(file_exists) read_HFB_chkp=.true.
  call read_scGWB_restart(nBas_twice,nfreqs,ntimes_twice,chem_pot,R_ao,R_ao_hfb,G_ao_iw_hfb,G_ao_itau,G_ao_itau_hfb,read_HFB_chkp)
  R_ao_iter=R_ao
  G_ao_itau_old(:,:,:)=G_ao_itau(:,:,:)
 endif

!-------------!
! scGF2B loop !
!-------------!

 iter=0
 iter_hfb=0
 do
  iter=iter+1
 
  ! Build Sigma_c(i w) [Eqs. 12-18 in PRB, 109, 255101 (2024)]
  ! and compute Galitskii-Migdal Ec energy using the time grid
  EcGM_itau=czero
  Sigma_c_w_ao=czero
  do itau=1,ntimes
   Sigma_c_plus=czero
   Sigma_c_minus=czero

   ! G(i tau) G(i tau) G(-i tau) -> Sigma_c(i tau)
    ! Sigma_he_2prime
    G_ao1(:,:) = G_ao_itau(2*itau-1,           1:nBas,           1:nBas)  ! he   itau
    G_ao2(:,:) = G_ao_itau(2*itau-1,           1:nBas,           1:nBas)  ! he   itau
    G_ao3(:,:) = G_ao_itau(2*itau  ,           1:nBas,           1:nBas)  ! he  -itau
    Ainter=czero;Binter=czero;Cinter=czero;
    call Sigma_c_GF2B_he_prime(nBas,Ainter,Binter,Cinter,G_ao1,G_ao2,G_ao3,ERI_AO,Sigma_c_plus(1:nBas,1:nBas)) 
    ! Sigma_he_2primeprime
    G_ao1(:,:) = G_ao_itau(2*itau-1,           1:nBas,           1:nBas)  ! he   itau
    G_ao2(:,:) = G_ao_itau(2*itau-1,           1:nBas,nBas+1:nBas_twice)  ! hh   itau
    G_ao3(:,:) = G_ao_itau(2*itau  ,nBas+1:nBas_twice,           1:nBas)  ! ee  -itau
    Ainter=czero;Binter=czero;Cinter=czero;
    call Sigma_c_GF2B_he_primeprime(nBas,Ainter,Binter,Cinter,G_ao1,G_ao2,G_ao3,ERI_AO,Sigma_c_plus(1:nBas,1:nBas)) 
    ! Sigma_eh_2prime
    G_ao1(:,:) = G_ao_itau(2*itau-1,nBas+1:nBas_twice,nBas+1:nBas_twice)  ! eh   itau
    G_ao2(:,:) = G_ao_itau(2*itau-1,nBas+1:nBas_twice,nBas+1:nBas_twice)  ! eh   itau
    G_ao3(:,:) = G_ao_itau(2*itau  ,nBas+1:nBas_twice,nBas+1:nBas_twice)  ! eh  -itau
    Ainter=czero;Binter=czero;Cinter=czero;
    call Sigma_c_GF2B_he_prime(nBas,Ainter,Binter,Cinter,G_ao1,G_ao2,G_ao3,ERI_AO,Sigma_c_plus(nBas+1:nBas_twice,nBas+1:nBas_twice)) ! Valid for real orbitals 
    ! Sigma_eh_2primeprime
    G_ao1(:,:) = G_ao_itau(2*itau-1,nBas+1:nBas_twice,nBas+1:nBas_twice)  ! eh   itau
    G_ao2(:,:) = G_ao_itau(2*itau-1,           1:nBas,nBas+1:nBas_twice)  ! hh   itau
    G_ao3(:,:) = G_ao_itau(2*itau  ,nBas+1:nBas_twice,           1:nBas)  ! ee  -itau
    Ainter=czero;Binter=czero;Cinter=czero;
    call Sigma_c_GF2B_eh_primeprime(nBas,Ainter,Binter,Cinter,G_ao1,G_ao2,G_ao3,ERI_AO,Sigma_c_plus(nBas+1:nBas_twice,nBas+1:nBas_twice))
    ! Sigma_hh_2prime
    G_ao1(:,:) = G_ao_itau(2*itau-1,           1:nBas,nBas+1:nBas_twice)  ! hh   itau
    G_ao2(:,:) = G_ao_itau(2*itau-1,           1:nBas,           1:nBas)  ! he   itau
    G_ao3(:,:) = G_ao_itau(2*itau  ,           1:nBas,           1:nBas)  ! he  -itau
    Ainter=czero;Binter=czero;Cinter=czero;
    call Sigma_c_GF2B_hh_prime(nBas,Ainter,Binter,Cinter,G_ao1,G_ao2,G_ao3,ERI_AO,Sigma_c_plus(1:nBas,nBas+1:nBas_twice))
    ! Sigma_hh_2primeprime
    G_ao1(:,:) = G_ao_itau(2*itau-1,           1:nBas,nBas+1:nBas_twice)  ! hh   itau
    G_ao2(:,:) = G_ao_itau(2*itau-1,           1:nBas,nBas+1:nBas_twice)  ! hh   itau
    G_ao3(:,:) = G_ao_itau(2*itau  ,nBas+1:nBas_twice,           1:nBas)  ! ee  -itau
    Ainter=czero;Binter=czero;Cinter=czero;
    call Sigma_c_GF2B_hh_primeprime(nBas,Ainter,Binter,Cinter,G_ao1,G_ao2,G_ao3,ERI_AO,Sigma_c_plus(1:nBas,nBas+1:nBas_twice))
    ! Sigma_ee_2prime
    G_ao1(:,:) = G_ao_itau(2*itau-1,nBas+1:nBas_twice,           1:nbas)  ! ee   itau
    G_ao2(:,:) = G_ao_itau(2*itau-1,           1:nBas,           1:nBas)  ! he   itau
    G_ao3(:,:) = G_ao_itau(2*itau  ,           1:nBas,           1:nBas)  ! he  -itau
    Ainter=czero;Binter=czero;Cinter=czero;
    call Sigma_c_GF2B_ee_prime(nBas,Ainter,Binter,Cinter,G_ao1,G_ao2,G_ao3,ERI_AO,Sigma_c_plus(nBas+1:nBas_twice,1:nBas))
    ! Sigma_ee_2primeprime
    G_ao1(:,:) = G_ao_itau(2*itau-1,nBas+1:nBas_twice,           1:nBas)  ! ee   itau
    G_ao2(:,:) = G_ao_itau(2*itau-1,           1:nBas,nBas+1:nBas_twice)  ! hh   itau
    G_ao3(:,:) = G_ao_itau(2*itau  ,nBas+1:nBas_twice,           1:nBas)  ! ee  -itau
    Ainter=czero;Binter=czero;Cinter=czero;
    call Sigma_c_GF2B_ee_primeprime(nBas,Ainter,Binter,Cinter,G_ao1,G_ao2,G_ao3,ERI_AO,Sigma_c_plus(nBas+1:nBas_twice,1:nBas))

   ! G(-i tau) G(-i tau) G(i tau) -> Sigma_c(-i tau)  
    ! Sigma_he_2prime
    G_ao1(:,:) = G_ao_itau(2*itau  ,           1:nBas,           1:nBas)  ! he  -itau
    G_ao2(:,:) = G_ao_itau(2*itau  ,           1:nBas,           1:nBas)  ! he  -itau
    G_ao3(:,:) = G_ao_itau(2*itau-1,           1:nBas,           1:nBas)  ! he   itau
    Ainter=czero;Binter=czero;Cinter=czero;
    call Sigma_c_GF2B_he_prime(nBas,Ainter,Binter,Cinter,G_ao1,G_ao2,G_ao3,ERI_AO,Sigma_c_minus(1:nBas,1:nBas)) 
    ! Sigma_he_2primeprime
    G_ao1(:,:) = G_ao_itau(2*itau  ,           1:nBas,           1:nBas)  ! he  -itau
    G_ao2(:,:) = G_ao_itau(2*itau  ,           1:nBas,nBas+1:nBas_twice)  ! hh  -itau
    G_ao3(:,:) = G_ao_itau(2*itau-1,nBas+1:nBas_twice,           1:nBas)  ! ee   itau
    Ainter=czero;Binter=czero;Cinter=czero;
    call Sigma_c_GF2B_he_primeprime(nBas,Ainter,Binter,Cinter,G_ao1,G_ao2,G_ao3,ERI_AO,Sigma_c_minus(1:nBas,1:nBas)) 
    ! Sigma_eh_2prime
    G_ao1(:,:) = G_ao_itau(2*itau  ,nBas+1:nBas_twice,nBas+1:nBas_twice)  ! eh  -itau
    G_ao2(:,:) = G_ao_itau(2*itau  ,nBas+1:nBas_twice,nBas+1:nBas_twice)  ! eh  -itau
    G_ao3(:,:) = G_ao_itau(2*itau-1,nBas+1:nBas_twice,nBas+1:nBas_twice)  ! eh   itau
    Ainter=czero;Binter=czero;Cinter=czero;
    call Sigma_c_GF2B_he_prime(nBas,Ainter,Binter,Cinter,G_ao1,G_ao2,G_ao3,ERI_AO,Sigma_c_minus(nBas+1:nBas_twice,nBas+1:nBas_twice)) ! Valid for real orbitals 
    ! Sigma_eh_2primeprime
    G_ao1(:,:) = G_ao_itau(2*itau  ,nBas+1:nBas_twice,nBas+1:nBas_twice)  ! eh  -itau
    G_ao2(:,:) = G_ao_itau(2*itau  ,           1:nBas,nBas+1:nBas_twice)  ! hh  -itau
    G_ao3(:,:) = G_ao_itau(2*itau-1,nBas+1:nBas_twice,           1:nBas)  ! ee   itau
    Ainter=czero;Binter=czero;Cinter=czero;
    call Sigma_c_GF2B_eh_primeprime(nBas,Ainter,Binter,Cinter,G_ao1,G_ao2,G_ao3,ERI_AO,Sigma_c_minus(nBas+1:nBas_twice,nBas+1:nBas_twice))
    ! Sigma_hh_2prime
    G_ao1(:,:) = G_ao_itau(2*itau  ,           1:nBas,nBas+1:nBas_twice)  ! hh  -itau
    G_ao2(:,:) = G_ao_itau(2*itau  ,           1:nBas,           1:nBas)  ! he  -itau
    G_ao3(:,:) = G_ao_itau(2*itau-1,           1:nBas,           1:nBas)  ! he   itau
    Ainter=czero;Binter=czero;Cinter=czero;
    call Sigma_c_GF2B_hh_prime(nBas,Ainter,Binter,Cinter,G_ao1,G_ao2,G_ao3,ERI_AO,Sigma_c_minus(1:nBas,nBas+1:nBas_twice))
    ! Sigma_hh_2primeprime
    G_ao1(:,:) = G_ao_itau(2*itau  ,           1:nBas,nBas+1:nBas_twice)  ! hh  -itau
    G_ao2(:,:) = G_ao_itau(2*itau  ,           1:nBas,nBas+1:nBas_twice)  ! hh  -itau
    G_ao3(:,:) = G_ao_itau(2*itau-1,nBas+1:nBas_twice,           1:nBas)  ! ee   itau
    Ainter=czero;Binter=czero;Cinter=czero;
    call Sigma_c_GF2B_hh_primeprime(nBas,Ainter,Binter,Cinter,G_ao1,G_ao2,G_ao3,ERI_AO,Sigma_c_minus(1:nBas,nBas+1:nBas_twice))
    ! Sigma_ee_2prime
    G_ao1(:,:) = G_ao_itau(2*itau  ,nBas+1:nBas_twice,           1:nbas)  ! ee  -itau
    G_ao2(:,:) = G_ao_itau(2*itau  ,           1:nBas,           1:nBas)  ! he  -itau
    G_ao3(:,:) = G_ao_itau(2*itau-1,           1:nBas,           1:nBas)  ! he   itau
    Ainter=czero;Binter=czero;Cinter=czero;
    call Sigma_c_GF2B_ee_prime(nBas,Ainter,Binter,Cinter,G_ao1,G_ao2,G_ao3,ERI_AO,Sigma_c_minus(nBas+1:nBas_twice,1:nBas))
    ! Sigma_ee_2primeprime
    G_ao1(:,:) = G_ao_itau(2*itau  ,nBas+1:nBas_twice,           1:nBas)  ! ee  -itau
    G_ao2(:,:) = G_ao_itau(2*itau  ,           1:nBas,nBas+1:nBas_twice)  ! hh  -itau
    G_ao3(:,:) = G_ao_itau(2*itau-1,nBas+1:nBas_twice,           1:nBas)  ! ee   itau
    Ainter=czero;Binter=czero;Cinter=czero;
    call Sigma_c_GF2B_ee_primeprime(nBas,Ainter,Binter,Cinter,G_ao1,G_ao2,G_ao3,ERI_AO,Sigma_c_minus(nBas+1:nBas_twice,1:nBas))

   ! Corrected Eqs. 17 and 18 in PRB, 109, 245101 (2024)
   Sigma_c_c= -im*(Sigma_c_plus+Sigma_c_minus)
   Sigma_c_s= -   (Sigma_c_plus-Sigma_c_minus)


   ! Corrected Eqs. 17 and 18 in PRB, 109, 245101 (2024)
   Sigma_c_c= -im*(Sigma_c_plus+Sigma_c_minus)
   Sigma_c_s= -   (Sigma_c_plus-Sigma_c_minus)
   ! Sigma_c(i tau) -> Sigma_c(i w)
   do ifreq=1,nfreqs
    Sigma_c_w_ao(ifreq,:,:) = Sigma_c_w_ao(ifreq,:,:)                        &
                            + 0.5d0*cost2w_weight(ifreq,itau)*Sigma_c_c(:,:) &
                            + 0.5d0*sint2w_weight(ifreq,itau)*Sigma_c_s(:,:)
   enddo
   ! Galitskii-Migdal energy [ PRB, 80, 041103R (2009) ]
    ! tau > 0
    Mat_gorkov_tmp(:,:) =G_ao_itau(2*itau-1,:,:) ! G_ao_itau (+)
    Mat_gorkov_tmp=matmul(Sigma_c_minus,Mat_gorkov_tmp)
    do ibas=1,nBas
     EcGM_itau=EcGM_itau+tweight(itau)*Mat_gorkov_tmp(ibas,ibas)
    enddo
    ! tau < 0
    Mat_gorkov_tmp(:,:) =G_ao_itau(2*itau,:,:)   ! G_ao_itau (-)
    Mat_gorkov_tmp=matmul(Sigma_c_plus,Mat_gorkov_tmp)
    do ibas=1,nBas
     EcGM_itau=EcGM_itau+tweight(itau)*Mat_gorkov_tmp(ibas,ibas)
    enddo
  enddo
  EcGM=-real(EcGM_itau) ! Including a factor 2 to sum over spin-channels  EcGM = - 1/2 \sum_spin \int Tr[ Sigma_c_spin(-it) G_spin(it) ]^he dt
                        !                                                      = - \int Tr[ Sigma_c_up(-it) G_up(it) ]^he dt for QP. restricted calcs.


  ! Converge with respect to the H_HFB operator (using only good R_ao matrices -> Tr[R_ao_block S_ao]=Nelectrons )
  if(.not.no_h_hfb) then ! Skiiping the opt w.r.t. the H_HFB operator to do later the linearized approximation on Go -> [ lin-G = Go + Go Sigma Go ]
   iter_hfb=0
   n_diisR=0
   rcondR=0d0
   err_diisR=0d0
   R_ao_old_diis=0d0
   do
    ! Build H_HFB
    iter_hfb=iter_hfb+1
    H_ao_hfb=0d0
    H_ao_hfb(1:nBas,1:nBas)=Hc(1:nBas,1:nBas)
    Ehfbl=0d0
    Ecore=0d0; Eh=0d0; Ex=0d0; Epair=0d0;
    do ibas=1,nBas
     do jbas=1,nBas
      obas=nBas+1+(jbas-1)
      Ehfbl=Ehfbl+2d0*R_ao(ibas,jbas)*Hc(ibas,jbas)
      Ecore=Ecore+2d0*R_ao(ibas,jbas)*Hc(ibas,jbas)
      do kbas=1,nBas
       do lbas=1,nBas
        qbas=nBas+1+(lbas-1)
        H_ao_hfb(ibas,jbas)=H_ao_hfb(ibas,jbas)+2d0*R_ao(kbas,lbas)*vMAT(1+(lbas-1)+(kbas-1)*nBas,1+(jbas-1)+(ibas-1)*nBas) &
                           -R_ao(kbas,lbas)*vMAT(1+(jbas-1)+(kbas-1)*nBas,1+(lbas-1)+(ibas-1)*nBas)
        H_ao_hfb(ibas,obas)=H_ao_hfb(ibas,obas)+sigma*R_ao(kbas,qbas)*vMAT(1+(ibas-1)+(kbas-1)*nBas,1+(jbas-1)+(lbas-1)*nBas)
        Ehfbl=Ehfbl+2d0*R_ao(kbas,lbas)*R_ao(ibas,jbas)*vMAT(1+(lbas-1)+(kbas-1)*nBas,1+(jbas-1)+(ibas-1)*nBas) &
             -R_ao(kbas,lbas)*R_ao(ibas,jbas)*vMAT(1+(jbas-1)+(kbas-1)*nBas,1+(lbas-1)+(ibas-1)*nBas)           &
             +sigma*R_ao(kbas,qbas)*R_ao(ibas,obas)*vMAT(1+(ibas-1)+(kbas-1)*nBas,1+(jbas-1)+(lbas-1)*nBas)
        Eh=Eh+2d0*R_ao(kbas,lbas)*R_ao(ibas,jbas)*vMAT(1+(lbas-1)+(kbas-1)*nBas,1+(jbas-1)+(ibas-1)*nBas)
        Ex=Ex-R_ao(kbas,lbas)*R_ao(ibas,jbas)*vMAT(1+(jbas-1)+(kbas-1)*nBas,1+(lbas-1)+(ibas-1)*nBas) 
        Epair=Epair+sigma*R_ao(kbas,qbas)*R_ao(ibas,obas)*vMAT(1+(ibas-1)+(kbas-1)*nBas,1+(jbas-1)+(lbas-1)*nBas)
       enddo
      enddo
     enddo
    enddo
    H_ao_hfb(nBas+1:nBas_twice,nBas+1:nBas_twice) = -H_ao_hfb(1:nBas,1:nBas           )
    H_ao_hfb(nBas+1:nBas_twice,1:nBas           ) =  H_ao_hfb(1:nBas,nBas+1:nBas_twice)
    ! Build G(i w) and R
    R_ao_old=R_ao
    call get_1rdm_scGWB(nBas,nBas_twice,nfreqs,chem_pot,S,H_ao_hfb,Sigma_c_w_ao,wcoord,wweight, &
                        Mat_gorkov_tmp,G_ao_iw_hfb,DeltaG_ao_iw,R_ao,R_ao_hfb,trace_1_rdm) 
    if(abs(trace_1_rdm-nElectrons)**2d0>thrs_N .and. chem_pot_scG) &
     call fix_chem_pot_scGWB_bisec(iter_hfb,nBas,nBas_twice,nfreqs,nElectrons,thrs_N,thrs_Ngrad,chem_pot,S,H_ao_hfb,Sigma_c_w_ao,   &
                                   wcoord,wweight,Mat_gorkov_tmp,G_ao_iw_hfb,DeltaG_ao_iw,R_ao,R_ao_hfb,trace_1_rdm,chem_pot_saved, &
                                   verbose_scGF2B)
    if(abs(trace_1_rdm-nElectrons)**2d0>thrs_N .and. .not.chem_pot_scG) &
     R_ao=nElectrons*R_ao/trace_1_rdm

    ! Check convergence of R_ao for fixed Sigma_c(i w)
    diff_Rao=0d0
    do ibas=1,nBas_twice
     do jbas=1,nBas_twice
      diff_Rao=diff_Rao+abs(R_ao(ibas,jbas)-R_ao_old(ibas,jbas))
     enddo
    enddo

    if(diff_Rao<=thrs_Rao) exit

    if(iter_hfb==maxSCF) exit

    ! Do mixing with previous R_ao to facilitate convergence
    if(maxDIIS>0) then
     n_diisR=min(n_diisR+1,maxDIIS)
     err_currentR=0d0
     idiis_indexR=1
     do ibas=1,nBas_twice
      do jbas=1,nBas_twice
       err_currentR(idiis_indexR)=R_ao(ibas,jbas)-R_ao_old(ibas,jbas)
       if(abs(err_currentR(idiis_indexR))<1e-12) err_currentR(idiis_indexR)=0d0
       R_ao_extrap(idiis_indexR)=R_ao(ibas,jbas)
       idiis_indexR=idiis_indexR+1
      enddo
     enddo
     call DIIS_extrapolation(rcondR,nBas_twiceSq,nBas_twiceSq,n_diisR,err_diisR,R_ao_old_diis,err_currentR,R_ao_extrap)
     idiis_indexR=1
     do ibas=1,nBas_twice
      do jbas=1,nBas_twice
       R_ao(ibas,jbas)=R_ao_extrap(idiis_indexR)
       idiis_indexR=idiis_indexR+1
      enddo
     enddo
    else
     R_ao(:,:)=alpha_mixing*R_ao(:,:)+(1d0-alpha_mixing)*R_ao_old(:,:)
    endif

   enddo
  endif

  ! Check convergence of R_ao after a scGF2B iteration
  diff_Rao=0d0
  do ibas=1,nBas_twice
   do jbas=1,nBas_twice
    diff_Rao=diff_Rao+abs(R_ao(ibas,jbas)-R_ao_iter(ibas,jbas))
   enddo 
  enddo
  R_ao_iter=R_ao

  ! Print iter info
  U_mo=-2d0*matmul(matmul(cHFBinv,R_ao(1:nBas,1:nBas)),transpose(cHFBinv)) ! Minus to order occ numbers
  call diagonalize_matrix(nOrb,U_mo,Occ)
  Occ=-Occ
  trace_1_rdm=sum(Occ)
  cNO=matmul(cHFB,U_mo)
  write(*,*)
  write(*,'(a,f15.8,a,i5,a,i5)') ' Trace scGF2B ',trace_1_rdm,' after ',iter_hfb,' HFB iterations at global iter ',iter
  write(*,'(a,f15.8)')        ' Change of R  ',diff_Rao
  write(*,'(a,f15.8)')        ' Chem. Pot.   ',chem_pot
  write(*,'(a,f15.8)')        ' Enuc         ',ENuc
  write(*,'(a,f15.8)')        ' Ehcore       ',Ecore
  write(*,'(a,f15.8)')        ' Hartree      ',Eh
  write(*,'(a,f15.8)')        ' Exchange     ',Ex
  write(*,'(a,f15.8)')        ' Epairing     ',Epair
  write(*,'(a,f15.8)')        ' Ehfbl        ',Ehfbl
  write(*,'(a,f15.8)')        ' EcGM         ',EcGM
  write(*,'(a,f15.8)')        ' Eelec        ',Ehfbl+EcGM
  write(*,'(a,f15.8)')        ' Etot         ',Ehfbl+EcGM+ENuc
  write(*,*)

  if(diff_Rao<=thrs_Rao) exit

  if(iter==maxSCF) exit

  ! Transform DeltaG(i w) -> DeltaG(i tau) [ i tau and -i tau ]
  !      [ the weights contain the 2 /(2 pi) = 1 / pi factor and the cos(tau w) or sin(tau w) ]
  G_ao_itau=czero
  do itau=1,ntimes
   Mat_gorkov_tmp(:,:)=czero
   Mat_gorkov_tmp2(:,:)=czero
   do ifreq=1,nfreqs
    Mat_gorkov_tmp(:,:) = Mat_gorkov_tmp(:,:)   + im*cosw2t_weight(itau,ifreq)*Real(DeltaG_ao_iw(ifreq,:,:))  &
                                                - im*sinw2t_weight(itau,ifreq)*Aimag(DeltaG_ao_iw(ifreq,:,:))
    Mat_gorkov_tmp2(:,:) = Mat_gorkov_tmp2(:,:) + im*cosw2t_weight(itau,ifreq)*Real(DeltaG_ao_iw(ifreq,:,:))  &
                                                + im*sinw2t_weight(itau,ifreq)*Aimag(DeltaG_ao_iw(ifreq,:,:))
   enddo
   ! Build G(i tau) = DeltaG(i tau) + Go(i tau)
   G_ao_itau(2*itau-1,:,:)=Mat_gorkov_tmp(:,:) +G_ao_itau_hfb(2*itau-1,:,:)
   G_ao_itau(2*itau  ,:,:)=Mat_gorkov_tmp2(:,:)+G_ao_itau_hfb(2*itau  ,:,:)
  enddo

  ! Do mixing with previous G(i tau) to facilitate convergence
  if(maxDIIS>0) then
   n_diis=min(n_diis+1,maxDIIS)
   err_current=czero
   idiis_index=1
   do itau=1,ntimes_twice
    do ibas=1,nBas_twice
     do jbas=1,nBas_twice
      err_current(idiis_index)=G_ao_itau(itau,ibas,jbas)-G_ao_itau_old(itau,ibas,jbas)
      if(abs(err_current(idiis_index))<1e-12) err_current(idiis_index)=czero
      G_itau_extrap(idiis_index)=G_ao_itau(itau,ibas,jbas)
      idiis_index=idiis_index+1
     enddo
    enddo
   enddo
   call complex_DIIS_extrapolation(rcond,nBas2Sqntimes2,nBas2Sqntimes2,n_diis,err_diis,G_itau_old_diis,err_current,G_itau_extrap)
   idiis_index=1
   do itau=1,ntimes_twice
    do ibas=1,nBas_twice
     do jbas=1,nBas_twice
      G_ao_itau(itau,ibas,jbas)=G_itau_extrap(idiis_index)
      idiis_index=idiis_index+1
     enddo
    enddo
   enddo
  else
   G_ao_itau(:,:,:)=alpha_mixing*G_ao_itau(:,:,:)+(1d0-alpha_mixing)*G_ao_itau_old(:,:,:)
  endif
  G_ao_itau_old(:,:,:)=G_ao_itau(:,:,:)

 enddo
 N_anom = trace_matrix(nBas,matmul(transpose(R_ao(1:nBas,nBas+1:nBas_twice)), &
          R_ao(1:nBas,nBas+1:nBas_twice)))
 write(*,*)
 write(*,'(A50)') '---------------------------------------'
 write(*,'(A50)') '     scGF2B calculation completed      '
 write(*,'(A50)') '---------------------------------------'
 write(*,*)
 write(*,'(a,f15.8,a,i5,a)') ' Trace scGF2B   ',trace_1_rdm,' after ',iter,' global iterations '
 write(*,'(a,f15.8)')        ' Change of R    ',diff_Rao
 write(*,'(a,f15.8)')        ' Chem. Pot.     ',chem_pot
 write(*,'(a,f15.8)')        ' N anomalus     ',N_anom
 write(*,'(a,f15.8)')        ' Enuc           ',ENuc
 write(*,'(a,f15.8)')        ' Ehcore         ',Ecore
 write(*,'(a,f15.8)')        ' Hartree        ',Eh
 write(*,'(a,f15.8)')        ' Exchange       ',Ex
 write(*,'(a,f15.8)')        ' Epairing       ',Epair
 write(*,'(a,f15.8)')        ' Ehfbl          ',Ehfbl
 write(*,'(a,f15.8)')        ' EcGM           ',EcGM
 write(*,'(a,f15.8)')        ' Eelec          ',Ehfbl+EcGM
 write(*,'(a,f15.8)')        ' scGF2B Energy  ',Ehfbl+EcGM+ENuc
 write(*,*)
 write(*,'(a,f15.8)')        ' EcPT2          ',EcGM/2d0
 write(*,*)
 write(*,*) ' Final occupation numbers'
 do ibas=1,nOrb
  write(*,'(I7,F15.8)') ibas,Occ(ibas)
 enddo
 if(verbose/=0) then
  write(*,*) ' Natural orbitals (columns)'
  do ibas=1,nBas
   write(*,'(*(f15.8))') cNO(ibas,:)
  enddo
 endif
 write(*,*)

 ! Write restart files
 call write_scGWB_restart(nBas_twice,ntimes,ntimes_twice,nfreqs,chem_pot,R_ao,R_ao_hfb,G_ao_itau,G_ao_itau_hfb, &
                         G_ao_iw_hfb,DeltaG_ao_iw)

 inquire(file='Print_Rao', exist=file_exists)
 if(file_exists) then
  write(*,*) 'R_scGF2B_ao'
  do ibas=1,nBas_twice
   write(*,'(*(f10.5))') R_ao(ibas,:)
  enddo
 endif

 ! Deallocate arrays
 deallocate(Occ)
 deallocate(cHFBinv)
 deallocate(U_QP_tmp)
 deallocate(cHFB_gorkov)
 deallocate(cNO)
 deallocate(U_mo)
 deallocate(R_ao)
 deallocate(R_ao_iter)
 deallocate(R_ao_hfb)
 deallocate(R_ao_old)
 deallocate(R_mo)
 deallocate(H_ao_hfb)
 deallocate(G_ao_tmp)
 deallocate(Mat_gorkov_tmp)
 deallocate(Mat_gorkov_tmp2)
 deallocate(DeltaG_ao_iw)
 deallocate(G_ao_iw_hfb)
 deallocate(G_ao_itau_hfb)
 deallocate(G_ao_itau)
 deallocate(G_ao_itau_old)
 deallocate(G_ao1)
 deallocate(G_ao2)
 deallocate(G_ao3)
 deallocate(Ainter)
 deallocate(Binter)
 deallocate(Cinter)
 deallocate(Sigma_c_w_ao)
 deallocate(Sigma_c_plus)
 deallocate(Sigma_c_minus)
 deallocate(Sigma_c_c)
 deallocate(Sigma_c_s)
 deallocate(tweight,tcoord)
 deallocate(sint2w_weight)
 deallocate(cost2w_weight)
 deallocate(cosw2t_weight)
 deallocate(sinw2t_weight)
 deallocate(Mat1,Mat2)
 deallocate(err_currentR)
 deallocate(R_ao_extrap)
 deallocate(err_diisR)
 deallocate(R_ao_old_diis)
 deallocate(G_itau_extrap)
 deallocate(err_current)
 deallocate(err_diis)
 deallocate(G_itau_old_diis)

 call wall_time(end_scGF2Bitauiw)
 
 t_scGF2Bitauiw = end_scGF2Bitauiw - start_scGF2Bitauiw
 write(*,'(A65,1X,F9.3,A8)') 'Total wall time for scGF2B = ',t_scGF2Bitauiw,' seconds'
 write(*,*)

end subroutine 
