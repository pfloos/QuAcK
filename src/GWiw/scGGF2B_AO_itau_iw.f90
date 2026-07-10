subroutine scGGF2B_AO_itau_iw(nBas2,nBas4,nOrb2,nOrb4,maxSCF,thresh_in,maxDIIS,restart_scGF2B,verbose_scGF2B,          &
                              chem_pot_scG,no_h_hfb,ENuc,Gen_Hc,Gen_S,Gen_R_in,Gen_cHFB,Gen_eQP_state,chem_pot,sigma,  &
                              nfreqs,wcoord,wweight,Gen_U_QP,db_ERI_AO)

! Generalized scGF2B

  implicit none
  include 'parameters.h'

! Input variables
 
  logical,intent(in)            :: no_h_hfb
  logical,intent(in)            :: restart_scGF2B
  logical,intent(in)            :: verbose_scGF2B
  logical,intent(in)            :: chem_pot_scG

  integer,intent(in)            :: nBas2
  integer,intent(in)            :: nBas4
  integer,intent(in)            :: nOrb2
  integer,intent(in)            :: nOrb4
  integer,intent(in)            :: maxSCF
  integer,intent(in)            :: maxDIIS

  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: thresh_in
  double precision,intent(in)   :: sigma
  double precision,intent(in)   :: Gen_Hc(nBas2,nBas2)
  double precision,intent(in)   :: Gen_R_in(nBas4,nBas4)
  double precision,intent(in)   :: Gen_S(nBas2,nBas2)
  double precision,intent(in)   :: db_ERI_AO(nBas2,nBas2,nBas2,nBas2)
  double precision,intent(in)   :: Gen_U_QP(nOrb4,nOrb4)

! Local variables

  logical                       :: file_exists
  logical                       :: read_HFB_chkp

  integer                       :: iunit=313
  integer                       :: ifreq,itau
  integer                       :: abas,bbas,cbas,dbas
  integer                       :: pbas,qbas,rbas,sbas
  integer                       :: iorb,jorb,korb,lorb
  integer                       :: iter,iter_hfb
  integer                       :: verbose
  integer                       :: ntimes
  integer                       :: nBas4Sq
  integer                       :: ntimes_twice
  integer                       :: imax_error_gw2gt
  integer                       :: imax_error_st2sw
  integer                       :: idiis_indexR,n_diisR
  integer                       :: idiis_index,n_diis
  integer                       :: nBas4Sqntimes2

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
  double precision              :: EcGM,Ehfbl,Ecore,Ehx,Epair
  double precision              :: trace_1_rdm
  double precision,external     :: trace_matrix
  double precision              :: start_scGF2Bitauiw,end_scGF2Bitauiw,t_scGF2Bitauiw
  double precision,allocatable  :: Occ(:)
  double precision,allocatable  :: Mat1(:,:)
  double precision,allocatable  :: Mat2(:,:)
  double precision,allocatable  :: Mat3(:,:)
  double precision,allocatable  :: Mat4(:,:)
  double precision,allocatable  :: Gen_cHFBinv(:,:)
  double precision,allocatable  :: Gen_cNO(:,:)
  double precision,allocatable  :: U_mo(:,:)
  double precision,allocatable  :: Gen_R_ao(:,:)
  double precision,allocatable  :: Gen_R_ao_iter(:,:)
  double precision,allocatable  :: Gen_R_ao_hfb(:,:)
  double precision,allocatable  :: Gen_R_ao_old(:,:)
  double precision,allocatable  :: Gen_H_ao_hfb(:,:)
  double precision,allocatable  :: err_currentR(:)
  double precision,allocatable  :: err_diisR(:,:)
  double precision,allocatable  :: Gen_R_ao_extrap(:)
  double precision,allocatable  :: Gen_R_ao_old_diis(:,:)
  double precision,allocatable  :: wcoord_int(:),wweight_int(:)
  double precision,allocatable  :: tweight(:),tcoord(:)
  double precision,allocatable  :: sint2w_weight(:,:)
  double precision,allocatable  :: cost2w_weight(:,:)
  double precision,allocatable  :: cosw2t_weight(:,:)
  double precision,allocatable  :: sinw2t_weight(:,:)

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
  complex*16,allocatable        :: Ainter(:,:,:,:)
  complex*16,allocatable        :: Binter(:,:,:,:)
  complex*16,allocatable        :: Cinter(:,:,:,:)

! Output variables
  integer,intent(inout)         :: nfreqs
  double precision,intent(inout):: chem_pot
  double precision,intent(inout):: Gen_eQP_state(nOrb2)
  double precision,intent(inout):: wcoord(nfreqs)
  double precision,intent(inout):: wweight(nfreqs)
  double precision,intent(inout):: Gen_cHFB(nBas2,nOrb2)
  
 call wall_time(start_scGF2Bitauiw)

 write(*,*)     
 write(*,*)'******************************************'
 write(*,*)'*   scGGF2B ( using it and iw grids )    *'
 write(*,*)'******************************************'
 write(*,*)

 write(*,*)
 write(*,'(A33,1X,F16.10,A3)') ' Initial chemical potential  = ',chem_pot,' au'
 if(chem_pot_scG) then
  write(*,'(A)') '   Adjusting the chemical potential is activated'
 else
  write(*,'(A)') '   Adjusting the chemical potential is deactivated'
 endif
 write(*,*)

 write(*,'(a,F15.8)') '  emin ',abs(Gen_eQP_state(nOrb2))
 write(*,'(a,F15.8)') '  emax ',abs(Gen_eQP_state(nOrb2)+Gen_eQP_state(1))

 ! Initialize variables
 ntimes=nfreqs
 ntimes_twice=ntimes*2
 nBas4Sq=nBas4*nBas4
 verbose=0
 if(verbose_scGF2B) verbose=1     
 eta=0d0
 thrs_N=1d-10
 thrs_Ngrad=1d-6
 inquire(file='threshold_N', exist=file_exists)
 if(file_exists) then
  open(unit=iunit,form='formatted',file='threshold_N',status='old')
  read(iunit,*) thrs_N
  close(iunit)
  write(*,*)
  write(*,'(a,f15.8)') ' using as thrs_N ',thrs_N
  write(*,*)
 endif
 thrs_Rao=thresh_in
 nElectrons=0d0
 do abas=1,nBas2
  do bbas=1,nBas2
   nElectrons=nElectrons+Gen_R_in(abas,bbas)*Gen_S(abas,bbas)
  enddo
 enddo
 nElectrons=nint(nElectrons) ! We use both spin channels
 chem_pot_saved=chem_pot
 alpha_mixing=0.6d0
 rcond=0d0
 n_diis=0
 nBas4Sqntimes2=nBas4Sq*ntimes_twice

 ! Allocate arrays
 allocate(Occ(nOrb2))
 allocate(Gen_cNO(nBas2,nOrb2))
 allocate(U_mo(nOrb2,nOrb2))
 allocate(Gen_cHFBinv(nOrb2,nBas2))
 allocate(Gen_R_ao(nBas4,nBas4))
 allocate(Gen_R_ao_iter(nBas4,nBas4))
 allocate(Gen_R_ao_hfb(nBas4,nBas4))
 allocate(Gen_R_ao_old(nBas4,nBas4))
 allocate(Gen_H_ao_hfb(nBas4,nBas4))
 allocate(G_ao_tmp(nBas2,nBas2))
 allocate(Mat_gorkov_tmp(nBas4,nBas4))
 allocate(Mat_gorkov_tmp2(nBas4,nBas4))
 allocate(G_ao_iw_hfb(nfreqs,nBas4,nBas4))
 allocate(DeltaG_ao_iw(nfreqs,nBas4,nBas4))
 allocate(G_ao_itau_hfb(ntimes_twice,nBas4,nBas4))
 allocate(G_ao_itau(ntimes_twice,nBas4,nBas4))
 allocate(G_ao_itau_old(ntimes_twice,nBas4,nBas4))
 allocate(Sigma_c_w_ao(nfreqs,nBas4,nBas4))
 allocate(Sigma_c_plus(nBas4,nBas4)) 
 allocate(Sigma_c_minus(nBas4,nBas4)) 
 allocate(Sigma_c_c(nBas4,nBas4))      
 allocate(Sigma_c_s(nBas4,nBas4))
 allocate(tweight(ntimes),tcoord(ntimes))
 allocate(sint2w_weight(nfreqs,ntimes))
 allocate(cost2w_weight(nfreqs,ntimes))
 allocate(cosw2t_weight(ntimes,nfreqs))
 allocate(sinw2t_weight(ntimes,nfreqs))
 allocate(Mat1(nOrb2,nOrb2),Mat2(nOrb2,nOrb2))
 allocate(Mat3(nOrb2,nOrb2),Mat4(nOrb2,nOrb2))
 allocate(Ainter(nBas2,nBas2,nBas2,nBas2))
 allocate(Binter(nBas2,nBas2,nBas2,nBas2))
 allocate(Cinter(nBas2,nBas2,nBas2,nBas2))
 allocate(G_ao1(nBas2,nBas2))
 allocate(G_ao2(nBas2,nBas2))
 allocate(G_ao3(nBas2,nBas2))
 allocate(err_currentR(1))
 allocate(err_diisR(1,1))
 allocate(Gen_R_ao_extrap(1))
 allocate(Gen_R_ao_old_diis(1,1))
 allocate(G_itau_extrap(1))
 allocate(err_current(1))
 allocate(err_diis(1,1))
 allocate(G_itau_old_diis(1,1))
 if(maxDIIS>0) then
  deallocate(err_currentR)
  deallocate(Gen_R_ao_extrap)
  deallocate(err_diisR)
  deallocate(Gen_R_ao_old_diis)
  deallocate(G_itau_extrap)
  deallocate(err_current)
  deallocate(err_diis)
  deallocate(G_itau_old_diis)
  allocate(err_currentR(nBas4Sq))
  allocate(Gen_R_ao_extrap(nBas4Sq))
  allocate(err_diisR(nBas4Sq,maxDIIS))
  allocate(Gen_R_ao_old_diis(nBas4Sq,maxDIIS))
  allocate(G_itau_extrap(nBas4Sqntimes2))
  allocate(err_current(nBas4Sqntimes2))
  allocate(err_diis(nBas4Sqntimes2,maxDIIS))
  allocate(G_itau_old_diis(nBas4Sqntimes2,maxDIIS))
 endif
 err_diis=czero
 G_itau_old_diis=czero

 ! Initialize arrays
 DeltaG_ao_iw(:,:,:)=czero
 Mat1(1:nOrb2,1:nOrb2)=Gen_U_QP(1:nOrb2      ,1:nOrb2)        ! \bar{V}
 Mat2(1:nOrb2,1:nOrb2)=Gen_U_QP(nOrb2+1:nOrb4,1:nOrb2)        ! \bar{U}
 do iorb=1,nOrb2
  Mat3(1:nOrb2,iorb)=Gen_U_QP(1:nOrb2      ,nOrb4-(iorb-1))   ! U
  Mat4(1:nOrb2,iorb)=Gen_U_QP(nOrb2+1:nOrb4,nOrb4-(iorb-1))   ! V
 enddo
 Gen_R_ao_hfb=Gen_R_in
 Gen_R_ao=Gen_R_ao_hfb
 Gen_R_ao_iter=Gen_R_ao_hfb
 Gen_cHFBinv=matmul(transpose(Gen_cHFB),Gen_S)
 Gen_H_ao_hfb=0d0
 Gen_H_ao_hfb(1:nBas2      ,1:nBas2      ) =   Gen_Hc(1:nBas2,1:nBas2) - chem_pot*Gen_S(1:nBas2,1:nBas2)
 Gen_H_ao_hfb(nBas2+1:nBas4,nBas2+1:nBas4) = -(Gen_Hc(1:nBas2,1:nBas2) - chem_pot*Gen_S(1:nBas2,1:nBas2))
 do abas=1,nBas2
  pbas=abas+nBas2
  do bbas=1,nBas2
   qbas=bbas+nBas2
   do cbas=1,nBas2
    rbas=cbas+nBas2
    do dbas=1,nBas2
     sbas=dbas+nBas2
     Gen_H_ao_hfb(abas,bbas)=Gen_H_ao_hfb(abas,bbas)+Gen_R_ao(cbas,dbas)*db_ERI_AO(abas,cbas,bbas,dbas) 
     Gen_H_ao_hfb(pbas,qbas)=Gen_H_ao_hfb(pbas,qbas)-Gen_R_ao(cbas,dbas)*db_ERI_AO(bbas,dbas,abas,cbas) 
     Gen_H_ao_hfb(abas,qbas)=Gen_H_ao_hfb(abas,qbas)+0.5d0*sigma*Gen_R_ao(cbas,sbas)*db_ERI_AO(abas,bbas,cbas,dbas)
     Gen_H_ao_hfb(pbas,bbas)=Gen_H_ao_hfb(pbas,bbas)+0.5d0*sigma*Gen_R_ao(rbas,dbas)*db_ERI_AO(abas,bbas,cbas,dbas)
    enddo
   enddo
  enddo
 enddo

 ! Read grids 
 call read_scGX_grids(ntimes,nfreqs,tcoord,tweight,wcoord,wweight,sint2w_weight,cost2w_weight, &
                      cosw2t_weight,sinw2t_weight,verbose)

 ! Build Go(i w)
 do ifreq=1,nfreqs
  weval_cpx=im*wcoord(ifreq)
  ! G_he(iw)
  call G_AO_GHFB_w(nBas2,nOrb2,nOrb4,eta,Gen_cHFB,Gen_eQP_state,weval_cpx,Mat1,Mat1,Mat3,Mat3,G_ao_tmp)
  G_ao_iw_hfb(ifreq,1:nBas2      ,1:nBas2      ) = G_ao_tmp(1:nBas2,1:nBas2)
  ! G_hh(iw)
  call G_AO_GHFB_w(nBas2,nOrb2,nOrb4,eta,Gen_cHFB,Gen_eQP_state,weval_cpx,Mat1,Mat2,Mat3,Mat4,G_ao_tmp)
  G_ao_iw_hfb(ifreq,1:nBas2      ,nBas2+1:nBas4) = G_ao_tmp(1:nBas2,1:nBas2)
  ! G_ee(iw)
  call G_AO_GHFB_w(nBas2,nOrb2,nOrb4,eta,Gen_cHFB,Gen_eQP_state,weval_cpx,Mat2,Mat1,Mat4,Mat3,G_ao_tmp)
  G_ao_iw_hfb(ifreq,nBas2+1:nBas4,1:nBas2      ) = G_ao_tmp(1:nBas2,1:nBas2)
  ! G_eh(iw)
  call G_AO_GHFB_w(nBas2,nOrb2,nOrb4,eta,Gen_cHFB,Gen_eQP_state,weval_cpx,Mat2,Mat2,Mat4,Mat4,G_ao_tmp)
  G_ao_iw_hfb(ifreq,nBas2+1:nBas4,nBas2+1:nBas4) = G_ao_tmp(1:nBas2,1:nBas2)
 enddo
 ! Build Go(i tau)
 do itau=1,ntimes
  ! tau > 0
  ! G_he(i tau)
  call G_AO_GHFB_itau(nBas2,nOrb2,nOrb4, tcoord(itau),G_ao_tmp,Gen_cHFB,Gen_eQP_state,Mat1,Mat1,Mat3,Mat3)
  G_ao_itau_hfb(2*itau-1,1:nBas2      ,1:nBas2      ) = G_ao_tmp(1:nBas2,1:nBas2)
  ! G_hh(i tau)
  call G_AO_GHFB_itau(nBas2,nOrb2,nOrb4, tcoord(itau),G_ao_tmp,Gen_cHFB,Gen_eQP_state,Mat1,Mat2,Mat3,Mat4)
  G_ao_itau_hfb(2*itau-1,1:nBas2      ,nBas2+1:nBas4) = G_ao_tmp(1:nBas2,1:nBas2)
  ! G_ee(i tau)
  call G_AO_GHFB_itau(nBas2,nOrb2,nOrb4, tcoord(itau),G_ao_tmp,Gen_cHFB,Gen_eQP_state,Mat2,Mat1,Mat4,Mat3)
  G_ao_itau_hfb(2*itau-1,nBas2+1:nBas4,1:nBas2      ) = G_ao_tmp(1:nBas2,1:nBas2)
  ! G_eh(i tau)
  call G_AO_GHFB_itau(nBas2,nOrb2,nOrb4, tcoord(itau),G_ao_tmp,Gen_cHFB,Gen_eQP_state,Mat2,Mat2,Mat4,Mat4)
  G_ao_itau_hfb(2*itau-1,nBas2+1:nBas4,nBas2+1:nBas4) = G_ao_tmp(1:nBas2,1:nBas2)
  ! tau < 0
  ! G_he(i tau)
  call G_AO_GHFB_itau(nBas2,nOrb2,nOrb4,-tcoord(itau),G_ao_tmp,Gen_cHFB,Gen_eQP_state,Mat1,Mat1,Mat3,Mat3)
  G_ao_itau_hfb(2*itau  ,1:nBas2      ,1:nBas2      ) = G_ao_tmp(1:nBas2,1:nBas2)
  ! G_hh(i tau)
  call G_AO_GHFB_itau(nBas2,nOrb2,nOrb4,-tcoord(itau),G_ao_tmp,Gen_cHFB,Gen_eQP_state,Mat1,Mat2,Mat3,Mat4)
  G_ao_itau_hfb(2*itau  ,1:nBas2      ,nBas2+1:nBas4) = G_ao_tmp(1:nBas2,1:nBas2)
  ! G_ee(i tau)
  call G_AO_GHFB_itau(nBas2,nOrb2,nOrb4,-tcoord(itau),G_ao_tmp,Gen_cHFB,Gen_eQP_state,Mat2,Mat1,Mat4,Mat3)
  G_ao_itau_hfb(2*itau  ,nBas2+1:nBas4,1:nBas2      ) = G_ao_tmp(1:nBas2,1:nBas2)
  ! G_eh(i tau)
  call G_AO_GHFB_itau(nBas2,nOrb2,nOrb4,-tcoord(itau),G_ao_tmp,Gen_cHFB,Gen_eQP_state,Mat2,Mat2,Mat4,Mat4)
  G_ao_itau_hfb(2*itau  ,nBas2+1:nBas4,nBas2+1:nBas4) = G_ao_tmp(1:nBas2,1:nBas2)
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
 write(*,'(a,*(f20.8))') ' MAE       ',sum_error_gw2gt/(2*ntimes*nBas4*nBas4)


 ! If required, read restart files
 if(restart_scGF2B) then
  inquire(file='read_HFB_scGWB', exist=file_exists)
  read_HFB_chkp=.false.
  if(file_exists) read_HFB_chkp=.true.
  call read_scGXB_restart(nBas4,nfreqs,ntimes_twice,chem_pot,Gen_R_ao,Gen_R_ao_hfb,G_ao_iw_hfb,G_ao_itau,G_ao_itau_hfb,read_HFB_chkp)
  Gen_R_ao_iter=Gen_R_ao
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
    G_ao1(1:nBas2,1:nBas2) = G_ao_itau(2*itau-1,      1:nBas2,      1:nBas2)  ! he   itau
    G_ao2(1:nBas2,1:nBas2) = G_ao_itau(2*itau-1,      1:nBas2,      1:nBas2)  ! he   itau
    G_ao3(1:nBas2,1:nBas2) = G_ao_itau(2*itau  ,      1:nBas2,      1:nBas2)  ! he  -itau
    Ainter=czero;Binter=czero;Cinter=czero;
    call Sigma_c_GGF2B_he_prime(nBas2,Ainter,Binter,Cinter,G_ao1,G_ao2,G_ao3,db_ERI_AO,Sigma_c_plus(1:nBas2,1:nBas2))
    ! Sigma_he_2primeprime
    G_ao1(1:nBas2,1:nBas2) = G_ao_itau(2*itau-1,      1:nBas2,      1:nBas2)  ! he   itau
    G_ao2(1:nBas2,1:nBas2) = G_ao_itau(2*itau-1,      1:nBas2,nBas2+1:nBas4)  ! hh   itau
    G_ao3(1:nBas2,1:nBas2) = G_ao_itau(2*itau  ,nBas2+1:nBas4,      1:nBas2)  ! ee  -itau
    Ainter=czero;Binter=czero;Cinter=czero;
    call Sigma_c_GGF2B_he_prime2(nBas2,Ainter,Binter,Cinter,G_ao1,G_ao2,G_ao3,db_ERI_AO,Sigma_c_plus(1:nBas2,1:nBas2))
    ! Sigma_eh_2prime
    G_ao1(1:nBas2,1:nBas2) = G_ao_itau(2*itau-1,nBas2+1:nBas4,nBas2+1:nBas4)  ! eh   itau
    G_ao2(1:nBas2,1:nBas2) = G_ao_itau(2*itau-1,nBas2+1:nBas4,nBas2+1:nBas4)  ! eh   itau
    G_ao3(1:nBas2,1:nBas2) = G_ao_itau(2*itau  ,nBas2+1:nBas4,nBas2+1:nBas4)  ! eh  -itau
    Ainter=czero;Binter=czero;Cinter=czero;
    call Sigma_c_GGF2B_eh_prime(nBas2,Ainter,Binter,Cinter,G_ao1,G_ao2,G_ao3,db_ERI_AO,Sigma_c_plus(nBas2+1:nBas4,nBas2+1:nBas4))
    ! Sigma_eh_2primeprime
    G_ao1(1:nBas2,1:nBas2) = G_ao_itau(2*itau-1,nBas2+1:nBas4,nBas2+1:nBas4)  ! eh   itau
    G_ao2(1:nBas2,1:nBas2) = G_ao_itau(2*itau-1,      1:nBas2,nBas2+1:nBas4)  ! hh   itau
    G_ao3(1:nBas2,1:nBas2) = G_ao_itau(2*itau  ,nBas2+1:nBas4,      1:nBas2)  ! ee  -itau
    Ainter=czero;Binter=czero;Cinter=czero;
    call Sigma_c_GGF2B_eh_prime2(nBas2,Ainter,Binter,Cinter,G_ao1,G_ao2,G_ao3,db_ERI_AO,Sigma_c_plus(nBas2+1:nBas4,nBas2+1:nBas4))
    ! Sigma_hh_2prime
    G_ao1(1:nBas2,1:nBas2) = G_ao_itau(2*itau-1,      1:nBas2,nBas2+1:nBas4)  ! hh   itau
    G_ao2(1:nBas2,1:nBas2) = G_ao_itau(2*itau-1,      1:nBas2,      1:nBas2)  ! he   itau
    G_ao3(1:nBas2,1:nBas2) = G_ao_itau(2*itau  ,      1:nBas2,      1:nBas2)  ! he  -itau
    Ainter=czero;Binter=czero;Cinter=czero;
    call Sigma_c_GGF2B_hh_prime(nBas2,Ainter,Binter,Cinter,G_ao1,G_ao2,G_ao3,db_ERI_AO,Sigma_c_plus(1:nBas2,nBas2+1:nBas4))
    ! Sigma_hh_2primeprime
    G_ao1(1:nBas2,1:nBas2) = G_ao_itau(2*itau-1,      1:nBas2,nBas2+1:nBas4)  ! hh   itau
    G_ao2(1:nBas2,1:nBas2) = G_ao_itau(2*itau-1,      1:nBas2,nBas2+1:nBas4)  ! hh   itau
    G_ao3(1:nBas2,1:nBas2) = G_ao_itau(2*itau  ,nBas2+1:nBas4,      1:nBas2)  ! ee  -itau
    Ainter=czero;Binter=czero;Cinter=czero;
    call Sigma_c_GGF2B_hh_prime2(nBas2,Ainter,Binter,Cinter,G_ao1,G_ao2,G_ao3,db_ERI_AO,Sigma_c_plus(1:nBas2,nBas2+1:nBas4))
    ! Sigma_ee_2prime
    G_ao1(1:nBas2,1:nBas2) = G_ao_itau(2*itau-1,nBas2+1:nBas4,      1:nBas2)  ! ee   itau
    G_ao2(1:nBas2,1:nBas2) = G_ao_itau(2*itau-1,      1:nBas2,      1:nBas2)  ! he   itau
    G_ao3(1:nBas2,1:nBas2) = G_ao_itau(2*itau  ,      1:nBas2,      1:nBas2)  ! he  -itau
    Ainter=czero;Binter=czero;Cinter=czero;
    call Sigma_c_GGF2B_ee_prime(nBas2,Ainter,Binter,Cinter,G_ao1,G_ao2,G_ao3,db_ERI_AO,Sigma_c_plus(nBas2+1:nBas4,1:nBas2))
    ! Sigma_ee_2primeprime
    G_ao1(1:nBas2,1:nBas2) = G_ao_itau(2*itau-1,nBas2+1:nBas4,      1:nBas2)  ! ee   itau
    G_ao2(1:nBas2,1:nBas2) = G_ao_itau(2*itau-1,      1:nBas2,nBas2+1:nBas4)  ! hh   itau
    G_ao3(1:nBas2,1:nBas2) = G_ao_itau(2*itau  ,nBas2+1:nBas4,      1:nBas2)  ! ee  -itau
    Ainter=czero;Binter=czero;Cinter=czero;
    call Sigma_c_GGF2B_ee_prime2(nBas2,Ainter,Binter,Cinter,G_ao1,G_ao2,G_ao3,db_ERI_AO,Sigma_c_plus(nBas2+1:nBas4,1:nBas2))

   ! G(-i tau) G(-i tau) G(i tau) -> Sigma_c(-i tau)
    ! Sigma_he_2prime
    G_ao1(1:nBas2,1:nBas2) = G_ao_itau(2*itau  ,      1:nBas2,      1:nBas2)  ! he  -itau
    G_ao2(1:nBas2,1:nBas2) = G_ao_itau(2*itau  ,      1:nBas2,      1:nBas2)  ! he  -itau
    G_ao3(1:nBas2,1:nBas2) = G_ao_itau(2*itau-1,      1:nBas2,      1:nBas2)  ! he   itau
    Ainter=czero;Binter=czero;Cinter=czero;
    call Sigma_c_GGF2B_he_prime(nBas2,Ainter,Binter,Cinter,G_ao1,G_ao2,G_ao3,db_ERI_AO,Sigma_c_minus(1:nBas2,1:nBas2))
    ! Sigma_he_2primeprime
    G_ao1(1:nBas2,1:nBas2) = G_ao_itau(2*itau  ,      1:nBas2,      1:nBas2)  ! he  -itau
    G_ao2(1:nBas2,1:nBas2) = G_ao_itau(2*itau  ,      1:nBas2,nBas2+1:nBas4)  ! hh  -itau
    G_ao3(1:nBas2,1:nBas2) = G_ao_itau(2*itau-1,nBas2+1:nBas4,      1:nBas2)  ! ee   itau
    Ainter=czero;Binter=czero;Cinter=czero;
    call Sigma_c_GGF2B_he_prime2(nBas2,Ainter,Binter,Cinter,G_ao1,G_ao2,G_ao3,db_ERI_AO,Sigma_c_minus(1:nBas2,1:nBas2))
    ! Sigma_eh_2prime
    G_ao1(1:nBas2,1:nBas2) = G_ao_itau(2*itau  ,nBas2+1:nBas4,nBas2+1:nBas4)  ! eh  -itau
    G_ao2(1:nBas2,1:nBas2) = G_ao_itau(2*itau  ,nBas2+1:nBas4,nBas2+1:nBas4)  ! eh  -itau
    G_ao3(1:nBas2,1:nBas2) = G_ao_itau(2*itau-1,nBas2+1:nBas4,nBas2+1:nBas4)  ! eh   itau
    Ainter=czero;Binter=czero;Cinter=czero;
    call Sigma_c_GGF2B_eh_prime(nBas2,Ainter,Binter,Cinter,G_ao1,G_ao2,G_ao3,db_ERI_AO,Sigma_c_minus(nBas2+1:nBas4,nBas2+1:nBas4))
    ! Sigma_eh_2primeprime
    G_ao1(1:nBas2,1:nBas2) = G_ao_itau(2*itau  ,nBas2+1:nBas4,nBas2+1:nBas4)  ! eh  -itau
    G_ao2(1:nBas2,1:nBas2) = G_ao_itau(2*itau  ,      1:nBas2,nBas2+1:nBas4)  ! hh  -itau
    G_ao3(1:nBas2,1:nBas2) = G_ao_itau(2*itau-1,nBas2+1:nBas4,      1:nBas2)  ! ee   itau
    Ainter=czero;Binter=czero;Cinter=czero;
    call Sigma_c_GGF2B_eh_prime2(nBas2,Ainter,Binter,Cinter,G_ao1,G_ao2,G_ao3,db_ERI_AO,Sigma_c_minus(nBas2+1:nBas4,nBas2+1:nBas4))
    ! Sigma_hh_2prime
    G_ao1(1:nBas2,1:nBas2) = G_ao_itau(2*itau  ,      1:nBas2,nBas2+1:nBas4)  ! hh  -itau
    G_ao2(1:nBas2,1:nBas2) = G_ao_itau(2*itau  ,      1:nBas2,      1:nBas2)  ! he  -itau
    G_ao3(1:nBas2,1:nBas2) = G_ao_itau(2*itau-1,      1:nBas2,      1:nBas2)  ! he   itau
    Ainter=czero;Binter=czero;Cinter=czero;
    call Sigma_c_GGF2B_hh_prime(nBas2,Ainter,Binter,Cinter,G_ao1,G_ao2,G_ao3,db_ERI_AO,Sigma_c_minus(1:nBas2,nBas2+1:nBas4))
    ! Sigma_hh_2primeprime
    G_ao1(1:nBas2,1:nBas2) = G_ao_itau(2*itau  ,      1:nBas2,nBas2+1:nBas4)  ! hh  -itau
    G_ao2(1:nBas2,1:nBas2) = G_ao_itau(2*itau  ,      1:nBas2,nBas2+1:nBas4)  ! hh  -itau
    G_ao3(1:nBas2,1:nBas2) = G_ao_itau(2*itau-1,nBas2+1:nBas4,      1:nBas2)  ! ee   itau
    Ainter=czero;Binter=czero;Cinter=czero;
    call Sigma_c_GGF2B_hh_prime2(nBas2,Ainter,Binter,Cinter,G_ao1,G_ao2,G_ao3,db_ERI_AO,Sigma_c_minus(1:nBas2,nBas2+1:nBas4))
    ! Sigma_ee_2prime
    G_ao1(1:nBas2,1:nBas2) = G_ao_itau(2*itau  ,nBas2+1:nBas4,      1:nBas2)  ! ee  -itau
    G_ao2(1:nBas2,1:nBas2) = G_ao_itau(2*itau  ,      1:nBas2,      1:nBas2)  ! he  -itau
    G_ao3(1:nBas2,1:nBas2) = G_ao_itau(2*itau-1,      1:nBas2,      1:nBas2)  ! he   itau
    Ainter=czero;Binter=czero;Cinter=czero;
    call Sigma_c_GGF2B_ee_prime(nBas2,Ainter,Binter,Cinter,G_ao1,G_ao2,G_ao3,db_ERI_AO,Sigma_c_minus(nBas2+1:nBas4,1:nBas2))
    ! Sigma_ee_2primeprime
    G_ao1(1:nBas2,1:nBas2) = G_ao_itau(2*itau  ,nBas2+1:nBas4,      1:nBas2)  ! ee  -itau
    G_ao2(1:nBas2,1:nBas2) = G_ao_itau(2*itau  ,      1:nBas2,nBas2+1:nBas4)  ! hh  -itau
    G_ao3(1:nBas2,1:nBas2) = G_ao_itau(2*itau-1,nBas2+1:nBas4,      1:nBas2)  ! ee   itau
    Ainter=czero;Binter=czero;Cinter=czero;
    call Sigma_c_GGF2B_ee_prime2(nBas2,Ainter,Binter,Cinter,G_ao1,G_ao2,G_ao3,db_ERI_AO,Sigma_c_minus(nBas2+1:nBas4,1:nBas2))

   inquire(file='Bog_Sigma_c_M8', exist=file_exists)
   if(file_exists) then
    call Sigma_c_GGF2B_brut(nBas2,nBas4,G_ao_itau(2*itau-1,:,:),G_ao_itau(2*itau,:,:),db_ERI_AO,Sigma_c_plus,Sigma_c_minus) 
   endif

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
    do abas=1,nBas2
     EcGM_itau=EcGM_itau+tweight(itau)*Mat_gorkov_tmp(abas,abas)
    enddo
    ! tau < 0
    Mat_gorkov_tmp(:,:) =G_ao_itau(2*itau,:,:)   ! G_ao_itau (-)
    Mat_gorkov_tmp=matmul(Sigma_c_plus,Mat_gorkov_tmp)
    do abas=1,nBas2
     EcGM_itau=EcGM_itau+tweight(itau)*Mat_gorkov_tmp(abas,abas)
    enddo
  enddo
  EcGM=-0.5d0*real(EcGM_itau) ! Including sum over spin-channels  EcGM = - 1/2 \sum_spin \int Tr[ Sigma_c_spin(-it) G_spin(it) ]^he dt

  ! Converge with respect to the H_HFB operator (using only good Gen_R_ao matrices -> Tr[R_ao_block S_ao]=Nelectrons )
  if(.not.no_h_hfb) then ! Skiiping the opt w.r.t. the H_HFB operator to do later the linearized approximation on Go -> [ lin-G = Go + Go Sigma Go ]
   iter_hfb=0
   n_diisR=0
   rcondR=0d0
   err_diisR=0d0
   Gen_R_ao_old_diis=0d0
   do
    ! Build H_HFB and energy components
    iter_hfb=iter_hfb+1
    Ehfbl=0d0
    Ecore=0d0; Ehx=0d0; Epair=0d0;
    Gen_H_ao_hfb=0d0
    Gen_H_ao_hfb(1:nBas2      ,1:nBas2      ) =  Gen_Hc(1:nBas2,1:nBas2) 
    Gen_H_ao_hfb(nBas2+1:nBas4,nBas2+1:nBas4) = -Gen_Hc(1:nBas2,1:nBas2) 
    do abas=1,nBas2
     pbas=abas+nBas2
     do bbas=1,nBas2
      qbas=bbas+nBas2
      do cbas=1,nBas2
       rbas=cbas+nBas2
       do dbas=1,nBas2
        sbas=dbas+nBas2
        Gen_H_ao_hfb(abas,bbas)=Gen_H_ao_hfb(abas,bbas)+Gen_R_ao(cbas,dbas)*db_ERI_AO(abas,cbas,bbas,dbas) 
        Gen_H_ao_hfb(pbas,qbas)=Gen_H_ao_hfb(pbas,qbas)-Gen_R_ao(cbas,dbas)*db_ERI_AO(bbas,dbas,abas,cbas) 
        Gen_H_ao_hfb(abas,qbas)=Gen_H_ao_hfb(abas,qbas)+0.5d0*sigma*Gen_R_ao(cbas,sbas)*db_ERI_AO(abas,bbas,cbas,dbas)
        Gen_H_ao_hfb(pbas,bbas)=Gen_H_ao_hfb(pbas,bbas)+0.5d0*sigma*Gen_R_ao(rbas,dbas)*db_ERI_AO(abas,bbas,cbas,dbas)
       enddo
      enddo
     enddo
    enddo
     ! Core
    Mat_gorkov_tmp=czero
    Mat_gorkov_tmp(1:nBas2      ,1:nBas2      ) =  Gen_Hc(1:nBas2,1:nBas2) 
    Mat_gorkov_tmp(nBas2+1:nBas4,nBas2+1:nBas4) = -Gen_Hc(1:nBas2,1:nBas2) 
    Mat_gorkov_tmp=matmul(Mat_gorkov_tmp,Gen_R_ao)
    do abas=1,nBas2
     Ecore=Ecore+real(Mat_gorkov_tmp(abas,abas))
    enddo
     ! Hx
    Mat_gorkov_tmp=Gen_H_ao_hfb
    Mat_gorkov_tmp(1:nBas2      ,nBas2+1:nBas4)=czero
    Mat_gorkov_tmp(nBas2+1:nBas4,1:nBas2      )=czero
    Mat_gorkov_tmp(1:nBas2      ,1:nBas2      )=Mat_gorkov_tmp(1:nBas2      ,1:nBas2      ) - Gen_Hc(1:nBas2,1:nBas2) 
    Mat_gorkov_tmp(nBas2+1:nBas4,nBas2+1:nBas4)=Mat_gorkov_tmp(nBas2+1:nBas4,nBas2+1:nBas4) + Gen_Hc(1:nBas2,1:nBas2) 
    Mat_gorkov_tmp=matmul(Mat_gorkov_tmp,Gen_R_ao)
    do abas=1,nBas2
     Ehx=Ehx+real(Mat_gorkov_tmp(abas,abas))
    enddo
    Ehx=0.5d0*Ehx
     ! Pair
    Mat_gorkov_tmp=Gen_H_ao_hfb
    Mat_gorkov_tmp(1:nBas2      ,1:nBas2      )=czero
    Mat_gorkov_tmp(nBas2+1:nBas4,nBas2+1:nBas4)=czero
    Mat_gorkov_tmp=matmul(Mat_gorkov_tmp,Gen_R_ao)
    do abas=1,nBas2
     Epair=Epair+real(Mat_gorkov_tmp(abas,abas))
    enddo
    Epair=0.5d0*Epair
     ! All
    Ehfbl=Ecore+Ehx+Epair
    ! Build G(i w) and R
    Gen_R_ao_old=Gen_R_ao
    call get_1rdm_scGXB(nBas2,nBas4,nfreqs,chem_pot,Gen_S,Gen_H_ao_hfb,Sigma_c_w_ao,wcoord,wweight, &
                        Mat_gorkov_tmp,G_ao_iw_hfb,DeltaG_ao_iw,Gen_R_ao,Gen_R_ao_hfb,trace_1_rdm) 
    if(abs(trace_1_rdm-nElectrons)**2d0>thrs_N .and. chem_pot_scG) &
     call fix_chem_pot_scGXB_bisec(iter_hfb,nBas2,nBas4,nfreqs,nElectrons,thrs_N,thrs_Ngrad,chem_pot,Gen_S,Gen_H_ao_hfb,Sigma_c_w_ao,   &
                                   wcoord,wweight,Mat_gorkov_tmp,G_ao_iw_hfb,DeltaG_ao_iw,Gen_R_ao,Gen_R_ao_hfb,trace_1_rdm,            &
                                   chem_pot_saved,verbose_scGF2B)
    if(abs(trace_1_rdm-nElectrons)**2d0>thrs_N .and. .not.chem_pot_scG) &
     Gen_R_ao=nElectrons*Gen_R_ao/trace_1_rdm

    ! Check convergence of Gen_R_ao for fixed Sigma_c(i w)
    diff_Rao=0d0
    do abas=1,nBas4
     do bbas=1,nBas4
      diff_Rao=diff_Rao+abs(Gen_R_ao(abas,bbas)-Gen_R_ao_old(abas,bbas))
     enddo
    enddo

    if(diff_Rao<=thrs_Rao) exit

    if(iter_hfb==maxSCF) exit

    ! Do mixing with previous Gen_R_ao to facilitate convergence
    if(maxDIIS>0) then
     n_diisR=min(n_diisR+1,maxDIIS)
     err_currentR=0d0
     idiis_indexR=1
     do abas=1,nBas4
      do bbas=1,nBas4
       err_currentR(idiis_indexR)=Gen_R_ao(abas,bbas)-Gen_R_ao_old(abas,bbas)
       if(abs(err_currentR(idiis_indexR))<1e-12) err_currentR(idiis_indexR)=0d0
       Gen_R_ao_extrap(idiis_indexR)=Gen_R_ao(abas,bbas)
       idiis_indexR=idiis_indexR+1
      enddo
     enddo
     call DIIS_extrapolation(rcondR,nBas4Sq,nBas4Sq,n_diisR,err_diisR,Gen_R_ao_old_diis,err_currentR,Gen_R_ao_extrap)
     idiis_indexR=1
     do abas=1,nBas4
      do bbas=1,nBas4
       Gen_R_ao(abas,bbas)=Gen_R_ao_extrap(idiis_indexR)
       idiis_indexR=idiis_indexR+1
      enddo
     enddo
    else
     Gen_R_ao(:,:)=alpha_mixing*Gen_R_ao(:,:)+(1d0-alpha_mixing)*Gen_R_ao_old(:,:)
    endif

   enddo
  endif

  ! Check convergence of Gen_R_ao after a scGF2B iteration
  diff_Rao=0d0
  do abas=1,nBas4
   do bbas=1,nBas4
    diff_Rao=diff_Rao+abs(Gen_R_ao(abas,bbas)-Gen_R_ao_iter(abas,bbas))
   enddo 
  enddo
  Gen_R_ao_iter=Gen_R_ao

  ! Print iter info
  U_mo=-matmul(matmul(Gen_cHFBinv,Gen_R_ao(1:nBas2,1:nBas2)),transpose(Gen_cHFBinv)) ! Minus to order occ numbers
  call diagonalize_matrix(nOrb2,U_mo,Occ)
  Occ=-Occ
  trace_1_rdm=sum(Occ)
  Gen_cNO=matmul(Gen_cHFB,U_mo)
  write(*,*)
  write(*,'(a,f15.8,a,i5,a,i5)') ' Trace scGGF2B ',trace_1_rdm,' after ',iter_hfb,' HFB iterations at global iter ',iter
  write(*,'(a,f15.8)')        ' Change of R   ',diff_Rao
  write(*,'(a,f15.8)')        ' Chem. Pot.    ',chem_pot
  write(*,'(a,f15.8)')        ' Enuc          ',ENuc
  write(*,'(a,f15.8)')        ' Ehcore        ',Ecore
  write(*,'(a,f15.8)')        ' Hx            ',Ehx
  write(*,'(a,f15.8)')        ' Epairing      ',Epair
  write(*,'(a,f15.8)')        ' Ehfbl         ',Ehfbl
  write(*,'(a,f15.8)')        ' EcGM          ',EcGM
  write(*,'(a,f15.8)')        ' Eelec         ',Ehfbl+EcGM
  write(*,'(a,f15.8)')        ' Etot          ',Ehfbl+EcGM+ENuc
  write(*,'(a,f15.8)')        ' EcPT2         ',EcGM/2d0
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
    do abas=1,nBas4
     do bbas=1,nBas4
      err_current(idiis_index)=G_ao_itau(itau,abas,bbas)-G_ao_itau_old(itau,abas,bbas)
      if(abs(err_current(idiis_index))<1e-12) err_current(idiis_index)=czero
      G_itau_extrap(idiis_index)=G_ao_itau(itau,abas,bbas)
      idiis_index=idiis_index+1
     enddo
    enddo
   enddo
   call complex_DIIS_extrapolation(rcond,nBas4Sqntimes2,nBas4Sqntimes2,n_diis,err_diis,G_itau_old_diis,err_current,G_itau_extrap)
   idiis_index=1
   do itau=1,ntimes_twice
    do abas=1,nBas4
     do bbas=1,nBas4
      G_ao_itau(itau,abas,bbas)=G_itau_extrap(idiis_index)
      idiis_index=idiis_index+1
     enddo
    enddo
   enddo
  else
   G_ao_itau(:,:,:)=alpha_mixing*G_ao_itau(:,:,:)+(1d0-alpha_mixing)*G_ao_itau_old(:,:,:)
  endif
  G_ao_itau_old(:,:,:)=G_ao_itau(:,:,:)

 enddo
 N_anom = trace_matrix(nBas2/2,matmul(transpose(Gen_R_ao(1:nBas2/2,nBas2/2+nBas2+1:nBas4)), &
          Gen_R_ao(1:nBas2/2,nBas2/2+nBas2+1:nBas4)))
 N_anom = N_anom + trace_matrix(nBas2/2,matmul(transpose(Gen_R_ao(nBas2/2+1:nBas2,nBas2+1:nBas2+nBas2/2)), &
          Gen_R_ao(nBas2/2+1:nBas2,nBas2+1:nBas2+nBas2/2)))
 write(*,*)
 write(*,'(A50)') '----------------------------------------'
 write(*,'(A50)') '     scGGF2B calculation completed      '
 write(*,'(A50)') '----------------------------------------'
 write(*,*)
 write(*,'(a,f15.8,a,i5,a)') ' Trace scGGF2B   ',trace_1_rdm,' after ',iter,' global iterations '
 write(*,'(a,f15.8)')        ' Change of R     ',diff_Rao
 write(*,'(a,f15.8)')        ' Chem. Pot.      ',chem_pot
 write(*,'(a,f15.8)')        ' N anomalus      ',N_anom
 write(*,'(a,f15.8)')        ' Enuc            ',ENuc
 write(*,'(a,f15.8)')        ' Ehcore          ',Ecore
 write(*,'(a,f15.8)')        ' Hartree-Exchange',Ehx
 write(*,'(a,f15.8)')        ' Epairing        ',Epair
 write(*,'(a,f15.8)')        ' Ehfbl           ',Ehfbl
 write(*,'(a,f15.8)')        ' EcGM            ',EcGM
 write(*,'(a,f15.8)')        ' Eelec           ',Ehfbl+EcGM
 write(*,'(a,f15.8)')        ' scGGF2B Energy  ',Ehfbl+EcGM+ENuc
 write(*,*)                                    
 write(*,'(a,f15.8)')        ' EcPT2           ',EcGM/2d0
 write(*,*)
 write(*,*) ' Final occupation numbers'
 do abas=1,nOrb2
  write(*,'(I7,F15.8)') abas,Occ(abas)
 enddo
 if(verbose/=0) then
  write(*,*) ' Natural orbitals (columns)'
  do abas=1,nBas2
   write(*,'(*(f15.8))') Gen_cNO(abas,:)
  enddo
 endif
 write(*,*)

 ! Write restart files
 call write_scGXB_restart(nBas4,ntimes,ntimes_twice,nfreqs,chem_pot,Gen_R_ao,Gen_R_ao_hfb,G_ao_itau,G_ao_itau_hfb, &
                         G_ao_iw_hfb,DeltaG_ao_iw)

 inquire(file='Print_Rao', exist=file_exists)
 if(file_exists) then
  write(*,*) 'R_scGGF2B_ao'
  do abas=1,nBas4
   write(*,'(*(f10.5))') Gen_R_ao(abas,:)
  enddo
 endif
 write(*,*)

 ! Deallocate arrays
 deallocate(Occ)
 deallocate(Gen_cHFBinv)
 deallocate(Gen_cNO)
 deallocate(U_mo)
 deallocate(Gen_R_ao)
 deallocate(Gen_R_ao_iter)
 deallocate(Gen_R_ao_hfb)
 deallocate(Gen_R_ao_old)
 deallocate(Gen_H_ao_hfb)
 deallocate(G_ao_tmp)
 deallocate(Mat_gorkov_tmp)
 deallocate(Mat_gorkov_tmp2)
 deallocate(DeltaG_ao_iw)
 deallocate(G_ao_iw_hfb)
 deallocate(G_ao_itau_hfb)
 deallocate(G_ao_itau)
 deallocate(G_ao_itau_old)
 deallocate(Sigma_c_w_ao)
 deallocate(Sigma_c_plus)
 deallocate(Sigma_c_minus)
 deallocate(Sigma_c_c)
 deallocate(Sigma_c_s)
 deallocate(Mat1,Mat2)
 deallocate(Mat3,Mat4)
 deallocate(Ainter)
 deallocate(Binter)
 deallocate(Cinter)
 deallocate(G_ao1)
 deallocate(G_ao2)
 deallocate(G_ao3)
 deallocate(tweight,tcoord)
 deallocate(sint2w_weight)
 deallocate(cost2w_weight)
 deallocate(cosw2t_weight)
 deallocate(sinw2t_weight)
 deallocate(err_currentR)
 deallocate(Gen_R_ao_extrap)
 deallocate(err_diisR)
 deallocate(Gen_R_ao_old_diis)
 deallocate(G_itau_extrap)
 deallocate(err_current)
 deallocate(err_diis)
 deallocate(G_itau_old_diis)

 call wall_time(end_scGF2Bitauiw)
 
 t_scGF2Bitauiw = end_scGF2Bitauiw - start_scGF2Bitauiw
 write(*,'(A65,1X,F9.3,A8)') 'Total wall time for scGF2B = ',t_scGF2Bitauiw,' seconds'
 write(*,*)

end subroutine 
