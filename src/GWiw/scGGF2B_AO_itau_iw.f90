subroutine scGGF2B_AO_itau_iw(nBas_twice,nBas_4times,nOrb_twice,nOrb_4times,maxSCF,thresh_in,maxDIIS,restart_scGF2B,verbose_scGF2B, &
                              chem_pot_scG,no_h_hfb,ENuc,Gen_Hc,Gen_S,Gen_R_in,Gen_cHFB,Gen_eQP_state,chem_pot,sigma,               &
                              nfreqs,wcoord,wweight,Gen_U_QP,db_ERI_AO)

! Generalized scGF2B

  implicit none
  include 'parameters.h'

! Input variables
 
  logical,intent(in)            :: no_h_hfb
  logical,intent(in)            :: restart_scGF2B
  logical,intent(in)            :: verbose_scGF2B
  logical,intent(in)            :: chem_pot_scG

  integer,intent(in)            :: nBas_twice
  integer,intent(in)            :: nBas_4times
  integer,intent(in)            :: nOrb_twice
  integer,intent(in)            :: nOrb_4times
  integer,intent(in)            :: maxSCF
  integer,intent(in)            :: maxDIIS

  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: thresh_in
  double precision,intent(in)   :: sigma
  double precision,intent(in)   :: Gen_Hc(nBas_twice,nBas_twice)
  double precision,intent(in)   :: Gen_R_in(nBas_4times,nBas_4times)
  double precision,intent(in)   :: Gen_S(nBas_twice,nBas_twice)
  double precision,intent(in)   :: db_ERI_AO(nBas_twice,nBas_twice,nBas_twice,nBas_twice)
  double precision,intent(in)   :: Gen_U_QP(nOrb_4times,nOrb_4times)

! Local variables

  logical                       :: file_exists
  logical                       :: read_HFB_chkp

  integer                       :: ifreq,itau
  integer                       :: ibas,jbas,kbas,lbas,obas,qbas
  integer                       :: iorb,jorb,korb,lorb
  integer                       :: iter,iter_hfb
  integer                       :: verbose
  integer                       :: ntimes
  integer                       :: nBas_4timesSq
  integer                       :: ntimes_twice
  integer                       :: nfreqs_int
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
  double precision              :: EcGM,Ehfbl,Ecore,Eh,Ex,Epair
  double precision              :: trace_1_rdm
  double precision,external     :: trace_matrix
  double precision              :: start_scGF2Bitauiw,end_scGF2Bitauiw,t_scGF2Bitauiw
  double precision,allocatable  :: Occ(:)
  double precision,allocatable  :: Gen_cHFBinv(:,:)
  double precision,allocatable  :: Gen_cNO(:,:)
  double precision,allocatable  :: U_mo(:,:)
  double precision,allocatable  :: Mat1(:,:)
  double precision,allocatable  :: Mat2(:,:)
  double precision,allocatable  :: R_ao(:,:)
  double precision,allocatable  :: R_ao_iter(:,:)
  double precision,allocatable  :: R_ao_hfb(:,:)
  double precision,allocatable  :: R_ao_old(:,:)
  double precision,allocatable  :: H_ao_hfb(:,:)
  double precision,allocatable  :: R_mo(:,:)
  double precision,allocatable  :: Gen_cHFB_gorkov(:,:)
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
  complex*16,allocatable        :: G_ao_itau(:,:,:)
  complex*16,allocatable        :: G_ao_itau_old(:,:,:)
  complex*16,allocatable        :: G_ao_itau_hfb(:,:,:)
  complex*16,allocatable        :: err_diis(:,:)
  complex*16,allocatable        :: G_itau_old_diis(:,:)
  complex*16,allocatable        :: G_ao_iw_hfb(:,:,:)
  complex*16,allocatable        :: DeltaG_ao_iw(:,:,:)

! Output variables
  integer,intent(inout)         :: nfreqs
  double precision,intent(inout):: chem_pot
  double precision,intent(inout):: Gen_eQP_state(nOrb_twice)
  double precision,intent(inout):: wcoord(nfreqs)
  double precision,intent(inout):: wweight(nfreqs)
  double precision,intent(inout):: Gen_cHFB(nBas_twice,nOrb_twice)
  
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

 write(*,'(a,F15.8)') '  emin ',abs(Gen_eQP_state(nOrb_twice))
 write(*,'(a,F15.8)') '  emax ',abs(Gen_eQP_state(nOrb_twice)+Gen_eQP_state(1))

 ! Initialize variables
 ntimes=nfreqs
 ntimes_twice=ntimes*2
 nBas_4timesSq=nBas_4times*nBas_4times
 verbose=0
 if(verbose_scGF2B) verbose=1     
 eta=0d0
 thrs_N=1d-10
 thrs_Ngrad=1d-6
 thrs_Rao=thresh_in
 nElectrons=0d0
 do ibas=1,nBas_twice
  do jbas=1,nBas_twice
   nElectrons=nElectrons+Gen_R_in(ibas,jbas)*Gen_S(ibas,jbas)
  enddo
 enddo
 nElectrons=nint(nElectrons) ! Here we prefer to use 1 spin-channel
 chem_pot_saved=chem_pot
 alpha_mixing=0.6d0
 rcond=0d0
 n_diis=0
 nBas4Sqntimes2=nBas_4timesSq*ntimes_twice

 ! Allocate arrays
 allocate(Occ(nOrb_twice))
 allocate(Gen_cNO(nBas_twice,nOrb_twice))
 allocate(U_mo(nOrb_twice,nOrb_twice))
 allocate(Gen_cHFBinv(nOrb_twice,nBas_twice))
 allocate(Gen_cHFB_gorkov(nBas_4times,nOrb_4times))
 allocate(R_ao(nBas_4times,nBas_4times))
 allocate(R_ao_iter(nBas_4times,nBas_4times))
 allocate(R_ao_hfb(nBas_4times,nBas_4times))
 allocate(R_ao_old(nBas_4times,nBas_4times))
 allocate(R_mo(nOrb_4times,nOrb_4times))
 allocate(H_ao_hfb(nBas_4times,nBas_4times))
 allocate(G_ao_tmp(nBas_twice,nBas_twice))
 allocate(Mat_gorkov_tmp(nBas_4times,nBas_4times))
 allocate(Mat_gorkov_tmp2(nBas_4times,nBas_4times))
 allocate(G_ao_iw_hfb(nfreqs,nBas_4times,nBas_4times))
 allocate(DeltaG_ao_iw(nfreqs,nBas_4times,nBas_4times))
 allocate(G_ao_itau_hfb(ntimes_twice,nBas_4times,nBas_4times))
 allocate(G_ao_itau(ntimes_twice,nBas_4times,nBas_4times))
 allocate(G_ao_itau_old(ntimes_twice,nBas_4times,nBas_4times))
 allocate(Sigma_c_w_ao(nfreqs,nBas_4times,nBas_4times))
 allocate(Sigma_c_plus(nBas_4times,nBas_4times)) 
 allocate(Sigma_c_minus(nBas_4times,nBas_4times)) 
 allocate(Sigma_c_c(nBas_4times,nBas_4times))      
 allocate(Sigma_c_s(nBas_4times,nBas_4times))
 allocate(tweight(ntimes),tcoord(ntimes))
 allocate(sint2w_weight(nfreqs,ntimes))
 allocate(cost2w_weight(nfreqs,ntimes))
 allocate(cosw2t_weight(ntimes,nfreqs))
 allocate(sinw2t_weight(ntimes,nfreqs))
 allocate(Mat1(nOrb_twice,nOrb_twice),Mat2(nOrb_twice,nOrb_twice))
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
  allocate(err_currentR(nBas_4timesSq))
  allocate(R_ao_extrap(nBas_4timesSq))
  allocate(err_diisR(nBas_4timesSq,maxDIIS))
  allocate(R_ao_old_diis(nBas_4timesSq,maxDIIS))
  allocate(G_itau_extrap(nBas4Sqntimes2))
  allocate(err_current(nBas4Sqntimes2))
  allocate(err_diis(nBas4Sqntimes2,maxDIIS))
  allocate(G_itau_old_diis(nBas4Sqntimes2,maxDIIS))
 endif
 err_diis=czero
 G_itau_old_diis=czero

 ! Initialize arrays
 ! If required, read restart files

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
    do ibas=1,nBas_twice
     EcGM_itau=EcGM_itau+tweight(itau)*Mat_gorkov_tmp(ibas,ibas)
    enddo
    ! tau < 0
    Mat_gorkov_tmp(:,:) =G_ao_itau(2*itau,:,:)   ! G_ao_itau (-)
    Mat_gorkov_tmp=matmul(Sigma_c_plus,Mat_gorkov_tmp)
    do ibas=1,nBas_twice
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
    H_ao_hfb(1:nBas_twice,1:nBas_twice)=Gen_Hc(1:nBas_twice,1:nBas_twice)
    Ehfbl=0d0
    Ecore=0d0; Eh=0d0; Ex=0d0; Epair=0d0;
    do ibas=1,nBas_twice
     do jbas=1,nBas_twice
      obas=nBas_twice+1+(jbas-1)
      Ehfbl=Ehfbl+2d0*R_ao(ibas,jbas)*Gen_Hc(ibas,jbas)
      Ecore=Ecore+2d0*R_ao(ibas,jbas)*Gen_Hc(ibas,jbas)
      do kbas=1,nBas_twice
       do lbas=1,nBas_twice
        qbas=nBas_twice+1+(lbas-1)
        H_ao_hfb(ibas,jbas)=H_ao_hfb(ibas,jbas)
        Ehfbl=Ehfbl+2d0*R_ao(kbas,lbas)*R_ao(ibas,jbas)
        Eh=Eh
        Ex=Ex 
        Epair=Epair+sigma 
       enddo
      enddo
     enddo
    enddo
    ! Build G(i w) and R
    R_ao_old=R_ao
    call get_1rdm_scGXB(nBas_twice,nBas_4times,nfreqs,chem_pot,Gen_S,H_ao_hfb,Sigma_c_w_ao,wcoord,wweight, &
                        Mat_gorkov_tmp,G_ao_iw_hfb,DeltaG_ao_iw,R_ao,R_ao_hfb,trace_1_rdm) 
    if(abs(trace_1_rdm-nElectrons)**2d0>thrs_N .and. chem_pot_scG) &
     call fix_chem_pot_scGXB_bisec(iter_hfb,nBas_twice,nBas_4times,nfreqs,nElectrons,thrs_N,thrs_Ngrad,chem_pot,Gen_S,H_ao_hfb,Sigma_c_w_ao,   &
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
     call DIIS_extrapolation(rcondR,nBas_4timesSq,nBas_4timesSq,n_diisR,err_diisR,R_ao_old_diis,err_currentR,R_ao_extrap)
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
  U_mo=-2d0*matmul(matmul(Gen_cHFBinv,R_ao(1:nBas_twice,1:nBas_twice)),transpose(Gen_cHFBinv)) ! Minus to order occ numbers
  call diagonalize_matrix(nOrb_twice,U_mo,Occ)
  Occ=-Occ
  trace_1_rdm=sum(Occ)
  Gen_cNO=matmul(Gen_cHFB,U_mo)
  write(*,*)
  write(*,'(a,f15.8,a,i5,a,i5)') ' Trace scGGF2B ',trace_1_rdm,' after ',iter_hfb,' HFB iterations at global iter ',iter
  write(*,'(a,f15.8)')        ' Change of R   ',diff_Rao
  write(*,'(a,f15.8)')        ' Chem. Pot.    ',chem_pot
  write(*,'(a,f15.8)')        ' Enuc          ',ENuc
  write(*,'(a,f15.8)')        ' Ehcore        ',Ecore
  write(*,'(a,f15.8)')        ' Hartree       ',Eh
  write(*,'(a,f15.8)')        ' Exchange      ',Ex
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
    do ibas=1,nBas_twice
     do jbas=1,nBas_twice
      err_current(idiis_index)=G_ao_itau(itau,ibas,jbas)-G_ao_itau_old(itau,ibas,jbas)
      if(abs(err_current(idiis_index))<1e-12) err_current(idiis_index)=czero
      G_itau_extrap(idiis_index)=G_ao_itau(itau,ibas,jbas)
      idiis_index=idiis_index+1
     enddo
    enddo
   enddo
   call complex_DIIS_extrapolation(rcond,nBas4Sqntimes2,nBas4Sqntimes2,n_diis,err_diis,G_itau_old_diis,err_current,G_itau_extrap)
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
! N_anom = trace_matrix(nBas,matmul(transpose(R_ao(1:nBas,nBas+1:nBas_twice)), &
!          R_ao(1:nBas,nBas+1:nBas_twice)))
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
 write(*,'(a,f15.8)')        ' Hartree         ',Eh
 write(*,'(a,f15.8)')        ' Exchange        ',Ex
 write(*,'(a,f15.8)')        ' Epairing        ',Epair
 write(*,'(a,f15.8)')        ' Ehfbl           ',Ehfbl
 write(*,'(a,f15.8)')        ' EcGM            ',EcGM
 write(*,'(a,f15.8)')        ' Eelec           ',Ehfbl+EcGM
 write(*,'(a,f15.8)')        ' scGGF2B Energy  ',Ehfbl+EcGM+ENuc
 write(*,*)                                    
 write(*,'(a,f15.8)')        ' EcPT2           ',EcGM/2d0
 write(*,*)
 write(*,*) ' Final occupation numbers'
 do ibas=1,nOrb_twice
  write(*,'(I7,F15.8)') ibas,Occ(ibas)
 enddo
 if(verbose/=0) then
  write(*,*) ' Natural orbitals (columns)'
  do ibas=1,nBas_twice
   write(*,'(*(f15.8))') Gen_cNO(ibas,:)
  enddo
 endif
 write(*,*)

 ! Write restart files
 call write_scGXB_restart(nBas_twice,ntimes,ntimes_twice,nfreqs,chem_pot,R_ao,R_ao_hfb,G_ao_itau,G_ao_itau_hfb, &
                         G_ao_iw_hfb,DeltaG_ao_iw)

 inquire(file='Print_Rao', exist=file_exists)
 if(file_exists) then
  write(*,*) 'R_scGGF2B_ao'
  do ibas=1,nBas_twice
   write(*,'(*(f10.5))') R_ao(ibas,:)
  enddo
 endif

 ! Deallocate arrays
 deallocate(Occ)
 deallocate(Gen_cHFBinv)
 deallocate(Gen_cHFB_gorkov)
 deallocate(Gen_cNO)
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
