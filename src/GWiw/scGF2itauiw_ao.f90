subroutine scGF2itauiw_ao(nBas,nOrb,nO,maxSCF,maxDIIS,dolinGF2,restart_scGF2,verbose_scGF2,no_fock,ENuc,Hc,S,P_in,cHF,eHF, &
                          nfreqs,wcoord,wweight,vMAT,ERI_AO)

! Restricted scGF2

  implicit none
  include 'parameters.h'

! Input variables
 
  logical,intent(in)            :: dolinGF2
  logical,intent(in)            :: no_fock
  logical,intent(in)            :: restart_scGF2
  logical,intent(in)            :: verbose_scGF2

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nO
  integer,intent(in)            :: maxSCF
  integer,intent(in)            :: maxDIIS

  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: Hc(nBas,nBas)
  double precision,intent(in)   :: P_in(nBas,nBas)
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: vMAT(nBas*nBas,nBas*nBas)
  double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas)

! Local variables
 
  logical                       :: file_exists
  logical                       :: read_SD_chkp

  integer                       :: n_diis
  integer                       :: n_diisP
  integer                       :: verbose
  integer                       :: nneg
  integer                       :: ntimes
  integer                       :: nBasSqntimes2
  integer                       :: ntimes_twice
  integer                       :: idiis_index
  integer                       :: idiis_indexP
  integer                       :: itau,ifreq
  integer                       :: ibas,jbas,kbas,lbas,nBas2
  integer                       :: mbas,sbas,pbas,qbas
  integer                       :: iter,iter_fock
  integer                       :: imax_error_sigma
  integer                       :: imax_error_gw2gt

  double precision              :: start_scGF2itauiw     ,end_scGF2itauiw       ,t_scGF2itauiw

  double precision              :: rcond
  double precision              :: rcondP
  double precision              :: alpha_mixing
  double precision              :: Ehfl,EcGM
!  double precision              :: EcGMw
  double precision              :: trace1,trace2
  double precision              :: eta,diff_Pao
  double precision              :: nElectrons
  double precision              :: trace_1_rdm
  double precision              :: err_EcGM
  double precision              :: thrs_N,thrs_Ngrad,thrs_Pao
  double precision              :: chem_pot,chem_pot_saved,chem_pot_align
  double precision              :: error_sigma
  double precision              :: max_error_sigma
  double precision              :: error_gw2gt
  double precision              :: max_error_gw2gt
  double precision              :: sum_error_gw2gt
  double precision              :: sd_dif
  double precision,allocatable  :: tweight(:),tcoord(:)
  double precision,allocatable  :: sint2w_weight(:,:)
  double precision,allocatable  :: cost2w_weight(:,:)
  double precision,allocatable  :: cosw2t_weight(:,:)
  double precision,allocatable  :: sinw2t_weight(:,:)
  double precision,allocatable  :: eSD(:)
  double precision,allocatable  :: eSD_old(:)
  double precision,allocatable  :: Occ(:)
  double precision,allocatable  :: Wp_ao_iw(:,:)
  double precision,allocatable  :: cNO(:,:)
  double precision,allocatable  :: cHFinv(:,:)
  double precision,allocatable  :: F_ao(:,:)
  double precision,allocatable  :: U_mo(:,:)
  double precision,allocatable  :: P_ao(:,:)
  double precision,allocatable  :: P_ao_hf(:,:)
  double precision,allocatable  :: P_ao_old(:,:)
  double precision,allocatable  :: P_ao_iter(:,:)
  double precision,allocatable  :: P_mo(:,:)
  double precision,allocatable  :: err_currentP(:)
  double precision,allocatable  :: err_diisP(:,:)
  double precision,allocatable  :: P_ao_extrap(:)
  double precision,allocatable  :: P_ao_old_diis(:,:)

  complex*16                    :: product
  complex*16                    :: weval_cpx
  complex*16                    :: EcGM_itau
!  complex*16                    :: EcGM_iw
  complex*16,allocatable        :: Sigma_c_w_ao(:,:,:)
  complex*16,allocatable        :: DeltaG_ao_iw(:,:,:)
  complex*16,allocatable        :: G_ao_itau(:,:,:)
  complex*16,allocatable        :: G_ao_itau_old(:,:,:)
  complex*16,allocatable        :: G_ao_itau_hf(:,:,:)
  complex*16,allocatable        :: G_ao_iw_hf(:,:,:)
  complex*16,allocatable        :: Sigma_c_c(:,:),Sigma_c_s(:,:)
  complex*16,allocatable        :: Sigma_c_plus(:,:),Sigma_c_minus(:,:)
  complex*16,allocatable        :: G_ao_1(:,:),G_ao_2(:,:)
  complex*16,allocatable        :: G_minus_itau(:,:),G_plus_itau(:,:)
  complex*16,allocatable        :: error_transf_mo(:,:,:)
  complex*16,allocatable        :: Sigma_c_w_mo(:,:)
  complex*16,allocatable        :: err_current(:)
  complex*16,allocatable        :: G_itau_extrap(:)
  complex*16,allocatable        :: err_diis(:,:)
  complex*16,allocatable        :: G_itau_old_diis(:,:)
  complex*16,allocatable        :: Chi0_ao_itau_vSq(:,:)
  complex*16,allocatable        :: Aimql(:,:,:,:)
  complex*16,allocatable        :: Bisql(:,:,:,:)
  complex*16,allocatable        :: Cispl(:,:,:,:)
!  complex*16,allocatable        :: Chi0_ao_iw_vSq(:,:,:)

! Output variables
  integer,intent(inout)         :: nfreqs
  double precision,intent(inout):: eHF(nOrb)
  double precision,intent(inout):: wcoord(nfreqs)
  double precision,intent(inout):: wweight(nfreqs)
  double precision,intent(inout):: cHF(nBas,nOrb)
  
!------------------------------------------------------------------------
! Build G(i tau) in AO basis and use it to build Xo (i tau) -> Xo (i w) !
!------------------------------------------------------------------------
 
 call wall_time(start_scGF2itauiw)

 write(*,*)     
 write(*,*)'*****************************************'
 write(*,*)'*    scGF2 ( using it and iw grids )    *'
 write(*,*)'*****************************************'
 write(*,*)

 read_SD_chkp=.false.
 n_diis=0
 verbose=0
 if(verbose_scGF2) verbose=1
 eta=0d0
 thrs_N=1d-10
 thrs_Ngrad=1d-6
 thrs_Pao=1d-6
 nElectrons=2d0*nO
 nBas2=nBas*nBas
 chem_pot_saved = 0.5d0*(eHF(nO)+eHF(nO+1))
 chem_pot = chem_pot_saved
 alpha_mixing=0.6d0
 rcond=0d0
 Ehfl=0d0
 write(*,*)
 write(*,'(A33,1X,F16.10,A3)') ' Initial chemical potential  = ',chem_pot,' au'
 write(*,*)
 eHF(:) = eHF(:)-chem_pot_saved
 if(verbose/=0) then
  write(*,*)
  write(*,*) ' Aligned HF energies from Go(iw) (a.u.) [ using HOMO-LUMO ]'
  do ibas=1,nOrb
   write(*,'(I7,F15.8)') ibas,eHF(ibas)
  enddo
  write(*,*)
 endif
 write(*,'(a,F15.8)') '  emin ',abs(eHF(nO))
 write(*,'(a,F15.8)') '  emax ',abs(eHF(nOrb)-eHF(1))
 write(*,*)
 inquire(file='read_SD_scGF2', exist=file_exists)
 if(file_exists) read_SD_chkp=.true.
  
 allocate(U_mo(nOrb,nOrb))
 allocate(P_ao(nBas,nBas),P_ao_old(nBas,nBas),P_ao_iter(nBas,nBas),P_ao_hf(nBas,nBas))
 allocate(F_ao(nBas,nBas),P_mo(nOrb,nOrb),cHFinv(nOrb,nBas),Occ(nOrb),eSD(nOrb),eSD_old(nOrb),cNO(nBas,nOrb))
 allocate(G_minus_itau(nBas,nBas),G_plus_itau(nBas,nBas)) 
 allocate(G_ao_1(nBas,nBas),G_ao_2(nBas,nBas)) 
 allocate(Sigma_c_c(nBas,nBas),Sigma_c_s(nBas,nBas)) 
 allocate(Sigma_c_plus(nBas,nBas),Sigma_c_minus(nBas,nBas))
 allocate(Aimql(nBas,nBas,nBas,nBAS))
 allocate(Bisql(nBas,nBas,nBas,nBAS))
 allocate(Cispl(nBas,nBas,nBas,nBAS))
 allocate(Chi0_ao_itau_vSq(nBas2,nBas2)) 
 cHFinv=matmul(transpose(cHF),S)
 P_ao_hf=P_in
 P_ao=P_in
 P_ao_iter=P_in
 F_ao=Hc
 Ehfl=0d0
 trace_1_rdm=0d0
 eSD_old(:)=eHF(:)
 do ibas=1,nBas
  do jbas=1,nBas
   Ehfl=Ehfl+P_ao(ibas,jbas)*Hc(ibas,jbas)
   trace_1_rdm=trace_1_rdm+P_ao(ibas,jbas)*S(ibas,jbas)
   do kbas=1,nBas
    do lbas=1,nBas
     F_ao(ibas,jbas)=F_ao(ibas,jbas)+P_ao(kbas,lbas)*ERI_AO(kbas,ibas,lbas,jbas) &
                    -0.5d0*P_ao(kbas,lbas)*ERI_AO(kbas,ibas,jbas,lbas)
     Ehfl=Ehfl+0.5d0*P_ao(kbas,lbas)*P_ao(ibas,jbas)*ERI_AO(kbas,ibas,lbas,jbas) &
         -0.25d0*P_ao(kbas,lbas)*P_ao(ibas,jbas)*ERI_AO(kbas,ibas,jbas,lbas)
    enddo
   enddo
  enddo
 enddo

!-----------------!
! Allocate arrays !
!-----------------!

 ntimes=nfreqs
 ntimes_twice=2*ntimes
 nBasSqntimes2=nBas2*ntimes_twice
 allocate(tweight(ntimes),tcoord(ntimes))
 allocate(sint2w_weight(nfreqs,ntimes))
 allocate(cost2w_weight(nfreqs,ntimes))
 allocate(cosw2t_weight(ntimes,nfreqs))
 allocate(sinw2t_weight(ntimes,nfreqs))
 allocate(Sigma_c_w_ao(nfreqs,nBas,nBas),DeltaG_ao_iw(nfreqs,nBas,nBas),G_ao_iw_hf(nfreqs,nBas,nBas))
 allocate(G_ao_itau(ntimes_twice,nBas,nBas),G_ao_itau_hf(ntimes_twice,nBas,nBas))
 allocate(G_ao_itau_old(ntimes_twice,nBas,nBas))
! allocate(Chi0_ao_iw_vSq(nfreqs,nBas2,nBas2))
 allocate(err_current(1))
 allocate(err_currentP(1))
 allocate(G_itau_extrap(1))
 allocate(P_ao_extrap(1))
 allocate(err_diis(1,1))
 allocate(err_diisP(1,1))
 allocate(G_itau_old_diis(1,1))
 allocate(P_ao_old_diis(1,1))
 if(maxDIIS>0) then
  deallocate(err_current)
  deallocate(err_currentP)
  deallocate(G_itau_extrap)
  deallocate(P_ao_extrap)
  deallocate(err_diis)
  deallocate(err_diisP)
  deallocate(G_itau_old_diis)
  deallocate(P_ao_old_diis)
  allocate(err_current(nBasSqntimes2))
  allocate(err_currentP(nBas2))
  allocate(G_itau_extrap(nBasSqntimes2))
  allocate(P_ao_extrap(nBas2))
  allocate(err_diis(nBasSqntimes2,maxDIIS))
  allocate(err_diisP(nBas2,maxDIIS))
  allocate(G_itau_old_diis(nBasSqntimes2,maxDIIS))
  allocate(P_ao_old_diis(nBas2,maxDIIS))
 endif

!---------------!
! Reading grids !
!---------------!

 call read_scGW_grids(ntimes,nfreqs,tcoord,tweight,wcoord,wweight,sint2w_weight,cost2w_weight, &
                      cosw2t_weight,sinw2t_weight,verbose)

!-------------------------------------------------------------------------!
! Test the quality of the grid for the Go(i w) -> G(i tau) transformation !
!-------------------------------------------------------------------------!

 ! Build Go(i w)
 write(*,*)
 write(*,'(a)') ' Error test for the Go(iw) -> G(it) transformation'
 write(*,*)
 do ifreq=1,nfreqs
  weval_cpx=im*wcoord(ifreq)
  call G_AO_RHF(nBas,nOrb,nO,eta,cHF,eHF,weval_cpx,G_ao_1)
  DeltaG_ao_iw(ifreq,:,:)=G_ao_1(:,:)
 enddo
 ! Fourier transform Go(i w) -> Go(i tau)
 G_ao_itau=czero
 do itau=1,ntimes
  G_plus_itau(:,:)=czero
  G_minus_itau(:,:)=czero
  do ifreq=1,nfreqs
   G_plus_itau(:,:) = G_plus_itau(:,:)   + im*cosw2t_weight(itau,ifreq)*Real(DeltaG_ao_iw(ifreq,:,:))  &
                                         - im*sinw2t_weight(itau,ifreq)*Aimag(DeltaG_ao_iw(ifreq,:,:)) 
   G_minus_itau(:,:) = G_minus_itau(:,:) + im*cosw2t_weight(itau,ifreq)*Real(DeltaG_ao_iw(ifreq,:,:))  &
                                         + im*sinw2t_weight(itau,ifreq)*Aimag(DeltaG_ao_iw(ifreq,:,:)) 
  enddo
  G_ao_itau(2*itau-1,:,:)=G_plus_itau(:,:) 
  G_ao_itau(2*itau  ,:,:)=G_minus_itau(:,:)
 enddo
 ! Check the error
 max_error_gw2gt=-1d0
 sum_error_gw2gt=0d0
 imax_error_gw2gt=1
 do itau=1,ntimes
  call G0itau_ao_RHF(nBas,nOrb,nO, tcoord(itau),G_plus_itau ,cHF,eHF)
  call G0itau_ao_RHF(nBas,nOrb,nO,-tcoord(itau),G_minus_itau,cHF,eHF)
  if(verbose/=0) then
   write(*,'(a,*(f20.8))') ' Fourier  ',im*tcoord(itau)
   do ibas=1,nBas
    write(*,'(*(f20.8))') G_ao_itau(2*itau-1,ibas,:)
   enddo
   write(*,'(a,*(f20.8))') ' Fourier  ',-im*tcoord(itau)
   do ibas=1,nBas
    write(*,'(*(f20.8))') G_ao_itau(2*itau  ,ibas,:)
   enddo
   write(*,'(a,*(f20.8))') ' Analytic  ',im*tcoord(itau)
   do ibas=1,nBas
    write(*,'(*(f20.8))') G_plus_itau(ibas,:)
   enddo
   write(*,'(a,*(f20.8))') ' Analytic  ',-im*tcoord(itau)
   do ibas=1,nBas
    write(*,'(*(f20.8))') G_minus_itau(ibas,:)
   enddo
  endif
  G_plus_itau(:,:) =abs(G_plus_itau(:,:) -G_ao_itau(2*itau-1,:,:))
  G_minus_itau(:,:)=abs(G_minus_itau(:,:)-G_ao_itau(2*itau  ,:,:))
  error_gw2gt=real(sum(G_plus_itau(:,:)))+real(sum(G_minus_itau(:,:)))
  sum_error_gw2gt=sum_error_gw2gt+error_gw2gt
  if(error_gw2gt>max_error_gw2gt) then
   imax_error_gw2gt=itau
   max_error_gw2gt=error_gw2gt
  endif
 enddo
 write(*,'(a,*(f20.8))') ' Sum error ',sum_error_gw2gt
 write(*,'(a,f20.8,a,2f20.8,a)') ' Max CAE   ',max_error_gw2gt,' is in the time +/-',0d0,tcoord(imax_error_gw2gt),'i'
 write(*,'(a,*(f20.8))') ' MAE       ',sum_error_gw2gt/(nfreqs*nBas*nBas)
 ! Reset to 0.0
 DeltaG_ao_iw=czero
 G_ao_itau=czero

!------------!
! scGF2 loop !
!----.-------!

 iter=0
 iter_fock=0
 do
  iter=iter+1

  ! For iter=1 we build G_ao_itau as the RHF one or read it from restart files
  ! [ we also initialize G_ao_iw_hf, G_ao_itau_hf, G_ao_itau_old, and (P_ao,P_ao_iter) ]
  if(iter==1) then
   G_ao_itau=czero
   do itau=1,ntimes
    call G0itau_ao_RHF(nBas,nOrb,nO, tcoord(itau),G_plus_itau ,cHF,eHF)
    call G0itau_ao_RHF(nBas,nOrb,nO,-tcoord(itau),G_minus_itau,cHF,eHF)
    G_ao_itau(2*itau-1,:,:)=G_plus_itau(:,:)
    G_ao_itau(2*itau  ,:,:)=G_minus_itau(:,:)
   enddo
   G_ao_itau_old(:,:,:)=G_ao_itau(:,:,:)
   G_ao_itau_hf(:,:,:)=G_ao_itau(:,:,:)
   do ifreq=1,nfreqs
    weval_cpx=im*wcoord(ifreq)
    call G_AO_RHF(nBas,nOrb,nO,eta,cHF,eHF,weval_cpx,G_ao_1)
    G_ao_iw_hf(ifreq,:,:)=G_ao_1(:,:)
   enddo
   ! Initialize DeltaG(i w) [ it will be G(i w) - Go(i w) ]
   DeltaG_ao_iw(:,:,:)=czero
   ! If required, read the restart files
   if(restart_scGF2) then
    call read_scGW_restart(nBas,nfreqs,ntimes_twice,chem_pot,P_ao,P_ao_hf,G_ao_iw_hf,G_ao_itau,G_ao_itau_hf,read_SD_chkp)
    P_ao_iter=P_ao
    G_ao_itau_old(:,:,:)=G_ao_itau(:,:,:)
   endif
  endif

  ! Build Sigma_c(i w) [Eqs. 12-18 in PRB, 109, 255101 (2024)] (for the M^8 algorithm see at the end of this file)
  Sigma_c_w_ao=czero
  do itau=1,ntimes
   G_plus_itau(:,:) =G_ao_itau(2*itau-1,:,:)
   G_minus_itau(:,:)=G_ao_itau(2*itau  ,:,:)
   Sigma_c_plus=czero
   Sigma_c_minus=czero
   ! Sigma_c_ij(i tau) = \sum_klmspq Gkl(i tau) Gms(i tau) Gpq(-i tau) v_iqmk (2 v_lspj - v_slpj)
   Aimql=czero
   Bisql=czero
   Cispl=czero
   do ibas=1,nBas
    do mbas=1,nBas
     do qbas=1,nBas
      do lbas=1,nBas
       do kbas=1,nBas
        Aimql(ibas,mbas,qbas,lbas)=Aimql(ibas,mbas,qbas,lbas)+G_plus_itau(kbas,lbas)*ERI_AO(ibas,qbas,mbas,kbas)
       enddo
      enddo
     enddo
    enddo
   enddo
   do ibas=1,nBas
    do sbas=1,nBas
     do qbas=1,nBas
      do lbas=1,nBas
       do mbas=1,nBas
        Bisql(ibas,sbas,qbas,lbas)=Bisql(ibas,sbas,qbas,lbas)+G_plus_itau(mbas,sbas)*Aimql(ibas,mbas,qbas,lbas)
       enddo
      enddo
     enddo
    enddo
   enddo
   do ibas=1,nBas
    do sbas=1,nBas
     do pbas=1,nBas
      do lbas=1,nBas
       do qbas=1,nBas
        Cispl(ibas,sbas,pbas,lbas)=Cispl(ibas,sbas,pbas,lbas)+G_minus_itau(pbas,qbas)*Bisql(ibas,sbas,qbas,lbas)
       enddo
      enddo
     enddo
    enddo
   enddo
   do ibas=1,nBas
    do jbas=1,nBas
     do sbas=1,nBas
      do pbas=1,nBas
       do lbas=1,nBas
        Sigma_c_plus(ibas,jbas) =Sigma_c_plus(ibas,jbas) +Cispl(ibas,sbas,pbas,lbas)*(2d0*ERI_AO(lbas,sbas,pbas,jbas)-ERI_AO(sbas,lbas,pbas,jbas))
       enddo
      enddo
     enddo
    enddo
   enddo
   ! Sigma_c_ij(-i tau) =  \sum_klmspq Gkl(-i tau) Gms(-i tau) Gpq(i tau) v_iqmk (2 v_lspj - v_slpj)
   Aimql=czero
   Bisql=czero
   Cispl=czero
   do ibas=1,nBas
    do mbas=1,nBas
     do qbas=1,nBas
      do lbas=1,nBas
       do kbas=1,nBas
        Aimql(ibas,mbas,qbas,lbas)=Aimql(ibas,mbas,qbas,lbas)+G_minus_itau(kbas,lbas)*ERI_AO(ibas,qbas,mbas,kbas)
       enddo
      enddo
     enddo
    enddo
   enddo
   do ibas=1,nBas
    do sbas=1,nBas
     do qbas=1,nBas
      do lbas=1,nBas
       do mbas=1,nBas
        Bisql(ibas,sbas,qbas,lbas)=Bisql(ibas,sbas,qbas,lbas)+G_minus_itau(mbas,sbas)*Aimql(ibas,mbas,qbas,lbas)
       enddo
      enddo
     enddo
    enddo
   enddo
   do ibas=1,nBas
    do sbas=1,nBas
     do pbas=1,nBas
      do lbas=1,nBas
       do qbas=1,nBas
        Cispl(ibas,sbas,pbas,lbas)=Cispl(ibas,sbas,pbas,lbas)+G_plus_itau(pbas,qbas)*Bisql(ibas,sbas,qbas,lbas)
       enddo
      enddo
     enddo
    enddo
   enddo
   do ibas=1,nBas
    do jbas=1,nBas
     do sbas=1,nBas
      do pbas=1,nBas
       do lbas=1,nBas
        Sigma_c_minus(ibas,jbas)=Sigma_c_minus(ibas,jbas) +Cispl(ibas,sbas,pbas,lbas)*(2d0*ERI_AO(lbas,sbas,pbas,jbas)-ERI_AO(sbas,lbas,pbas,jbas))
       enddo
      enddo
     enddo
    enddo
   enddo
   ! Corrected Eqs. 17 and 18 in PRB, 109, 245101 (2024)
   Sigma_c_c= -im*(Sigma_c_plus+Sigma_c_minus)
   Sigma_c_s= -   (Sigma_c_plus-Sigma_c_minus)
   ! Sigma_c(i tau) -> Sigma_c(i w)
   do ifreq=1,nfreqs
    Sigma_c_w_ao(ifreq,:,:) = Sigma_c_w_ao(ifreq,:,:)                        &
                            + 0.5d0*cost2w_weight(ifreq,itau)*Sigma_c_c(:,:) &
                            + 0.5d0*sint2w_weight(ifreq,itau)*Sigma_c_s(:,:)
   enddo 
  enddo

  ! Compute EcGM [ from Sigma_c(iw) is BAD. We use Eq. 10 from J. Chem. Theory Comput., 10, 2498 to compute it from Xo ]
  EcGM=0d0
!  EcGMw=0d0
!  Chi0_ao_iw_vSq(:,:,:)=czero
  ! Build using the time grid Xo(i tau) = -2i G(i tau) G(-i tau)
  do itau=1,ntimes
   ! Xo(i tau) = -2i G(i tau) G(-i tau)
   do ibas=1,nBas
    do jbas=1,nBas
     do kbas=1,nBas
      do lbas=1,nBas                       
                                   ! r1   r2'                    r2   r1'
       product = G_ao_itau(2*itau-1,ibas,jbas)*G_ao_itau(2*itau,kbas,lbas)
       if(abs(product)<1e-12) product=czero
       Chi0_ao_itau_vSq(1+(lbas-1)+(ibas-1)*nBas,1+(kbas-1)+(jbas-1)*nBas) = product
      enddo
     enddo
    enddo
   enddo
   Chi0_ao_itau_vSq=-2d0*im*Chi0_ao_itau_vSq ! The 2 factor is added to account for both spin contributions [ i.e., (up,up,up,up) and (down,down,down,down) ]
!   ! Xo(i tau) -> Xo(i w) [ the weight already contains the cos(tau w) and a factor 2 because int_-Infty ^Infty -> 2 int_0 ^Infty ]
!   do ifreq=1,nfreqs
!    Chi0_ao_iw_vSq(ifreq,:,:) = Chi0_ao_iw_vSq(ifreq,:,:) - im*cost2w_weight(ifreq,itau)*Chi0_ao_itau_vSq(:,:)
!   enddo 
   Chi0_ao_itau_vSq=matmul(Chi0_ao_itau_vSq,vMAT)              ! Xo(i tau) v
   Chi0_ao_itau_vSq=matmul(Chi0_ao_itau_vSq,Chi0_ao_itau_vSq)  ! [ Xo(i tau) v ]^2
   ! EcGM = 1/4 int Tr{ [ Xo(i tau) v ]^2 }
   EcGM_itau=czero
   do ibas=1,nBas2
    EcGM_itau=EcGM_itau+Chi0_ao_itau_vSq(ibas,ibas)
   enddo
   EcGM=EcGM+0.25d0*tweight(itau)*real(EcGM_itau) 
  enddo
!  ! Complete the Xo(i tau) -> Xo(i w)
!  Chi0_ao_iw_vSq(:,:,:) = Real(Chi0_ao_iw_vSq(:,:,:)) ! The factor 2 is stored in the weight [ and we just retain the real part ]
!  do ifreq=1,nfreqs
!   Chi0_ao_iw_vSq(ifreq,:,:)=matmul(Chi0_ao_iw_vSq(ifreq,:,:),vMAT)                       ! Xo(i w) v
!   Chi0_ao_iw_vSq(ifreq,:,:)=matmul(Chi0_ao_iw_vSq(ifreq,:,:),Chi0_ao_iw_vSq(ifreq,:,:))  ! [ Xo(i w) v ]^2
!   ! EcGM = - 1/8pi int Tr{ [ Xo(i tau) v ]^2 }
!   EcGM_iw=czero
!   do ibas=1,nBas2
!    EcGM_iw=EcGM_iw+Chi0_ao_iw_vSq(ifreq,ibas,ibas)
!   enddo
!   EcGMw=EcGMw-wweight(ifreq)*real(EcGM_iw)
!  enddo
!  EcGMw=EcGMw/(8d0*pi)

  ! Check the error in Sigma_c(i w) at iter=1 [ if this is calc. is not with restart ]
  if(iter==1 .and. .not.restart_scGF2) then
   write(*,*)
   write(*,'(a)') ' Error test for the Sigma_c(iw) construction at iter 1 [ compared with the analytic Sigma_c(iw) obtained from HF ] '
   write(*,*)
   max_error_sigma=-1d0;imax_error_sigma=1;
   allocate(error_transf_mo(nfreqs,nOrb,nOrb),Sigma_c_w_mo(nOrb,nOrb))
   ! Build the analytic Sigma_c(iw)
   call build_analityc_rhf_Sigma_c_iw_GF2(nBas,nOrb,nO,verbose,cHF,eHF,nfreqs,wcoord,ERI_AO,error_transf_mo,err_EcGM) ! error_transf_mo set to Sigma_c_mo(iw)
   do ifreq=1,nfreqs
    Sigma_c_w_mo=matmul(matmul(transpose(cHF(:,:)),Sigma_c_w_ao(ifreq,:,:)),cHF(:,:)) ! Fourier: Sigma_c_ao(iw) -> Sigma_c_mo(iw)
    !Sigma_c_w_ao(ifreq,:,:)=matmul(transpose(cHFinv),matmul(error_transf_mo(ifreq,:,:),cHFinv)) ! Analytic: Sigma_c_mo(iw) -> Sigma_c_ao(iw)
    if(verbose/=0) then
     write(*,'(a,*(f20.8))') ' Fourier  ',im*wcoord(ifreq)
     do ibas=1,nOrb
      write(*,'(*(f20.8))') Sigma_c_w_mo(ibas,:)
     enddo
    endif
    error_transf_mo(ifreq,:,:)=abs(error_transf_mo(ifreq,:,:)-Sigma_c_w_mo(:,:))
    error_sigma=real(sum(error_transf_mo(ifreq,:,:)))
    if(error_sigma>max_error_sigma) then
     imax_error_sigma=ifreq
     max_error_sigma=error_sigma
    endif
   enddo
   write(*,'(a,*(f20.8))') ' Sum error ',sum(error_transf_mo)
   write(*,'(a,f20.8,a,2f20.8,a)') ' Max CAE   ',max_error_sigma,' is in the frequency ',0d0,wcoord(imax_error_sigma),'i'
   write(*,'(a,*(f20.8))') ' MAE       ',sum(error_transf_mo)/(nfreqs*nBas*nBas)
   write(*,'(a)')         ' Using Xo(it) ' 
   write(*,'(a,f17.8,a)') ' EcGM analytic',err_EcGM,' a.u.'
   write(*,'(a,f18.8,a)') ' EcGM numeric',EcGM,' a.u.'
   write(*,'(a,f20.8,a)') ' EcGM error',abs(err_EcGM-EcGM),' a.u.'
!   write(*,'(a)')         ' Using Xo(iw) ' 
!   write(*,'(a,f18.8,a)') ' EcGM numeric',EcGMw,' a.u.'
!   write(*,'(a,f20.8,a)') ' EcGM error',abs(err_EcGM-EcGMw),' a.u.'
   deallocate(error_transf_mo,Sigma_c_w_mo)
  endif

  ! Converge with respect to the Fock operator (using only good P_ao matrices -> Tr[P_ao S_ao]=Nelectrons )
  if(.not.no_fock) then ! Skiiping the opt w.r.t. the Fock operator to do later the linearized approximation on Go -> [ lin-G = Go + Go Sigma Go ]
   iter_fock=0
   n_diisP=0
   rcondP=0d0
   err_diisP=0d0
   P_ao_old_diis=0d0
   do
    ! Build F
    iter_fock=iter_fock+1
    F_ao=Hc
    Ehfl=0d0
    do ibas=1,nBas
     do jbas=1,nBas
      Ehfl=Ehfl+P_ao(ibas,jbas)*Hc(ibas,jbas)
      do kbas=1,nBas
       do lbas=1,nBas
        F_ao(ibas,jbas)=F_ao(ibas,jbas)+P_ao(kbas,lbas)*ERI_AO(kbas,ibas,lbas,jbas) &
                       -0.5d0*P_ao(kbas,lbas)*ERI_AO(kbas,ibas,jbas,lbas)
        Ehfl=Ehfl+0.5d0*P_ao(kbas,lbas)*P_ao(ibas,jbas)*ERI_AO(kbas,ibas,lbas,jbas) &
            -0.25d0*P_ao(kbas,lbas)*P_ao(ibas,jbas)*ERI_AO(kbas,ibas,jbas,lbas)
       enddo
      enddo
     enddo
    enddo
    ! Build G(i w) and n(r)
    P_ao_old=P_ao
    call get_1rdm_scGW(nBas,nfreqs,chem_pot,S,F_ao,Sigma_c_w_ao,wcoord,wweight, &
                       G_ao_1,G_ao_iw_hf,DeltaG_ao_iw,P_ao,P_ao_hf,trace_1_rdm) 
    if(abs(trace_1_rdm-nElectrons)**2d0>thrs_N) &
     call fix_chem_pot_scGW_bisec(iter_fock,nBas,nfreqs,nElectrons,thrs_N,thrs_Ngrad,chem_pot,S,F_ao,Sigma_c_w_ao,wcoord,wweight, &
                                  G_ao_1,G_ao_iw_hf,DeltaG_ao_iw,P_ao,P_ao_hf,trace_1_rdm,chem_pot_saved,verbose_scGF2)
    ! Check convergence of P_ao for fixed Sigma_c(i w)
    diff_Pao=0d0
    do ibas=1,nBas
     do jbas=1,nBas
      diff_Pao=diff_Pao+abs(P_ao(ibas,jbas)-P_ao_old(ibas,jbas))
     enddo
    enddo

    if(diff_Pao<=thrs_Pao) exit
   
    if(iter_fock==maxSCF) exit
 
    ! Do mixing with previous P_ao to facilitate convergence
    if(maxDIIS>0) then
     n_diisP=min(n_diisP+1,maxDIIS)
     err_currentP=0d0
     idiis_indexP=1
     do ibas=1,nBas
      do jbas=1,nBas
       err_currentP(idiis_indexP)=P_ao(ibas,jbas)-P_ao_old(ibas,jbas)
       P_ao_extrap(idiis_indexP)=P_ao(ibas,jbas)
       idiis_indexP=idiis_indexP+1
      enddo
     enddo
     call DIIS_extrapolation(rcondP,nBas2,nBas2,n_diisP,err_diisP,P_ao_old_diis,err_currentP,P_ao_extrap)
     idiis_indexP=1
     do ibas=1,nBas
      do jbas=1,nBas
       P_ao(ibas,jbas)=P_ao_extrap(idiis_indexP) 
       idiis_indexP=idiis_indexP+1
      enddo
     enddo
    endif
   
   enddo
  endif

  ! Check convergence of P_ao after a scGF2 iteration
  diff_Pao=0d0
  do ibas=1,nBas
   do jbas=1,nBas
    diff_Pao=diff_Pao+abs(P_ao(ibas,jbas)-P_ao_iter(ibas,jbas))
   enddo
  enddo
  P_ao_iter=P_ao

  ! Print iter info
  P_mo=-matmul(matmul(cHFinv,P_ao),transpose(cHFinv)) ! Minus to order occ numbers
  call diagonalize_matrix(nOrb,P_mo,Occ)
  Occ=-Occ
  cNO=matmul(cHF,P_mo)
  write(*,*)
  write(*,'(a,f15.8,a,i5,a,i5)') ' Trace scGF2 ',trace_1_rdm,' after ',iter_fock,' Fock iterations at global iter ',iter
  write(*,'(a,f15.8)')        ' Change of P ',diff_Pao
  write(*,'(a,f15.8)')        ' Chem. Pot.  ',chem_pot
  write(*,'(a,f15.8)')        ' Enuc        ',ENuc
  write(*,'(a,f15.8)')        ' Ehfl        ',Ehfl
  write(*,'(a,f15.8)')        ' EcGM(it)    ',EcGM
  write(*,'(a,f15.8)')        ' Eelec(it)   ',Ehfl+EcGM
  write(*,'(a,f15.8)')        ' Etot(it)    ',Ehfl+EcGM+ENuc
!  write(*,'(a,f15.8)')        ' EcGM(iw)    ',EcGMw
!  write(*,'(a,f15.8)')        ' Eelec(iw)   ',Ehfl+EcGMw
!  write(*,'(a,f15.8)')        ' Etot(iw)    ',Ehfl+EcGMw+ENuc
  write(*,*)

  if(diff_Pao<=thrs_Pao) exit

  if(iter==maxSCF) exit

  ! Build the new G_ao_iw_hf, G_ao_itau_hf, and P_ao_hf
  U_mo=matmul(transpose(cHF),matmul(F_ao,cHF))
  do ibas=1,nBas
   U_mo(ibas,ibas)=U_mo(ibas,ibas)-chem_pot
  enddo
  call diagonalize_matrix(nOrb,U_mo,eSD)
  chem_pot_align=0.5d0*(eSD(nO)+eSD(nO+1))
  if(verbose/=0) then
   write(*,*) '    orb       Occ        SD energies  Aligned SD energies [ from Go(iw) (a.u.) ]'
   do ibas=1,nOrb
    write(*,'(I7,3F15.8)') ibas,Occ(ibas),eSD(ibas),eSD(ibas)-chem_pot_align
   enddo
  endif
  eSD(:)=eSD(:)-chem_pot_align
  nneg=0
  do ibas=1,nOrb
   if(eSD(ibas)<0d0) nneg=nneg+1
  enddo
  sd_dif=sum(abs(eSD(:)-eSD_old(:)))
  write(*,'(a,f15.8)') '     | eSD,i - eSD,i-1 | ',sd_dif
  if(nneg==nO .and. sd_dif>1d-2) then
   write(*,*)
   write(*,'(a,i5)') ' Computing new Go(iw), Go(it), and P_HF matrices at global iter ',iter
   write(*,*)
   eSD_old(:)=eSD(:)
   ! Compute new MO coefs
   cHF=matmul(cHF,U_mo)
   cHFinv=matmul(transpose(cHF),S)
   ! New P_ao_hf
   P_ao_hf(:,:) = 2d0*matmul(cHF(:,1:nO),transpose(cHF(:,1:nO)))
   ! New G_ao_itau_hf
   G_ao_itau_hf=czero
   do itau=1,ntimes
    call G0itau_ao_RHF(nBas,nOrb,nO, tcoord(itau),G_plus_itau ,cHF,eSD)
    call G0itau_ao_RHF(nBas,nOrb,nO,-tcoord(itau),G_minus_itau,cHF,eSD)
    G_ao_itau_hf(2*itau-1,:,:)=G_plus_itau(:,:)
    G_ao_itau_hf(2*itau  ,:,:)=G_minus_itau(:,:)
   enddo
   ! New G_ao_iw_hf [ Go_new(iw) ]
   DeltaG_ao_iw(:,:,:)=G_ao_iw_hf(:,:,:)+DeltaG_ao_iw(:,:,:) ! Saving G(iw) in DeltaG_ao_iw
   do ifreq=1,nfreqs
    weval_cpx=im*wcoord(ifreq)
    call G_AO_RHF(nBas,nOrb,nO,eta,cHF,eSD,weval_cpx,G_ao_1)
    G_ao_iw_hf(ifreq,:,:)=G_ao_1(:,:)
   enddo
   DeltaG_ao_iw(:,:,:)=DeltaG_ao_iw(:,:,:)-G_ao_iw_hf(:,:,:) ! Setting back DeltaG(iw) = G(iw) - Go_new(iw)
  endif

  ! Transform DeltaG(i w) -> DeltaG(i tau) [ i tau and -i tau ]
  !      [ the weights contain the 2 /(2 pi) = 1 / pi factor and the cos(tau w) or sin(tau w) ]
  G_ao_itau=czero
  do itau=1,ntimes
   G_plus_itau(:,:)=czero
   G_minus_itau(:,:)=czero
   do ifreq=1,nfreqs
    G_plus_itau(:,:) = G_plus_itau(:,:)   + im*cosw2t_weight(itau,ifreq)*Real(DeltaG_ao_iw(ifreq,:,:))  &
                                          - im*sinw2t_weight(itau,ifreq)*Aimag(DeltaG_ao_iw(ifreq,:,:)) 
    G_minus_itau(:,:) = G_minus_itau(:,:) + im*cosw2t_weight(itau,ifreq)*Real(DeltaG_ao_iw(ifreq,:,:))  &
                                          + im*sinw2t_weight(itau,ifreq)*Aimag(DeltaG_ao_iw(ifreq,:,:)) 
   enddo
   ! Build G(i tau) = DeltaG(i tau) + Go(i tau)
   G_ao_itau(2*itau-1,:,:)=G_plus_itau(:,:) +G_ao_itau_hf(2*itau-1,:,:)
   G_ao_itau(2*itau  ,:,:)=G_minus_itau(:,:)+G_ao_itau_hf(2*itau  ,:,:)
  enddo
 
  ! Do mixing with previous G(i tau) to facilitate convergence
  if(maxDIIS>0) then
   n_diis=min(n_diis+1,maxDIIS)
   err_current=czero
   idiis_index=1
   do itau=1,ntimes_twice
    do ibas=1,nBas
     do jbas=1,nBas
      err_current(idiis_index)=G_ao_itau(itau,ibas,jbas)-G_ao_itau_old(itau,ibas,jbas)
      G_itau_extrap(idiis_index)=G_ao_itau(itau,ibas,jbas)
      idiis_index=idiis_index+1
     enddo
    enddo
   enddo
   call complex_DIIS_extrapolation(rcond,nBasSqntimes2,nBasSqntimes2,n_diis,err_diis,G_itau_old_diis,err_current,G_itau_extrap)
   idiis_index=1
   do itau=1,ntimes_twice
    do ibas=1,nBas
     do jbas=1,nBas
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
 write(*,*)
 write(*,'(A50)') '---------------------------------------'
 write(*,'(A50)') '     scGF2 calculation completed       '
 write(*,'(A50)') '---------------------------------------'
 write(*,*)
 write(*,'(a,f15.8,a,i5,a)') ' Trace scGF2      ',trace_1_rdm,' after ',iter,' global iterations '
 write(*,'(a,f15.8)')        ' Change of P      ',diff_Pao
 write(*,'(a,f15.8)')        ' Chem. Pot.       ',chem_pot
 write(*,'(a,f15.8)')        ' Enuc             ',ENuc
 write(*,'(a,f15.8)')        ' Ehfl             ',Ehfl
 write(*,'(a,f15.8)')        ' EcGM(it)         ',EcGM
 write(*,'(a,f15.8)')        ' Eelec(it)        ',Ehfl+EcGM
 write(*,'(a,f15.8)')        ' scGF2 Energy (it)',Ehfl+EcGM+ENuc
! write(*,'(a,f15.8)')        ' EcGM(iw)         ',EcGMw
! write(*,'(a,f15.8)')        ' Eelec(iw)        ',Ehfl+EcGMw
! write(*,'(a,f15.8)')        ' Etot(iw)         ',Ehfl+EcGMw+ENuc
! write(*,'(a,f15.8)')        ' scGF2 Energy (iw)',Ehfl+EcGMw+ENuc
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
 call write_scGW_restart(nBas,ntimes,ntimes_twice,nfreqs,chem_pot,P_ao,P_ao_hf,G_ao_itau,G_ao_itau_hf, &
                         G_ao_iw_hf,DeltaG_ao_iw)
 
 ! Using the correlated G and Sigma_c to test the linearized density matrix approximation
 if(dolinGF2) then
  write(*,*)
  write(*,*) ' -------------------------------------------'
  write(*,*) ' Testing the linearized approximation with G'
  write(*,*) '         G^lin = G + G Sigma_c G'
  write(*,*) ' -------------------------------------------'
  P_ao_old=0d0
  G_ao_1(:,:)=czero
  do ifreq=1,nfreqs
   G_ao_1(:,:)=G_ao_iw_hf(ifreq,:,:)+DeltaG_ao_iw(ifreq,:,:)
   G_ao_1(:,:)=matmul(matmul(G_ao_1(:,:),Sigma_c_w_ao(ifreq,:,:)),G_ao_1(:,:))
   P_ao_old(:,:) = P_ao_old(:,:) + wweight(ifreq)*real(G_ao_1(:,:)+conjg(G_ao_1(:,:))) ! Integrate along iw
  enddo
  P_ao_old=P_ao_old/pi
  P_ao_old=P_ao+P_ao_old
  trace_1_rdm=0d0
  do ibas=1,nBas
   do jbas=1,nBas
    trace_1_rdm=trace_1_rdm+P_ao_old(ibas,jbas)*S(ibas,jbas)
   enddo
  enddo
  Ehfl=0d0
  do ibas=1,nBas
   do jbas=1,nBas
    Ehfl=Ehfl+P_ao_old(ibas,jbas)*Hc(ibas,jbas)
    do kbas=1,nBas
     do lbas=1,nBas
      Ehfl=Ehfl+0.5d0*P_ao(kbas,lbas)*P_ao(ibas,jbas)*ERI_AO(kbas,ibas,lbas,jbas) &
          -0.25d0*P_ao(kbas,lbas)*P_ao(ibas,jbas)*ERI_AO(kbas,ibas,jbas,lbas)
     enddo
    enddo
   enddo
  enddo
  P_mo=-matmul(matmul(cHFinv,P_ao_old),transpose(cHFinv)) ! Minus to order occ numbers
  call diagonalize_matrix(nOrb,P_mo,Occ)
  Occ=-Occ
  cNO=matmul(cHF,P_mo)
  write(*,'(a,f15.8)')        ' Enuc        ',ENuc
  write(*,'(a,f15.8)')        ' Ehfl        ',Ehfl
  write(*,'(a,f15.8)')        ' EcGM        ',EcGM
  write(*,'(a,f15.8)')        ' Eelec       ',Ehfl+EcGM
  write(*,'(a,f15.8)')        ' lin-G Energy',Ehfl+EcGM+ENuc
  write(*,*)
  write(*,'(a,f15.8,a,i5,a)') ' Trace lin-scGF2 ',trace_1_rdm
  write(*,*)
  write(*,*) ' Lin-G occupation numbers'
  do ibas=1,nOrb
   write(*,'(I7,F15.8)') ibas,Occ(ibas)
  enddo
  if(verbose/=0) then
   write(*,*) ' Natural orbitals (columns)'
   do ibas=1,nBas
    write(*,'(*(f15.8))') cNO(ibas,:)
   enddo
  endif
 endif

 call wall_time(end_scGF2itauiw)
 
 t_scGF2itauiw = end_scGF2itauiw - start_scGF2itauiw
 write(*,'(A65,1X,F9.3,A8)') 'Total wall time for scGF2 =',t_scGF2itauiw,' seconds'
 write(*,*)

 ! Restore values and deallocate dyn arrays
 eHF(:) = eHF(:)+chem_pot_saved
 deallocate(tcoord,tweight) 
 deallocate(sint2w_weight)
 deallocate(cost2w_weight)
 deallocate(cosw2t_weight)
 deallocate(sinw2t_weight)
 deallocate(G_ao_itau_old)
 deallocate(G_ao_itau,G_ao_itau_hf)
 deallocate(Sigma_c_w_ao,DeltaG_ao_iw,G_ao_iw_hf)
 deallocate(P_ao,P_ao_old,P_ao_iter,P_ao_hf,F_ao,P_mo,cHFinv,cNO,U_mo,Occ,eSD,eSD_old) 
 deallocate(Sigma_c_plus,Sigma_c_minus) 
 deallocate(Sigma_c_c,Sigma_c_s) 
 deallocate(G_minus_itau,G_plus_itau) 
 deallocate(G_ao_1,G_ao_2) 
 deallocate(err_current)
 deallocate(err_currentP)
 deallocate(G_itau_extrap)
 deallocate(P_ao_extrap)
 deallocate(err_diis)
 deallocate(err_diisP)
 deallocate(G_itau_old_diis)
 deallocate(P_ao_old_diis)
 deallocate(Aimql)
 deallocate(Bisql)
 deallocate(Cispl)
 deallocate(Chi0_ao_itau_vSq) 
! deallocate(Chi0_ao_iw_vSq)

end subroutine 

!  double precision              :: ERI_contrib
   ! M^8 brut force
!   do ibas=1,nBas
!    do jbas=1,nBas
!     do kbas=1,nBas
!      do lbas=1,nBas 
!       do mbas=1,nBas 
!        do sbas=1,nBas 
!         do pbas=1,nBas 
!          do qbas=1,nBas
!           ERI_contrib=ERI_AO(ibas,qbas,mbas,kbas)*(2d0*ERI_AO(lbas,sbas,pbas,jbas)-ERI_AO(sbas,lbas,pbas,jbas)) 
!           Sigma_c_plus(ibas,jbas) =Sigma_c_plus(ibas,jbas) + G_plus_itau(kbas,lbas)* G_plus_itau(mbas,sbas) &
!                                                            *G_minus_itau(pbas,qbas)*ERI_contrib
!           Sigma_c_minus(ibas,jbas)=Sigma_c_minus(ibas,jbas)+G_minus_itau(kbas,lbas)*G_minus_itau(mbas,sbas) &
!                                                            * G_plus_itau(pbas,qbas)*ERI_contrib
!          enddo
!         enddo
!        enddo
!       enddo
!      enddo
!     enddo
!    enddo
!   enddo

