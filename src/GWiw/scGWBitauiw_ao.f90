subroutine scGWBitauiw_ao(nBas,nOrb,nOrb_twice,maxSCF,maxDIIS,dolinGW,restart_scGWB,verbose_scGWB,no_fock,ENuc,Hc,S,X_in, &
                          P_in,Pan_in,cHFB,eQP_state,chem_pot,sigma,nfreqs,wcoord,wweight,U_QP,vMAT,ERI_AO)

! Restricted scGWB

  implicit none
  include 'parameters.h'

! Input variables
 
  logical,intent(in)            :: dolinGW
  logical,intent(in)            :: no_fock
  logical,intent(in)            :: restart_scGWB
  logical,intent(in)            :: verbose_scGWB

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nOrb_twice
  integer,intent(in)            :: maxSCF
  integer,intent(in)            :: maxDIIS

  double precision,intent(in)   :: ENuc
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

  integer                       :: ifreq,itau
  integer                       :: ibas,jbas,kbas,lbas,mbas,obas,pbas,qbas
  integer                       :: verbose
  integer                       :: ntimes
  integer                       :: nBasSq
  integer                       :: nBas_twice
  integer                       :: ntimes_twice
  integer                       :: imax_error_gw2gt

  double precision              :: eta
  double precision              :: error_gw2gt
  double precision              :: max_error_gw2gt
  double precision              :: sum_error_gw2gt
  double precision              :: EcGM,Ehfbl,Ecore,Eh,Ex,Epair
  double precision              :: trace1,trace2
  double precision              :: trace_1_rdm
  double precision              :: start_scGWBitauiw     ,end_scGWBitauiw       ,t_scGWBitauiw
  double precision,allocatable  :: Occ(:)
  double precision,allocatable  :: cHFBinv(:,:)
  double precision,allocatable  :: cNO(:,:)
  double precision,allocatable  :: U_mo(:,:)
  double precision,allocatable  :: Mat1(:,:)
  double precision,allocatable  :: Mat2(:,:)
  double precision,allocatable  :: Wp_ao_iw(:,:)
  double precision,allocatable  :: Wp_ao_itau(:,:,:)
  double precision,allocatable  :: R_ao(:,:)
  double precision,allocatable  :: R_ao_old(:,:)
  double precision,allocatable  :: tweight(:),tcoord(:)
  double precision,allocatable  :: sint2w_weight(:,:)
  double precision,allocatable  :: cost2w_weight(:,:)
  double precision,allocatable  :: cosw2t_weight(:,:)
  double precision,allocatable  :: sinw2t_weight(:,:)

  complex*16                    :: product
  complex*16                    :: weval_cpx
  complex*16,allocatable        :: Sigma_c_w_ao(:,:,:)
  complex*16,allocatable        :: Sigma_c_c(:,:),Sigma_c_s(:,:)
  complex*16,allocatable        :: Sigma_c_plus(:,:),Sigma_c_minus(:,:)
  complex*16,allocatable        :: G_ao_tmp(:,:)
  complex*16,allocatable        :: G_gorkov_tmp(:,:)
  complex*16,allocatable        :: G_ao_itau(:,:,:)
  complex*16,allocatable        :: G_ao_itau_hfb(:,:,:)
  complex*16,allocatable        :: G_ao_iw_hfb(:,:,:)
  complex*16,allocatable        :: DeltaG_ao_iw(:,:,:)
  complex*16,allocatable        :: Chi0_ao_iw(:,:,:)
  complex*16,allocatable        :: Chi0_ao_itau(:,:)

! Output variables
  integer,intent(inout)         :: nfreqs
  double precision,intent(inout):: chem_pot
  double precision,intent(inout):: eQP_state(nOrb_twice)
  double precision,intent(inout):: wcoord(nfreqs)
  double precision,intent(inout):: wweight(nfreqs)
  double precision,intent(inout):: cHFB(nBas,nOrb)
  
 call wall_time(start_scGWBitauiw)

 write(*,*)     
 write(*,*)'*****************************************'
 write(*,*)'*    scGWB ( using it and iw grids )    *'
 write(*,*)'*****************************************'
 write(*,*)

 write(*,*)
 write(*,'(A33,1X,F16.10,A3)') ' Initial chemical potential  = ',chem_pot,' au'
 write(*,*)

 write(*,'(a,F15.8)') '  emin ',abs(eQP_state(nOrb))
 write(*,'(a,F15.8)') '  emax ',abs(eQP_state(nOrb)+eQP_state(1))

 ! Initialize variables
 ntimes=nfreqs
 nBasSq=nBas*nBas
 ntimes_twice=2*ntimes
 nBas_twice=2*nBas
 verbose=0
 if(verbose_scGWB) verbose=1     
 eta=0d0


 ! Allocate arrays
 allocate(Occ(nOrb))
 allocate(cNO(nBas,nOrb))
 allocate(cHFBinv(nOrb,nBas))
 allocate(U_mo(nOrb,nOrb))
 allocate(R_ao(nBas_twice,nBas_twice))
 allocate(R_ao_old(nBas_twice,nBas_twice))
 allocate(G_ao_tmp(nBas,nBas))
 allocate(G_gorkov_tmp(nBas_twice,nBas_twice))
 allocate(G_ao_iw_hfb(nfreqs,nBas_twice,nBas_twice))
 allocate(DeltaG_ao_iw(nfreqs,nBas_twice,nBas_twice))
 allocate(G_ao_itau_hfb(ntimes_twice,nBas_twice,nBas_twice))
 allocate(G_ao_itau(ntimes_twice,nBas_twice,nBas_twice))
 allocate(Chi0_ao_itau(nBasSq,nBasSq))
 allocate(Chi0_ao_iw(nfreqs,nBasSq,nBasSq),Wp_ao_iw(nBasSq,nBasSq))
 allocate(Wp_ao_itau(ntimes,nBasSq,nBasSq))
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

 ! Initialize arrays
 DeltaG_ao_iw(:,:,:)=czero
 Mat1(1:nOrb,1:nOrb)=U_QP(1:nOrb,1:nOrb)
 Mat2(1:nOrb,1:nOrb)=U_QP(nOrb+1:nOrb_twice,1:nOrb)
 R_ao(1:nBas           ,1:nBas           ) = 0.5d0*P_in(1:nBas,1:nBas)
 R_ao(nBas+1:nBas_twice,nBas+1:nBas_twice) = matmul(X_in(1:nBas,1:nOrb), transpose(X_in(1:nBas,1:nOrb)))-0.5d0*P_in(1:nBas,1:nBas)
 R_ao(1:nBas           ,nBas+1:nBas_twice) = Pan_in(1:nBas,1:nBas)
 R_ao(nBas+1:nBas_twice,1:nBas           ) = Pan_in(1:nBas,1:nBas)
 cHFBinv=matmul(transpose(cHFB),S)

 ! Read grids 
 call read_scGW_grids(ntimes,nfreqs,tcoord,tweight,wcoord,wweight,sint2w_weight,cost2w_weight, &
                      cosw2t_weight,sinw2t_weight,verbose)

 ! Build Go(i w)
 do ifreq=1,nfreqs
  weval_cpx=im*wcoord(ifreq)
  ! G_he(iw)
  call G_AO_RHFB(nBas,nOrb,nOrb_twice,eta,cHFB,eQP_state,weval_cpx, Mat1, Mat1, Mat2, Mat2,G_ao_tmp)
  G_ao_iw_hfb(ifreq,1:nBas           ,1:nBas           ) = G_ao_tmp(1:nBas,1:nBas)
  ! G_hh(iw)
  call G_AO_RHFB(nBas,nOrb,nOrb_twice,eta,cHFB,eQP_state,weval_cpx, Mat1, Mat2,-Mat2, Mat1,G_ao_tmp)
  G_ao_iw_hfb(ifreq,1:nBas           ,nBas+1:nBas_twice) = G_ao_tmp(1:nBas,1:nBas)
  ! G_ee(iw)
  call G_AO_RHFB(nBas,nOrb,nOrb_twice,eta,cHFB,eQP_state,weval_cpx, Mat2, Mat1, Mat1,-Mat2,G_ao_tmp)
  G_ao_iw_hfb(ifreq,nBas+1:nBas_twice,1:nBas           ) = G_ao_tmp(1:nBas,1:nBas)
  ! G_eh(iw)
  call G_AO_RHFB(nBas,nOrb,nOrb_twice,eta,cHFB,eQP_state,weval_cpx, Mat2, Mat2, Mat1, Mat1,G_ao_tmp)
  G_ao_iw_hfb(ifreq,nBas+1:nBas_twice,nBas+1:nBas_twice) = G_ao_tmp(1:nBas,1:nBas)
 enddo
 ! Build Go(i tau)
 do itau=1,ntimes
  ! tau > 0
  ! G_he(i tau)
  call G0itau_ao_RHFB(nBas,nOrb,nOrb_twice, tcoord(itau),G_ao_tmp,cHFB,eQP_state, Mat1, Mat1, Mat2, Mat2)
  G_ao_itau_hfb(2*itau-1,1:nBas           ,1:nBas           ) = G_ao_tmp(1:nBas,1:nBas)
  ! G_hh(i tau)
  call G0itau_ao_RHFB(nBas,nOrb,nOrb_twice, tcoord(itau),G_ao_tmp,cHFB,eQP_state, Mat1, Mat2,-Mat2, Mat1)
  G_ao_itau_hfb(2*itau-1,1:nBas           ,nBas+1:nBas_twice) = G_ao_tmp(1:nBas,1:nBas)
  ! G_ee(i tau)
  call G0itau_ao_RHFB(nBas,nOrb,nOrb_twice, tcoord(itau),G_ao_tmp,cHFB,eQP_state, Mat2, Mat1, Mat1,-Mat2)
  G_ao_itau_hfb(2*itau-1,nBas+1:nBas_twice,1:nBas           ) = G_ao_tmp(1:nBas,1:nBas)
  ! G_eh(i tau)
  call G0itau_ao_RHFB(nBas,nOrb,nOrb_twice, tcoord(itau),G_ao_tmp,cHFB,eQP_state, Mat2, Mat2, Mat1, Mat1)
  G_ao_itau_hfb(2*itau-1,nBas+1:nBas_twice,nBas+1:nBas_twice) = G_ao_tmp(1:nBas,1:nBas)
  ! tau < 0
  ! G_he(i tau)
  call G0itau_ao_RHFB(nBas,nOrb,nOrb_twice,-tcoord(itau),G_ao_tmp,cHFB,eQP_state, Mat1, Mat1, Mat2, Mat2)
  G_ao_itau_hfb(2*itau  ,1:nBas           ,1:nBas           ) = G_ao_tmp(1:nBas,1:nBas)
  ! G_hh(i tau)
  call G0itau_ao_RHFB(nBas,nOrb,nOrb_twice,-tcoord(itau),G_ao_tmp,cHFB,eQP_state, Mat1, Mat2,-Mat2, Mat1)
  G_ao_itau_hfb(2*itau  ,1:nBas           ,nBas+1:nBas_twice) = G_ao_tmp(1:nBas,1:nBas)
  ! G_ee(i tau)
  call G0itau_ao_RHFB(nBas,nOrb,nOrb_twice,-tcoord(itau),G_ao_tmp,cHFB,eQP_state, Mat2, Mat1, Mat1,-Mat2)
  G_ao_itau_hfb(2*itau  ,nBas+1:nBas_twice,1:nBas           ) = G_ao_tmp(1:nBas,1:nBas)
  ! G_eh(i tau)
  call G0itau_ao_RHFB(nBas,nOrb,nOrb_twice,-tcoord(itau),G_ao_tmp,cHFB,eQP_state, Mat2, Mat2, Mat1, Mat1)
  G_ao_itau_hfb(2*itau  ,nBas+1:nBas_twice,nBas+1:nBas_twice) = G_ao_tmp(1:nBas,1:nBas)
 enddo
 ! Initialize G(i tau)
 G_ao_itau(:,:,:)=G_ao_itau_hfb(:,:,:)
 ! Check error in the Fourier transformation
 write(*,*)
 write(*,'(a)') ' Error test for the Go(iw) -> G(it) transformation'
 write(*,*)
 max_error_gw2gt=-1d0
 sum_error_gw2gt=0d0
 imax_error_gw2gt=1
 do itau=1,ntimes
  ! tau > 0
  G_gorkov_tmp=czero
  do ifreq=1,nfreqs
   G_gorkov_tmp(:,:) = G_gorkov_tmp(:,:) + im*cosw2t_weight(itau,ifreq)*Real(G_ao_iw_hfb(ifreq,:,:))  &
                                         - im*sinw2t_weight(itau,ifreq)*Aimag(G_ao_iw_hfb(ifreq,:,:))
  enddo
  G_gorkov_tmp(:,:)=abs(G_gorkov_tmp(:,:)-G_ao_itau_hfb(2*itau-1,:,:))
  error_gw2gt=real(sum(G_gorkov_tmp(:,:)))
  ! tau < 0
  G_gorkov_tmp=czero
  do ifreq=1,nfreqs
   G_gorkov_tmp(:,:) = G_gorkov_tmp(:,:) + im*cosw2t_weight(itau,ifreq)*Real(G_ao_iw_hfb(ifreq,:,:))  &
                                         + im*sinw2t_weight(itau,ifreq)*Aimag(G_ao_iw_hfb(ifreq,:,:))
  enddo
  G_gorkov_tmp(:,:)=abs(G_gorkov_tmp(:,:)-G_ao_itau_hfb(2*itau  ,:,:))
  error_gw2gt=error_gw2gt+real(sum(G_gorkov_tmp(:,:)))
  sum_error_gw2gt=sum_error_gw2gt+error_gw2gt
  if(error_gw2gt>max_error_gw2gt) then
   imax_error_gw2gt=itau
   max_error_gw2gt=error_gw2gt
  endif
 enddo
 write(*,'(a,*(f20.8))') ' Sum error ',sum_error_gw2gt
 write(*,'(a,f20.8,a,2f20.8,a)') ' Max CAE   ',max_error_gw2gt,' is in the time +/-',0d0,tcoord(imax_error_gw2gt),'i'
 write(*,'(a,*(f20.8))') ' MAE       ',sum_error_gw2gt/(nfreqs*nBasSq*nBasSq)

 
  ! Build using the time grid Xo(i tau) = -2i [ G_he(i tau) G_he(-i tau) + G_hh(i tau) G_ee(-i tau) ]
  !  then Fourier transform Xo(i tau) -> Xo(i w)
  Chi0_ao_iw(:,:,:)=czero
  do itau=1,ntimes
   ! Xo(i tau) = -2i [ G_he(i tau) G_he(-i tau) + G_hh(i tau) G_ee(-i tau) ]
   do ibas=1,nBas
    do jbas=1,nBas
     mbas=nBas+1+(jbas-1)
     do kbas=1,nBas
      obas=nBas+1+(kbas-1)
      do lbas=1,nBas
                                   ! r1   r2'                    r2   r1'
       product = G_ao_itau(2*itau-1,ibas,jbas)*G_ao_itau(2*itau,kbas,lbas) &
               + G_ao_itau(2*itau-1,ibas,mbas)*G_ao_itau(2*itau,obas,lbas)  ! This is the right sign when G_hh and G_ee are built with  - Mat2
       if(abs(product)<1e-12) product=czero
       Chi0_ao_itau(1+(lbas-1)+(ibas-1)*nBas,1+(kbas-1)+(jbas-1)*nBas) = product
      enddo
     enddo
    enddo
   enddo
   Chi0_ao_itau=-2d0*im*Chi0_ao_itau ! The 2 factor is added to account for both spin contributions [ i.e., (up,up,up,up) and (down,down,down,down) ]
   ! Xo(i tau) -> Xo(i w)            [ the weight already contains the cos(tau w) and a factor 2 because int_-Infty ^Infty -> 2 int_0 ^Infty ]
   do ifreq=1,nfreqs
    Chi0_ao_iw(ifreq,:,:) = Chi0_ao_iw(ifreq,:,:) - im*cost2w_weight(ifreq,itau)*Chi0_ao_itau(:,:)
   enddo
  enddo
  ! Complete the Xo(i tau) -> Xo(i w)
  Chi0_ao_iw(:,:,:) = Real(Chi0_ao_iw(:,:,:)) ! The factor 2 is stored in the weight [ and we just retain the real part ]

  ! Build Wp(i w) and compute Ec Galitskii-Migdal 
  !  and Wp(i w) -> Wp(i tau)
  EcGM=0d0
  Wp_ao_itau=0d0
  do ifreq=1,nfreqs
   trace1=0d0; trace2=0d0;
   ! Xo(i w) -> Wp_ao_iw(i w)
   Wp_ao_iw(:,:)=-matmul(Real(Chi0_ao_iw(ifreq,:,:)),vMAT(:,:))
   do ibas=1,nBasSq
    trace1=trace1+Wp_ao_iw(ibas,ibas)
    Wp_ao_iw(ibas,ibas)=Wp_ao_iw(ibas,ibas)+1d0
   enddo
   call inverse_matrix(nBasSq,Wp_ao_iw,Wp_ao_iw)
   Wp_ao_iw(:,:)=matmul(Wp_ao_iw(:,:),Real(Chi0_ao_iw(ifreq,:,:)))
   Wp_ao_iw(:,:)=matmul(Wp_ao_iw(:,:),vMAT(:,:))
   do ibas=1,nBasSq
    trace2=trace2+Wp_ao_iw(ibas,ibas)
   enddo
   EcGM=EcGM-wweight(ifreq)*(trace2+trace1)/(2d0*pi) ! iw contribution to EcGM
   Wp_ao_iw(:,:)=matmul(vMAT(:,:),Wp_ao_iw(:,:))     ! Now Wp_ao_iw is on the iw grid
   ! Wp(i w) -> Wp(i tau) [ this transformation misses that Fourier[ Wp(i tau) ] is imaginary because of the factor i / 2pi ]
   !                      [ However, the weight contains a 2 /(2 pi) = 1 / pi factor and the cos(tau w).                    ]
   do itau=1,ntimes
    Wp_ao_itau(itau,:,:) = Wp_ao_itau(itau,:,:) + cosw2t_weight(itau,ifreq)*Wp_ao_iw(:,:)
   enddo
  enddo

  ! Build Sigma_c(i w) [Eqs. 12-18 in PRB, 109, 255101 (2024) ]
  Sigma_c_w_ao=czero
  do itau=1,ntimes
   Sigma_c_plus=czero
   Sigma_c_minus=czero
   ! Sigma_c(i tau) = i G(i tau) Wp(i tau)
   do ibas=1,nBas
    mbas=nBas+1+(ibas-1)
    do jbas=1,nBas
     obas=nBas+1+(jbas-1)
     do kbas=1,nBas
      pbas=nBas+1+(kbas-1)
      do lbas=1,nBas
       qbas=nBas+1+(lbas-1)
       ! Sigma_c_he(i tau) = i G_he Wp(i tau)
       Sigma_c_plus(ibas,jbas) =Sigma_c_plus(ibas,jbas)+im*G_ao_itau(2*itau-1,kbas,lbas)    &
                               *im*Wp_ao_itau(itau,1+(kbas-1)+(ibas-1)*nBas,1+(jbas-1)+(lbas-1)*nBas) ! Adding i to Wp that was missing
       Sigma_c_minus(ibas,jbas)=Sigma_c_minus(ibas,jbas)+im*G_ao_itau(2*itau ,kbas,lbas)  &
                               *im*Wp_ao_itau(itau,1+(kbas-1)+(ibas-1)*nBas,1+(jbas-1)+(lbas-1)*nBas) ! Adding i to Wp that was missing
       ! Sigma_c_hh(i tau) = - i G_hh Wp(i tau)
       Sigma_c_plus(ibas,obas) =Sigma_c_plus(ibas,obas)+im*G_ao_itau(2*itau-1,kbas,qbas)    &
                               *im*Wp_ao_itau(itau,1+(kbas-1)+(ibas-1)*nBas,1+(lbas-1)+(jbas-1)*nBas) ! Adding i to Wp that was missing
       Sigma_c_minus(ibas,obas)=Sigma_c_minus(ibas,obas)+im*G_ao_itau(2*itau ,kbas,qbas)  &
                               *im*Wp_ao_itau(itau,1+(kbas-1)+(ibas-1)*nBas,1+(lbas-1)+(jbas-1)*nBas) ! Adding i to Wp that was missing
       ! Sigma_c_ee(i tau) = - i G_ee Wp(i tau)
       Sigma_c_plus(mbas,jbas) =Sigma_c_plus(mbas,jbas)+im*G_ao_itau(2*itau-1,pbas,lbas)    &
                               *im*Wp_ao_itau(itau,1+(ibas-1)+(kbas-1)*nBas,1+(jbas-1)+(lbas-1)*nBas) ! Adding i to Wp that was missing
       Sigma_c_minus(mbas,jbas)=Sigma_c_minus(mbas,jbas)+im*G_ao_itau(2*itau ,pbas,lbas)  &
                               *im*Wp_ao_itau(itau,1+(ibas-1)+(kbas-1)*nBas,1+(jbas-1)+(lbas-1)*nBas) ! Adding i to Wp that was missing
       ! Sigma_c_eh(i tau) = i G_eh Wp(i tau)
       Sigma_c_plus(mbas,obas) =Sigma_c_plus(mbas,obas)+im*G_ao_itau(2*itau-1,pbas,qbas)    &
                               *im*Wp_ao_itau(itau,1+(ibas-1)+(kbas-1)*nBas,1+(lbas-1)+(jbas-1)*nBas) ! Adding i to Wp that was missing
       Sigma_c_minus(mbas,obas)=Sigma_c_minus(mbas,obas)+im*G_ao_itau(2*itau ,pbas,qbas)  &
                               *im*Wp_ao_itau(itau,1+(ibas-1)+(kbas-1)*nBas,1+(lbas-1)+(jbas-1)*nBas) ! Adding i to Wp that was missing
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



! Using the correlated G and Sigma_c to test the linearized density matrix approximation
 if(dolinGW) then
  write(*,*)
  write(*,*) ' -----------------------------------------------------'
  write(*,*) ' Testing the linearized approximation with G^Gorkov'
  write(*,*) '  G^lin,Gorkov = G^Gorkov + G^Gorkov Sigma_c G^Gorkov'
  write(*,*) ' -----------------------------------------------------'
  R_ao_old=0d0
  G_gorkov_tmp(:,:)=czero
  do ifreq=1,nfreqs
   G_gorkov_tmp(:,:)=G_ao_iw_hfb(ifreq,:,:)+DeltaG_ao_iw(ifreq,:,:)
   G_gorkov_tmp(:,:)=matmul(matmul(G_gorkov_tmp(:,:),Sigma_c_w_ao(ifreq,:,:)),G_gorkov_tmp(:,:))
   R_ao_old(:,:) = R_ao_old(:,:) + wweight(ifreq)*real(G_gorkov_tmp(:,:)+conjg(G_gorkov_tmp(:,:))) ! Integrate along iw
  enddo
  R_ao_old=R_ao_old/pi
  R_ao_old(1:nBas,1:nBas)=2d0*R_ao(1:nBas,1:nBas)+R_ao_old(1:nBas,1:nBas)       ! Sum both spin channels
  R_ao_old(1:nBas,nBas+1:)=R_ao(1:nBas,nBas+1:)+0.5d0*R_ao_old(1:nBas,nBas+1:)  ! We only need one spin-channel
  trace_1_rdm=0d0
  do ibas=1,nBas
   do jbas=1,nBas
    trace_1_rdm=trace_1_rdm+R_ao_old(ibas,jbas)*S(ibas,jbas)
   enddo
  enddo
  Ecore=0d0; Eh=0d0; Ex=0d0; Epair=0d0;
  do ibas=1,nBas
   do jbas=1,nBas
    obas=nBas+1+(jbas-1)
    Ecore=Ecore+R_ao_old(ibas,jbas)*Hc(ibas,jbas)
    do kbas=1,nBas
     do lbas=1,nBas
      qbas=nBas+1+(lbas-1)
      Eh=Eh+0.5d0*R_ao_old(kbas,lbas)*R_ao_old(ibas,jbas)*vMAT(1+(lbas-1)+(kbas-1)*nBas,1+(jbas-1)+(ibas-1)*nBas)
      Ex=Ex-0.25d0*R_ao_old(kbas,lbas)*R_ao_old(ibas,jbas)*vMAT(1+(jbas-1)+(kbas-1)*nBas,1+(lbas-1)+(ibas-1)*nBas)
      Epair=Epair+sigma*R_ao_old(ibas,obas)*R_ao_old(kbas,qbas)*vMAT(1+(ibas-1)+(kbas-1)*nBas,1+(jbas-1)+(lbas-1)*nBas) 
     enddo
    enddo
   enddo
  enddo
  Ehfbl=Ecore+Eh+Ex+Epair
  U_mo=-matmul(matmul(cHFBinv,R_ao_old(1:nBas,1:nBas)),transpose(cHFBinv)) ! Minus to order occ numbers
  call diagonalize_matrix(nOrb,U_mo,Occ)
  Occ=-Occ
  cNO=matmul(cHFB,U_mo)
  write(*,'(a,f15.8)')        ' Enuc        ',ENuc
  write(*,'(a,f15.8)')        ' Hcore       ',Ecore
  write(*,'(a,f15.8)')        ' Hartree     ',Eh
  write(*,'(a,f15.8)')        ' Exchange    ',Ex
  write(*,'(a,f15.8)')        ' Epairing    ',Epair
  write(*,'(a,f15.8)')        ' Ehfbl       ',Ehfbl
  write(*,'(a,f15.8)')        ' EcGM        ',EcGM
  write(*,'(a,f15.8)')        ' Eelec       ',Ehfbl+EcGM
  write(*,'(a,f15.8)')        ' lin-G Energy',Ehfbl+EcGM+ENuc
  write(*,*)
  write(*,'(a,f15.8,a,i5,a)') ' Trace lin-scGWB  ',trace_1_rdm
  write(*,*)
  write(*,*) ' Lin-G Bogoliubov occupation numbers'
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



 ! Deallocate arrays
 deallocate(Occ)
 deallocate(cHFBinv)
 deallocate(cNO)
 deallocate(U_mo)
 deallocate(R_ao)
 deallocate(R_ao_old)
 deallocate(G_ao_tmp)
 deallocate(G_gorkov_tmp)
 deallocate(DeltaG_ao_iw)
 deallocate(G_ao_iw_hfb)
 deallocate(G_ao_itau_hfb)
 deallocate(G_ao_itau)
 deallocate(Chi0_ao_itau)
 deallocate(Chi0_ao_iw,Wp_ao_iw,Wp_ao_itau)
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


 call wall_time(end_scGWBitauiw)
 
 t_scGWBitauiw = end_scGWBitauiw - start_scGWBitauiw
 write(*,'(A65,1X,F9.3,A8)') 'Total wall time for scGWB = ',t_scGWBitauiw,' seconds'
 write(*,*)

end subroutine 
