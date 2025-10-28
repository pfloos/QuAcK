subroutine scGWitauiw_ao(nBas,nOrb,nO,maxSCF,dolinGW,restart_scGW,no_fock,ENuc,Hc,S,P_in,cHF,eHF, &
                        nfreqs,wcoord,wweight,vMAT,ERI_AO)

! Restricted scGW

  implicit none
  include 'parameters.h'

! Input variables
 
  logical,intent(in)            :: dolinGW
  logical,intent(in)            :: no_fock
  logical,intent(in)            :: restart_scGW

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nO
  integer,intent(in)            :: maxSCF

  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: Hc(nBas,nBas)
  double precision,intent(in)   :: P_in(nBas,nBas)
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: cHF(nBas,nOrb)
  double precision,intent(in)   :: vMAT(nBas*nBas,nBas*nBas)
  double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas)

! Local variables
 
  logical                       :: file_exists

  integer                       :: iunit=312
  integer                       :: verbose
  integer                       :: ntimes
  integer                       :: ntimes_twice
  integer                       :: itau,ifreq
  integer                       :: ibas,jbas,kbas,lbas,nBas2
  integer                       :: iter,iter_fock
  integer                       :: imax_error_sigma

  double precision              :: start_scGWitauiw     ,end_scGWitauiw       ,t_scGWitauiw

  double precision              :: alpha_mixing
  double precision              :: val_print_r
  double precision              :: Ehfl,EcGM
  double precision              :: trace1,trace2
  double precision              :: eta,diff_Pao
  double precision              :: nElectrons
  double precision              :: trace_1_rdm
  double precision              :: thrs_N,thrs_Pao
  double precision              :: chem_pot,chem_pot_saved
  double precision              :: error_I
  double precision              :: error_sigma
  double precision              :: max_error_sigma
  double precision,allocatable  :: tweight(:),tcoord(:)
  double precision,allocatable  :: I_weight(:,:)
  double precision,allocatable  :: sint2w_weight(:,:)
  double precision,allocatable  :: cost2w_weight(:,:)
  double precision,allocatable  :: cosw2t_weight(:,:)
  double precision,allocatable  :: sinw2t_weight(:,:)
  double precision,allocatable  :: Occ(:)
  double precision,allocatable  :: Wp_ao_iw(:,:)
  double precision,allocatable  :: cHFinv(:,:)
  double precision,allocatable  :: F_ao(:,:)
  double precision,allocatable  :: P_ao(:,:)
  double precision,allocatable  :: P_ao_cte(:,:)
  double precision,allocatable  :: P_ao_old(:,:)
  double precision,allocatable  :: P_ao_iter(:,:)
  double precision,allocatable  :: P_mo(:,:)
  double precision,allocatable  :: Wp_ao_itau(:,:,:)

  complex*16                    :: product
  complex*16                    :: val_print_c
  complex*16                    :: weval_cpx
  complex*16,allocatable        :: Sigma_c_w_ao(:,:,:)
  complex*16,allocatable        :: DeltaG_ao_iw(:,:,:)
  complex*16,allocatable        :: G_ao_itau(:,:,:)
  complex*16,allocatable        :: G_ao_itau_old(:,:,:)
  complex*16,allocatable        :: G_ao_itau_cte(:,:,:)
  complex*16,allocatable        :: G_ao_iw_cte(:,:,:)
  complex*16,allocatable        :: Sigma_c_c(:,:),Sigma_c_s(:,:)
  complex*16,allocatable        :: Sigma_c_plus(:,:),Sigma_c_minus(:,:)
  complex*16,allocatable        :: G_ao_1(:,:),G_ao_2(:,:)
  complex*16,allocatable        :: G_minus_itau(:,:),G_plus_itau(:,:)
  complex*16,allocatable        :: Chi0_ao_itau(:,:)
  complex*16,allocatable        :: Chi0_ao_iw(:,:,:)
  complex*16,allocatable        :: error_transf_mo(:,:,:)
  complex*16,allocatable        :: Sigma_c_w_mo(:,:)

! Output variables
  integer,intent(inout)         :: nfreqs
  double precision,intent(inout):: eHF(nOrb)
  double precision,intent(inout):: wcoord(nfreqs)
  double precision,intent(inout):: wweight(nfreqs)
  
!------------------------------------------------------------------------
! Build G(i tau) in AO basis and use it to build Xo (i tau) -> Xo (i w)
!------------------------------------------------------------------------
 
 call wall_time(start_scGWitauiw)

 write(*,*)     
 write(*,*)'*******************************************'
 write(*,*)'*     scGW ( using itau and iw grids )    *'
 write(*,*)'*******************************************'
 write(*,*)

 verbose=0
 eta=0d0
 thrs_N=1d-8
 thrs_Pao=1d-6
 nElectrons=2d0*nO
 nBas2=nBas*nBas
 chem_pot_saved = 0.5d0*(eHF(nO)+eHF(nO+1))
 chem_pot = chem_pot_saved
 alpha_mixing=0.6d0
 Ehfl=0d0
 write(*,*)
 write(*,'(A33,1X,F16.10,A3)') ' Initial chemical potential  = ',chem_pot,' au'
 write(*,*)
 eHF(:) = eHF(:)-chem_pot_saved
   
 allocate(Chi0_ao_iw(nfreqs,nBas2,nBas2))
 allocate(P_ao(nBas,nBas),P_ao_old(nBas,nBas),P_ao_iter(nBas,nBas),P_ao_cte(nBas,nBas))
 allocate(F_ao(nBas,nBas),P_mo(nOrb,nOrb),cHFinv(nOrb,nBas),Occ(nOrb))
 allocate(G_minus_itau(nBas,nBas),G_plus_itau(nBas,nBas)) 
 allocate(G_ao_1(nBas,nBas),G_ao_2(nBas,nBas)) 
 allocate(Sigma_c_c(nBas,nBas),Sigma_c_s(nBas,nBas)) 
 allocate(Sigma_c_plus(nBas,nBas),Sigma_c_minus(nBas,nBas)) 
 allocate(Chi0_ao_itau(nBas2,nBas2),Wp_ao_iw(nBas2,nBas2)) 
 cHFinv=matmul(transpose(cHF),S)
 P_ao_cte=P_in
 P_ao=P_in
 P_ao_iter=P_in
 F_ao=Hc
 Ehfl=0d0
 trace_1_rdm=0d0
 do ibas=1,nBas
  do jbas=1,nBas
   Ehfl=Ehfl+P_ao(ibas,jbas)*Hc(ibas,jbas)
   trace_1_rdm=trace_1_rdm+P_ao(ibas,jbas)*S(ibas,jbas)
   do kbas=1,nBas
    do lbas=1,nBas
     F_ao(ibas,jbas)=F_ao(ibas,jbas)+P_ao(kbas,lbas)*vMAT(1+(lbas-1)+(kbas-1)*nBas,1+(jbas-1)+(ibas-1)*nBas) &
                    -0.5d0*P_ao(kbas,lbas)*vMAT(1+(jbas-1)+(kbas-1)*nBas,1+(lbas-1)+(ibas-1)*nBas)
     Ehfl=Ehfl+0.5d0*P_ao(kbas,lbas)*P_ao(ibas,jbas)*vMAT(1+(lbas-1)+(kbas-1)*nBas,1+(jbas-1)+(ibas-1)*nBas) &
         -0.25d0*P_ao(kbas,lbas)*P_ao(ibas,jbas)*vMAT(1+(jbas-1)+(kbas-1)*nBas,1+(lbas-1)+(ibas-1)*nBas)
    enddo
   enddo
  enddo
 enddo

!---------------!
! Reading grids !
!---------------!

 ntimes=nfreqs
 ntimes_twice=2*ntimes
 allocate(tweight(ntimes),tcoord(ntimes))
 allocate(I_weight(nfreqs,ntimes))
 allocate(sint2w_weight(nfreqs,ntimes))
 allocate(cost2w_weight(nfreqs,ntimes))
 allocate(cosw2t_weight(ntimes,nfreqs))
 allocate(sinw2t_weight(ntimes,nfreqs))
 allocate(Sigma_c_w_ao(nfreqs,nBas,nBas),DeltaG_ao_iw(nfreqs,nBas,nBas),G_ao_iw_cte(nfreqs,nBas,nBas))
 allocate(G_ao_itau(ntimes_twice,nBas,nBas),G_ao_itau_cte(ntimes_twice,nBas,nBas))
 allocate(G_ao_itau_old(ntimes_twice,nBas,nBas))
 allocate(Wp_ao_itau(ntimes,nBas2,nBas2))
 write(*,*)
 inquire(file='./grids/tcoord.txt', exist=file_exists)
 if(file_exists) then
  write(*,*) 'Reading tcoord from tcoord.txt'
  tcoord=0d0
  open(unit=iunit, form='formatted', file='./grids/tcoord.txt', status='old')
  read(iunit,*) ntimes
  do itau=1,ntimes
   read(iunit,*) tcoord(itau)
  enddo
  close(iunit)
 else
  write(*,*) ' Error! Could not find the tcoord.txt file'
  return
 endif
 inquire(file='./grids/tweight.txt', exist=file_exists)
 if(file_exists) then
  write(*,*) 'Reading tweight from tweight.txt'
  tweight=0d0
  open(unit=iunit, form='formatted', file='./grids/tweight.txt', status='old')
  read(iunit,*) ntimes
  do itau=1,ntimes
   read(iunit,*) tweight(itau)
  enddo
  close(iunit)
 else
  write(*,*) ' Error! Could not find the tweight.txt file'
  return
 endif
 inquire(file='./grids/wcoord.txt', exist=file_exists)
 if(file_exists) then
  write(*,*) 'Reading wcoord from wcoord.txt'
  wcoord=0d0
  open(unit=iunit, form='formatted', file='./grids/wcoord.txt', status='old')
  read(iunit,*) nfreqs
  do ifreq=1,nfreqs
   read(iunit,*) wcoord(ifreq)
  enddo
  close(iunit)
 else
  write(*,*) ' Error! Could not find the wcoord.txt file'
  return
 endif
 inquire(file='./grids/wweight.txt', exist=file_exists)
 if(file_exists) then
  write(*,*) 'Reading wweight from wweight.txt'
  wweight=0d0
  open(unit=iunit, form='formatted', file='./grids/wweight.txt', status='old')
  read(iunit,*) nfreqs
  do ifreq=1,nfreqs
   read(iunit,*) wweight(ifreq)
  enddo
  close(iunit)
 else
  write(*,*) ' Error! Could not find the wweight.txt file'
  return
 endif
 inquire(file='./grids/cost2w_weight.txt', exist=file_exists)
 if(file_exists) then
  write(*,*) 'Reading cost2w_weight from cost2w_weight.txt'
  cost2w_weight=0d0
  open(unit=iunit, form='formatted', file='./grids/cost2w_weight.txt', status='old')
  read(iunit,*) ifreq
  do ifreq=1,nfreqs
   do itau=1,ntimes
    read(iunit,*) cost2w_weight(ifreq,itau)
   enddo
  enddo
  close(iunit)
 else
  write(*,*) ' Error! Could not find the cost2w_weight.txt file'
  return
 endif
 inquire(file='./grids/cosw2t_weight.txt', exist=file_exists)
 if(file_exists) then
  write(*,*) 'Reading cosw2t_weight from cosw2t_weight.txt'
  cosw2t_weight=0d0
  open(unit=iunit, form='formatted', file='./grids/cosw2t_weight.txt', status='old')
  read(iunit,*) ifreq
  do itau=1,ntimes
   do ifreq=1,nfreqs
    read(iunit,*) cosw2t_weight(itau,ifreq)
   enddo
  enddo
  close(iunit)
 else
  write(*,*) ' Error! Could not find the cosw2t_weight.txt file'
  return
 endif
 inquire(file='./grids/sint2w_weight.txt', exist=file_exists)
 if(file_exists) then
  write(*,*) 'Reading sint2w_weight from sint2w_weight.txt'
  sint2w_weight=0d0
  open(unit=iunit, form='formatted', file='./grids/sint2w_weight.txt', status='old')
  read(iunit,*) ifreq
  do ifreq=1,nfreqs
   do itau=1,ntimes
    read(iunit,*) sint2w_weight(ifreq,itau)
   enddo
  enddo
  close(iunit)
 else
  write(*,*) ' Error! Could not find the sint2w_weight.txt file'
  return
 endif
 inquire(file='./grids/sinw2t_weight.txt', exist=file_exists)
 if(file_exists) then
  write(*,*) 'Reading sinw2t_weight from sinw2t_weight.txt'
  sinw2t_weight=0d0
  open(unit=iunit, form='formatted', file='./grids/sinw2t_weight.txt', status='old')
  read(iunit,*) ifreq
  do ifreq=1,nfreqs
   do itau=1,ntimes
    read(iunit,*) sinw2t_weight(ifreq,itau)
   enddo
  enddo
  close(iunit)
 else
  write(*,*) ' Error! Could not find the sinw2t_weight.txt file'
  return
 endif
 write(*,'(a,i5,a)') ' Using ',nfreqs,' frequencies and times'
 I_weight=matmul(cost2w_weight,cosw2t_weight)
 do ifreq=1,nfreqs
  I_weight(ifreq,ifreq)=I_weight(ifreq,ifreq)-1d0
 enddo
 error_I=0d0
 do ifreq=1,nfreqs
  do itau=1,ntimes
   error_I=error_I+abs(I_weight(ifreq,itau))
  enddo
 enddo
 write(*,'(a,f20.5)') ' Deviation from Identity in Cos-Cos ',error_I
 I_weight=matmul(sint2w_weight,sinw2t_weight)
 do ifreq=1,nfreqs
  I_weight(ifreq,ifreq)=I_weight(ifreq,ifreq)-1d0
 enddo
 error_I=0d0
 do ifreq=1,nfreqs
  do itau=1,ntimes
   error_I=error_I+abs(I_weight(ifreq,itau))
  enddo
 enddo
 write(*,'(a,f20.5)') ' Deviation from Identity in Sin-Sin ',error_I
 write(*,*)
 deallocate(I_weight)

!-----------!
! scGW loop !
!-----------!

 iter=0
 iter_fock=0
 do
  iter=iter+1

  ! For iter=1 we build G_ao_itau as the RHF one or read it from restart files
  ! [ we also initialize G_ao_iw_cte, G_ao_itau_cte, G_ao_itau_old, and (P_ao,P_ao_iter) ]
  if(iter==1) then
   G_ao_itau=czero
   do itau=1,ntimes
    call G0itau_ao_RHF(nBas,nOrb,nO, tcoord(itau),G_plus_itau ,cHF,eHF)
    call G0itau_ao_RHF(nBas,nOrb,nO,-tcoord(itau),G_minus_itau,cHF,eHF)
    G_ao_itau(2*itau-1,:,:)=G_plus_itau(:,:)
    G_ao_itau(2*itau  ,:,:)=G_minus_itau(:,:)
   enddo
   G_ao_itau_old(:,:,:)=G_ao_itau(:,:,:)
   G_ao_itau_cte(:,:,:)=G_ao_itau(:,:,:)
   do ifreq=1,nfreqs
    weval_cpx=im*wcoord(ifreq)
    call G_AO_RHF(nBas,nOrb,nO,eta,cHF,eHF,weval_cpx,G_ao_1)
    G_ao_iw_cte(ifreq,:,:)=G_ao_1(:,:)
   enddo
   ! Initialize DeltaG(i w) [ that will be G(i w) - Go(i w)/Gcte(i w) ]
   DeltaG_ao_iw(:,:,:)=czero
   ! If required, read the restart files
   if(restart_scGW) then
    write(*,*)
    write(*,'(a)') ' Reading restart files'
    write(*,*)
    open(unit=iunit,form='unformatted',file='scGW_Gitau_bin',status='old')
    write(*,'(a)') ' Reading scGW_Gitau_bin'
    read(iunit) ibas
    read(iunit) ibas
    do itau=1,ntimes_twice
     do ibas=1,nBas
      do jbas=1,nBas
       read(iunit) val_print_c
       G_ao_itau(itau,ibas,jbas)=val_print_c
      enddo
     enddo
    enddo
    close(iunit)
    open(unit=iunit,form='unformatted',file='scGW_Pao_bin',status='old')
    write(*,'(a)') ' Reading scGW_Pao_bin'
    read(iunit) ibas
    do ibas=1,nBas
     do jbas=1,nBas
      read(iunit) val_print_r
      P_ao(ibas,jbas)=val_print_r
     enddo
    enddo
    close(iunit)
    P_ao_iter=P_ao
    G_ao_itau_old(:,:,:)=G_ao_itau(:,:,:)
   endif
  endif

  ! Build using the time grid Xo(i tau) = -2i G(i tau) G(-i tau)
  !  then Fourier transform Xo(i tau) -> Xo(i w)
  Chi0_ao_iw(:,:,:)=czero
  do itau=1,ntimes
   ! Xo(i tau) = -2i G(i tau) G(-i tau)
   do ibas=1,nBas
    do jbas=1,nBas
     do kbas=1,nBas
      do lbas=1,nBas                       
                                   ! r1   r2'                    r2   r1'
       product = G_ao_itau(2*itau-1,ibas,jbas)*G_ao_itau(2*itau,kbas,lbas)
       if(abs(product)<1e-12) product=czero
       Chi0_ao_itau(1+(lbas-1)+(ibas-1)*nBas,1+(kbas-1)+(jbas-1)*nBas) = product
      enddo
     enddo
    enddo
   enddo
   Chi0_ao_itau=-2d0*im*Chi0_ao_itau ! The 2 factor is added to account for both spin contributions [ i.e., (up,up,up,up) and (down,down,down,down) ]
   ! Xo(i tau) -> Xo(i w) [ the weight already contains the cos(tau w) and a factor 2 because int_-Infty ^Infty -> 2 int_0 ^Infty ]
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
   do ibas=1,nBas2
    trace1=trace1+Wp_ao_iw(ibas,ibas)
    Wp_ao_iw(ibas,ibas)=Wp_ao_iw(ibas,ibas)+1d0
   enddo
   call inverse_matrix(nBas2,Wp_ao_iw,Wp_ao_iw)
   Wp_ao_iw(:,:)=matmul(Wp_ao_iw(:,:),Real(Chi0_ao_iw(ifreq,:,:)))
   Wp_ao_iw(:,:)=matmul(Wp_ao_iw(:,:),vMAT(:,:))
   do ibas=1,nBas2
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

  ! Build Sigma_c(i w)
  Sigma_c_w_ao=czero
  do itau=1,ntimes
   G_plus_itau(:,:) =G_ao_itau(2*itau-1,:,:)
   G_minus_itau(:,:)=G_ao_itau(2*itau  ,:,:)
   Sigma_c_plus=czero
   Sigma_c_minus=czero
   ! Sigma_c(i tau) = i G(i tau) Wp(i tau)
   do ibas=1,nBas
    do jbas=1,nBas
     do kbas=1,nBas
      do lbas=1,nBas 
       Sigma_c_plus(ibas,jbas) =Sigma_c_plus(ibas,jbas)+im*G_plus_itau(kbas,lbas)    &
                               *im*Wp_ao_itau(itau,1+(kbas-1)+(ibas-1)*nBas,1+(jbas-1)+(lbas-1)*nBas) ! Adding i to Wp that was missing
       Sigma_c_minus(ibas,jbas)=Sigma_c_minus(ibas,jbas)+im*G_minus_itau(kbas,lbas)  &
                               *im*Wp_ao_itau(itau,1+(kbas-1)+(ibas-1)*nBas,1+(jbas-1)+(lbas-1)*nBas) ! Adding i to Wp that was missing
      enddo
     enddo
    enddo
   enddo
   Sigma_c_c= -im*(Sigma_c_plus+Sigma_c_minus)
   Sigma_c_s= -   (Sigma_c_plus-Sigma_c_minus)
   ! Sigma_c(i tau) -> Sigma_c(i w)
   do ifreq=1,nfreqs
    Sigma_c_w_ao(ifreq,:,:) = Sigma_c_w_ao(ifreq,:,:)                        &
                            + 0.5d0*cost2w_weight(ifreq,itau)*Sigma_c_c(:,:) &
                            + 0.5d0*sint2w_weight(ifreq,itau)*Sigma_c_s(:,:)
   enddo 
  enddo

  ! Check the error in Sigma_c(i w) at iter=1 when G_...cte=Go [ i.e., it is the HF one ]
  if(iter==1 .and. .not.restart_scGW) then
   write(*,*)
   write(*,'(a)') ' Error test for the Sigma_c(iw) construction at iter 1 [ compared with the analytic Sigma_c(iw) obtained from HF ] '
   write(*,*)
   max_error_sigma=-1d0;imax_error_sigma=1;
   allocate(error_transf_mo(nfreqs,nOrb,nOrb),Sigma_c_w_mo(nOrb,nOrb))
   ! Build the analytic Sigma_c(iw)
   call build_analityc_rhf_Sigma_c_iw(nBas,nOrb,nO,verbose,cHF,eHF,nfreqs,wcoord,ERI_AO,error_transf_mo) ! error_transf_mo set to Sigma_c_mo(iw)
   do ifreq=1,nfreqs
    Sigma_c_w_mo=matmul(matmul(transpose(cHF(:,:)),Sigma_c_w_ao(ifreq,:,:)),cHF(:,:)) ! Fourier: Sigma_c_ao(iw) -> Sigma_c_wo(iw)
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
   write(*,'(a,f20.8,a,2f20.8,a)') ' Max error ',max_error_sigma,' in the frequency ',0d0,wcoord(imax_error_sigma),'i'
   write(*,'(a,*(f20.8))') ' MAE       ',sum(error_transf_mo)/(nfreqs*nBas*nBas)
   deallocate(error_transf_mo,Sigma_c_w_mo)
  endif

  ! Converge with respect to the Fock operator (using only good P_ao matrices)
  if(.not.no_fock) then ! Skiiping the opt w.r.t. the Fock operator we will do linearized approximation on Go -> [ lin-G = Go + Go Sigma Go ]
   iter_fock=0
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
        F_ao(ibas,jbas)=F_ao(ibas,jbas)+P_ao(kbas,lbas)*vMAT(1+(lbas-1)+(kbas-1)*nBas,1+(jbas-1)+(ibas-1)*nBas) &
                       -0.5d0*P_ao(kbas,lbas)*vMAT(1+(jbas-1)+(kbas-1)*nBas,1+(lbas-1)+(ibas-1)*nBas)
        Ehfl=Ehfl+0.5d0*P_ao(kbas,lbas)*P_ao(ibas,jbas)*vMAT(1+(lbas-1)+(kbas-1)*nBas,1+(jbas-1)+(ibas-1)*nBas) &
            -0.25d0*P_ao(kbas,lbas)*P_ao(ibas,jbas)*vMAT(1+(jbas-1)+(kbas-1)*nBas,1+(lbas-1)+(ibas-1)*nBas)
       enddo
      enddo
     enddo
    enddo
    ! Build G(i w) and n(r)
    P_ao_old=P_ao
    call get_1rdm_scGW(nBas,nfreqs,nElectrons,thrs_N,chem_pot,S,F_ao,Sigma_c_w_ao,wcoord,wweight, &
                       G_ao_1,G_ao_iw_cte,DeltaG_ao_iw,P_ao,P_ao_cte,trace_1_rdm) 
    if(abs(trace_1_rdm-nElectrons)>thrs_N) &
     call fix_chem_pot_scGW(nBas,nfreqs,nElectrons,thrs_N,chem_pot,S,F_ao,Sigma_c_w_ao,wcoord,wweight, &
                           G_ao_1,G_ao_iw_cte,DeltaG_ao_iw,P_ao,P_ao_cte,trace_1_rdm)
    ! Check convergence of P_ao for fixed Sigma_c(i w)
    diff_Pao=0d0
    do ibas=1,nBas
     do jbas=1,nBas
      diff_Pao=diff_Pao+abs(P_ao(ibas,jbas)-P_ao_old(ibas,jbas))
     enddo
    enddo

    if(diff_Pao<=thrs_Pao) exit
   
    if(iter_fock==maxSCF) exit
   
   enddo
  endif

  ! Check convergence of P_ao after a scGW iteration
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
  write(*,*)
  write(*,'(a,f15.8,a,i5,a,i5)') ' Trace scGW  ',trace_1_rdm,' after ',iter_fock,' Fock iterations at global iter ',iter
  write(*,'(a,f15.8)')        ' Change of P ',diff_Pao
  write(*,'(a,f15.8)')        ' Chem. Pot.  ',chem_pot
  write(*,'(a,f15.8)')        ' EcGM        ',EcGM
  write(*,'(a,f15.8)')        ' Eelec       ',Ehfl+EcGM
  write(*,'(a,f15.8)')        ' Etot        ',Ehfl+EcGM+ENuc
  write(*,*)
  write(*,*) ' Occupation numbers'
  Occ=-Occ
  do ibas=1,nOrb
   write(*,'(I7,F15.8)') ibas,Occ(ibas)
  enddo

  if(diff_Pao<=thrs_Pao) exit

  if(iter==maxSCF) exit

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
   ! Build G(i tau) = DeltaG(i tau) + Go(i tau)/Gcte(i tau)
   G_ao_itau(2*itau-1,:,:)=G_plus_itau(:,:) +G_ao_itau_cte(2*itau-1,:,:)
   G_ao_itau(2*itau  ,:,:)=G_minus_itau(:,:)+G_ao_itau_cte(2*itau  ,:,:)
  enddo
 
  ! Do mixing with previous G(i tau) to facilitate convergence
  G_ao_itau(:,:,:)=alpha_mixing*G_ao_itau(:,:,:)+(1d0-alpha_mixing)*G_ao_itau_old(:,:,:)
  G_ao_itau_old(:,:,:)=G_ao_itau(:,:,:)

 enddo
 write(*,*)
 write(*,'(A50)') '---------------------------------------'
 write(*,'(A50)') '      scGW calculation completed       '
 write(*,'(A50)') '---------------------------------------'
 write(*,*)
 write(*,'(a,f15.8,a,i5,a)') ' Trace scGW  ',trace_1_rdm,' after ',iter,' global iterations '
 write(*,'(a,f15.8)')        ' Change of P ',diff_Pao
 write(*,'(a,f15.8)')        ' Chem. Pot.  ',chem_pot
 write(*,'(a,f15.8)')        ' Hcore+Hx    ',Ehfl
 write(*,'(a,f15.8)')        ' EcGM        ',EcGM
 write(*,'(a,f15.8)')        ' Eelec       ',Ehfl+EcGM
 write(*,'(a,f15.8)')        ' scGW Energy ',Ehfl+EcGM+ENuc
 write(*,*)
 write(*,*) ' Final occupation numbers'
 do ibas=1,nOrb
  write(*,'(I7,F15.8)') ibas,Occ(ibas)
 enddo

  ! Write restart files
  open(unit=iunit,form='unformatted',file='scGW_Gitau_bin')
  write(iunit) nBas 
  write(iunit) ntimes
  do itau=1,ntimes_twice
   do ibas=1,nBas
    do jbas=1,nBas
     val_print_c = G_ao_itau(itau,ibas,jbas) 
     if(abs(val_print_c)<1d-8) val_print_c=czero
     write(iunit) val_print_c
    enddo
   enddo
  enddo
  write(iunit) iunit
  close(iunit)
  open(unit=iunit,form='unformatted',file='scGW_Giw_bin')
  write(iunit) nBas 
  write(iunit) ntimes
  do ifreq=1,nfreqs
   do ibas=1,nBas
    do jbas=1,nBas
     val_print_c = G_ao_iw_cte(ifreq,ibas,jbas)+DeltaG_ao_iw(ifreq,ibas,jbas) 
     if(abs(val_print_c)<1d-8) val_print_c=czero
     write(iunit) val_print_c
    enddo
   enddo
  enddo
  write(iunit) iunit
  close(iunit)
  open(unit=iunit,form='unformatted',file='scGW_Pao_bin')
  write(iunit) nBas 
  do ibas=1,nBas
   do jbas=1,nBas
    val_print_r = P_ao(ibas,jbas) 
    if(abs(val_print_r)<1d-8) val_print_r=czero
    write(iunit) val_print_r
   enddo
  enddo
  write(iunit) iunit
  close(iunit)

 ! Using the correlated G and Sigma_c to test the linearized density matrix approximation
 if(dolinGW) then
  write(*,*)
  write(*,*) ' -------------------------------------------'
  write(*,*) ' Testing the linearized approximation with G'
  write(*,*) '         G^lin = G + G Sigma_c G'
  write(*,*) ' -------------------------------------------'
  P_ao_old=0d0
  G_ao_1(:,:)=czero
  do ifreq=1,nfreqs
   G_ao_1(:,:)=G_ao_iw_cte(ifreq,:,:)+DeltaG_ao_iw(ifreq,:,:)
   G_ao_1(:,:)=matmul(matmul(G_ao_1(:,:),Sigma_c_w_ao(ifreq,:,:)),G_ao_1(:,:))
   P_ao_old(:,:) = P_ao_old(:,:) + wweight(ifreq)*real(G_ao_1(:,:)+conjg(G_ao_1(:,:))) ! Integrate along iw
  enddo
  P_ao_old=P_ao_old/pi
  P_ao_old=P_ao+P_ao_old
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
      Ehfl=Ehfl+0.5d0*P_ao_old(kbas,lbas)*P_ao_old(ibas,jbas)*vMAT(1+(lbas-1)+(kbas-1)*nBas,1+(jbas-1)+(ibas-1)*nBas) &
          -0.25d0*P_ao_old(kbas,lbas)*P_ao_old(ibas,jbas)*vMAT(1+(jbas-1)+(kbas-1)*nBas,1+(lbas-1)+(ibas-1)*nBas)
     enddo
    enddo
   enddo
  enddo
  P_mo=-matmul(matmul(cHFinv,P_ao_old),transpose(cHFinv)) ! Minus to order occ numbers
  call diagonalize_matrix(nOrb,P_mo,Occ)
  write(*,'(a,f15.8)')        ' Hcore+Hx    ',Ehfl
  write(*,'(a,f15.8)')        ' EcGM        ',EcGM
  write(*,'(a,f15.8)')        ' Eelec       ',Ehfl+EcGM
  write(*,'(a,f15.8)')        ' lin-G Energy',Ehfl+EcGM+ENuc
  write(*,*)
  write(*,'(a,f15.8,a,i5,a)') ' Trace lin-scGW  ',trace_1_rdm
  write(*,*)
  write(*,*) ' Lin-G occupation numbers'
  Occ=-Occ
  do ibas=1,nOrb
   write(*,'(I7,F15.8)') ibas,Occ(ibas)
  enddo
 endif

 call wall_time(end_scGWitauiw)
 
 t_scGWitauiw = end_scGWitauiw - start_scGWitauiw
 write(*,'(A65,1X,F9.3,A8)') 'Total wall time for scGW = ',t_scGWitauiw,' seconds'
 write(*,*)

 ! Restore values and deallocate dyn arrays
 eHF(:) = eHF(:)+chem_pot_saved
 deallocate(Wp_ao_itau)
 deallocate(Chi0_ao_iw,Wp_ao_iw)
 deallocate(tcoord,tweight) 
 deallocate(sint2w_weight)
 deallocate(cost2w_weight)
 deallocate(cosw2t_weight)
 deallocate(sinw2t_weight)
 deallocate(G_ao_itau_old)
 deallocate(G_ao_itau,G_ao_itau_cte)
 deallocate(Sigma_c_w_ao,DeltaG_ao_iw,G_ao_iw_cte)
 deallocate(P_ao,P_ao_old,P_ao_iter,P_ao_cte,F_ao,P_mo,cHFinv,Occ) 
 deallocate(Sigma_c_plus,Sigma_c_minus) 
 deallocate(Sigma_c_c,Sigma_c_s) 
 deallocate(G_minus_itau,G_plus_itau) 
 deallocate(G_ao_1,G_ao_2) 
 deallocate(Chi0_ao_itau) 

end subroutine 

subroutine fix_chem_pot_scGW(nBas,nfreqs,nElectrons,thrs_N,chem_pot,S,F_ao,Sigma_c_w_ao,wcoord,wweight, &
                             G_ao,G_ao_iw_cte,DeltaG_ao_iw,P_ao,P_ao_cte,trace_1_rdm) 

! Fix the chemical potential for scGW 

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nfreqs
  double precision,intent(in)   :: nElectrons
  double precision,intent(in)   :: thrs_N
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: F_ao(nBas,nBas)
  double precision,intent(in)   :: wcoord(nfreqs)
  double precision,intent(in)   :: wweight(nfreqs)
  double precision,intent(in)   :: P_ao_cte(nBas,nBas)
  complex*16,intent(in)         :: Sigma_c_w_ao(nfreqs,nBas,nBas)
  complex*16,intent(in)         :: G_ao_iw_cte(nfreqs,nBas,nBas)

! Local variables

  integer                       :: isteps
  double precision              :: thrs_closer
  double precision              :: delta_chem_pot
  double precision              :: chem_pot_change
  double precision              :: chem_pot_old
  double precision              :: grad_electrons
  double precision              :: trace_2up
  double precision              :: trace_up
  double precision              :: trace_down
  double precision              :: trace_2down
  double precision              :: trace_old

! Output variables

  double precision,intent(inout):: chem_pot
  double precision,intent(out)  :: trace_1_rdm
  double precision,intent(out)  :: P_ao(nBas,nBas)
  complex*16,intent(out)        :: G_ao(nBas,nBas)
  complex*16,intent(out)        :: DeltaG_ao_iw(nfreqs,nBas,nBas)

  !  Initialize 

  isteps = 0
  delta_chem_pot  = 2d-1
  thrs_closer     = 2d-1
  chem_pot_change = 0d0
  grad_electrons  = 1d0
  trace_1_rdm      = -1d0

  write(*,*)
  write(*,'(a)') ' Fixing the Tr[1D] at scGW '
  write(*,*)
  write(*,*)'------------------------------------------------------'
  write(*,'(1X,A1,1X,A15,1X,A1,1X,A15,1X,A1A15,2X,A1)') &
          '|','Tr[1D]','|','Chem. Pot.','|','Grad N','|'
  write(*,*)'------------------------------------------------------'

  ! First approach close the value with an error lower than 1

  trace_old = 1d2
  do while( abs(trace_old-nElectrons) > thrs_closer .and. isteps <= 100 )
   isteps = isteps + 1
   call get_1rdm_scGW(nBas,nfreqs,nElectrons,thrs_N,chem_pot,S,F_ao,Sigma_c_w_ao, &
                      wcoord,wweight,G_ao,G_ao_iw_cte,DeltaG_ao_iw,P_ao,P_ao_cte,trace_old) 
   call get_1rdm_scGW(nBas,nfreqs,nElectrons,thrs_N,chem_pot-delta_chem_pot,S,F_ao,Sigma_c_w_ao, &
                      wcoord,wweight,G_ao,G_ao_iw_cte,DeltaG_ao_iw,P_ao,P_ao_cte,trace_down) 
   call get_1rdm_scGW(nBas,nfreqs,nElectrons,thrs_N,chem_pot+delta_chem_pot,S,F_ao,Sigma_c_w_ao, &
                      wcoord,wweight,G_ao,G_ao_iw_cte,DeltaG_ao_iw,P_ao,P_ao_cte,trace_up) 
   if( abs(trace_up-nElectrons) > abs(trace_old-nElectrons) .and. abs(trace_down-nElectrons) > abs(trace_old-nElectrons) ) then
     write(*,'(1X,A1,F16.10,1X,A1,F16.10,1X,A1F16.10,1X,A1)') &
     '|',trace_old,'|',chem_pot,'|',grad_electrons,'|'
     delta_chem_pot = 0.5d0*delta_chem_pot
     thrs_closer = 0.5d0*thrs_closer
     write(*,*) "| contracting ...                                     |"
     if(delta_chem_pot<1d-2) exit
   else
     if( abs(trace_up-nElectrons) < abs(trace_old-nElectrons) ) then
      chem_pot=chem_pot+delta_chem_pot
      write(*,'(1X,A1,F16.10,1X,A1,F16.10,1X,A1F16.10,1X,A1)') &
      '|',trace_up,'|',chem_pot,'|',grad_electrons,'|'
     else
      if( abs(trace_down-nElectrons) < abs(trace_old-nElectrons) ) then
       chem_pot=chem_pot-delta_chem_pot
       write(*,'(1X,A1,F16.10,1X,A1,F16.10,1X,A1F16.10,1X,A1)') &
       '|',trace_down,'|',chem_pot,'|',grad_electrons,'|'
      endif
     endif
   endif
  enddo

  ! Do  final search

  write(*,*)'------------------------------------------------------'
  isteps = 0
  delta_chem_pot  = 1.0d-3
  do while( abs(trace_1_rdm-nElectrons) > thrs_N .and. isteps <= 1000 )
   isteps = isteps + 1
   chem_pot = chem_pot + chem_pot_change
   call get_1rdm_scGW(nBas,nfreqs,nElectrons,thrs_N,chem_pot,S,F_ao,Sigma_c_w_ao, &
                      wcoord,wweight,G_ao,G_ao_iw_cte,DeltaG_ao_iw,P_ao,P_ao_cte,trace_1_rdm)
   call get_1rdm_scGW(nBas,nfreqs,nElectrons,thrs_N,chem_pot+2d0*delta_chem_pot,S,F_ao,Sigma_c_w_ao, &
                      wcoord,wweight,G_ao,G_ao_iw_cte,DeltaG_ao_iw,P_ao,P_ao_cte,trace_2up)
   call get_1rdm_scGW(nBas,nfreqs,nElectrons,thrs_N,chem_pot+delta_chem_pot,S,F_ao,Sigma_c_w_ao, &
                      wcoord,wweight,G_ao,G_ao_iw_cte,DeltaG_ao_iw,P_ao,P_ao_cte,trace_up)
   call get_1rdm_scGW(nBas,nfreqs,nElectrons,thrs_N,chem_pot-delta_chem_pot,S,F_ao,Sigma_c_w_ao, &
                      wcoord,wweight,G_ao,G_ao_iw_cte,DeltaG_ao_iw,P_ao,P_ao_cte,trace_down)
   call get_1rdm_scGW(nBas,nfreqs,nElectrons,thrs_N,chem_pot-2d0*delta_chem_pot,S,F_ao,Sigma_c_w_ao, &
                      wcoord,wweight,G_ao,G_ao_iw_cte,DeltaG_ao_iw,P_ao,P_ao_cte,trace_2down)
!   grad_electrons = (trace_up-trace_down)/(2d0*delta_chem_pot)
   grad_electrons = (-trace_2up+8d0*trace_up-8d0*trace_down+trace_2down)/(12d0*delta_chem_pot)
   chem_pot_change = -(trace_1_rdm-nElectrons)/(grad_electrons+1d-10)
   ! Maximum change is bounded within +/- 0.10
   chem_pot_change = max( min( chem_pot_change , 0.1d0 / real(isteps) ), -0.1d0 / real(isteps) )
   write(*,'(1X,A1,F16.10,1X,A1,F16.10,1X,A1F16.10,1X,A1)') &
   '|',trace_1_rdm,'|',chem_pot,'|',grad_electrons,'|'
  enddo
  write(*,*)'------------------------------------------------------'
  write(*,*)
  call get_1rdm_scGW(nBas,nfreqs,nElectrons,thrs_N,chem_pot,S,F_ao,Sigma_c_w_ao, &
                     wcoord,wweight,G_ao,G_ao_iw_cte,DeltaG_ao_iw,P_ao,P_ao_cte,trace_old) 

end subroutine
    
subroutine get_1rdm_scGW(nBas,nfreqs,nElectrons,thrs_N,chem_pot,S,F_ao,Sigma_c_w_ao,wcoord,wweight, &
                         G_ao,G_ao_iw_cte,DeltaG_ao_iw,P_ao,P_ao_cte,trace_1_rdm) 

! Compute the scGW 1RDM

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nfreqs
  double precision,intent(in)   :: chem_pot
  double precision,intent(in)   :: nElectrons
  double precision,intent(in)   :: thrs_N
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: F_ao(nBas,nBas)
  double precision,intent(in)   :: wcoord(nfreqs)
  double precision,intent(in)   :: wweight(nfreqs)
  double precision,intent(in)   :: P_ao_cte(nBas,nBas)
  complex*16,intent(in)         :: G_ao_iw_cte(nfreqs,nBas,nBas)
  complex*16,intent(in)         :: Sigma_c_w_ao(nfreqs,nBas,nBas)

! Local variables

  integer                       :: ifreq
  integer                       :: ibas,jbas
  complex*16                    :: weval_cpx

! Output variables

  double precision,intent(out)  :: trace_1_rdm
  double precision,intent(out)  :: P_ao(nBas,nBas)
  complex*16,intent(out)        :: G_ao(nBas,nBas)
  complex*16,intent(out)        :: DeltaG_ao_iw(nfreqs,nBas,nBas)

  P_ao=0d0
  DeltaG_ao_iw=czero
  do ifreq=1,nfreqs
   weval_cpx=im*wcoord(ifreq)
   ! Setting G(w) = [ (w+chem_pot)S - F - Sigma_c(w) ]^-1
   G_ao(:,:)= (weval_cpx + chem_pot)*S(:,:) - F_ao(:,:) - Sigma_c_w_ao(ifreq,:,:) ! G(iw)^-1
   call complex_inverse_matrix(nBas,G_ao,G_ao)                                    ! G(iw)
   G_ao(:,:)=G_ao(:,:)-G_ao_iw_cte(ifreq,:,:)                                     ! G_corr(iw) = G(iw) - Go(iw)/Gcte(iw) 
   DeltaG_ao_iw(ifreq,:,:)=G_ao(:,:)
   P_ao(:,:) = P_ao(:,:) + wweight(ifreq)*real(G_ao(:,:))      ! P_corr = 1/(2 pi) int_-Infty ^Infty G_corr(iw) dw = 1/pi int_0 ^Infty Re[ G_corr(iw) ] dw
  enddo
  P_ao(:,:) = 2d0*P_ao(:,:)/pi + P_ao_cte(:,:)                 ! Times 2 to sum both spin channels
  trace_1_rdm=0d0
  do ibas=1,nBas
   do jbas=1,nBas
    trace_1_rdm=trace_1_rdm+P_ao(ibas,jbas)*S(ibas,jbas)
   enddo
  enddo

end subroutine

