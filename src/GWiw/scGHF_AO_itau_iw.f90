subroutine scGHF_AO_itau_iw(nBas,nOrb,nO,maxSCF,maxDIIS,verbose_scGHF,restart_scGHF,chem_pot_scG, &
                            ENuc,Hc,S,X,P_in,cHF,eHF,nfreqs,wcoord,wweight,vMAT)

! Restricted scGHF

  implicit none
  include 'parameters.h'

! Input variables
 
  logical,intent(in)            :: verbose_scGHF
  logical,intent(in)            :: restart_scGHF
  logical,intent(in)            :: chem_pot_scG

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nO
  integer,intent(in)            :: maxSCF
  integer,intent(in)            :: maxDIIS

  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: Hc(nBas,nBas)
  double precision,intent(in)   :: P_in(nBas,nBas)
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: X(nBas,nOrb)
  double precision,intent(in)   :: vMAT(nBas*nBas,nBas*nBas)

! Local variables
 
  logical                       :: file_exists
  logical                       :: read_SD_chkp

  integer                       :: n_diis
  integer                       :: n_diisP
  integer                       :: verbose
  integer                       :: ntimes
  integer                       :: nBasSqnfreqs
  integer                       :: ntimes_twice
  integer                       :: idiis_index
  integer                       :: idiis_indexP
  integer                       :: ifreq
  integer                       :: ibas,jbas,kbas,lbas,nBasSq
  integer                       :: iter,iter_fock
  integer                       :: imax_error_sigma
  integer                       :: imax_error_gw2gt

  double precision              :: rcond
  double precision              :: rcondP
  double precision              :: alpha_mixing
  double precision              :: Ehfl
  double precision              :: eta,diff_Pao
  double precision              :: nElectrons
  double precision              :: err_EcGM
  double precision              :: trace_1_rdm
  double precision              :: thrs_N,thrs_Ngrad,thrs_Pao
  double precision              :: chem_pot,chem_pot_saved,chem_pot_align
  double precision              :: error_sigma
  double precision              :: max_error_sigma
  double precision              :: error_gw2gt
  double precision              :: max_error_gw2gt
  double precision              :: sum_error_gw2gt
  double precision              :: start_scGHFitauiw,end_scGHFitauiw,t_scGHFitauiw
  double precision,allocatable  :: tweight(:),tcoord(:)
  double precision,allocatable  :: sint2w_weight(:,:)
  double precision,allocatable  :: cost2w_weight(:,:)
  double precision,allocatable  :: cosw2t_weight(:,:)
  double precision,allocatable  :: sinw2t_weight(:,:)
  double precision,allocatable  :: Occ(:)
  double precision,allocatable  :: cHFinv(:,:)
  double precision,allocatable  :: cNO(:,:)
  double precision,allocatable  :: F_ao(:,:)
  double precision,allocatable  :: U_mo(:,:)
  double precision,allocatable  :: P_ao(:,:)
  double precision,allocatable  :: P_ao_hf(:,:)
  double precision,allocatable  :: P_ao_old(:,:)
  double precision,allocatable  :: P_ao_iter(:,:)
  double precision,allocatable  :: err_currentP(:)
  double precision,allocatable  :: err_diisP(:,:)
  double precision,allocatable  :: P_ao_extrap(:)
  double precision,allocatable  :: P_ao_old_diis(:,:)

  complex*16                    :: product
  complex*16                    :: weval_cpx
  complex*16,allocatable        :: Sigma_c_w_ao(:,:,:)
  complex*16,allocatable        :: DeltaG_ao_iw(:,:,:)
  complex*16,allocatable        :: G_ao_itau(:,:,:)
  complex*16,allocatable        :: G_ao_iw_hf(:,:,:)
  complex*16,allocatable        :: G_ao_iw(:,:,:)
  complex*16,allocatable        :: G_ao_iw_old(:,:,:)
  complex*16,allocatable        :: Mat_ao_tmp(:,:)
  complex*16,allocatable        :: err_current(:)
  complex*16,allocatable        :: G_iw_extrap(:)
  complex*16,allocatable        :: err_diis(:,:)
  complex*16,allocatable        :: G_iw_old_diis(:,:)

! Output variables
  integer,intent(inout)         :: nfreqs
  double precision,intent(inout):: eHF(nOrb)
  double precision,intent(inout):: wcoord(nfreqs)
  double precision,intent(inout):: wweight(nfreqs)
  double precision,intent(inout):: cHF(nBas,nOrb)
  
!------------------------------------------------------------------------
! Build G(i tau) in AO basis and use it to build Xo (i tau) -> Xo (i w) !
!------------------------------------------------------------------------
 
 call wall_time(start_scGHFitauiw)

 write(*,*)     
 write(*,*)'*****************************************'
 write(*,*)'*    scGHF ( using it and iw grids )    *'
 write(*,*)'*****************************************'
 write(*,*)

 ! Initialize variables
 n_diis=0
 verbose=0
 if(verbose_scGHF) verbose=1
 eta=0d0
 thrs_N=1d-10
 thrs_Ngrad=1d-6
 thrs_Pao=1d-6
 nElectrons=2d0*nO
 nBasSq=nBas*nBas
 chem_pot_saved = 0.5d0*(eHF(nO)+eHF(nO+1))
 chem_pot = chem_pot_saved
 alpha_mixing=0.6d0
 rcond=0d0
 Ehfl=0d0
 ntimes=nfreqs
 ntimes_twice=2*ntimes
 nBasSqnfreqs=nBasSq*nfreqs
 read_SD_chkp=.false.
 write(*,*)
 write(*,'(A33,1X,F16.10,A3)') ' Initial chemical potential  = ',chem_pot,' au'
 if(chem_pot_scG) then
  write(*,'(A)') '   Adjusting the chemical potential is activated'
 else
  write(*,'(A)') '   Adjusting the chemical potential is deactivated'
 endif
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

 ! Allocate arrays  
 allocate(U_mo(nOrb,nOrb))
 allocate(P_ao(nBas,nBas),P_ao_old(nBas,nBas),P_ao_iter(nBas,nBas),P_ao_hf(nBas,nBas))
 allocate(F_ao(nBas,nBas),cHFinv(nOrb,nBas),Occ(nOrb),cNO(nBas,nOrb))
 allocate(Mat_ao_tmp(nBas,nBas)) 
 allocate(tweight(ntimes),tcoord(ntimes))
 allocate(sint2w_weight(nfreqs,ntimes))
 allocate(cost2w_weight(nfreqs,ntimes))
 allocate(cosw2t_weight(ntimes,nfreqs))
 allocate(sinw2t_weight(ntimes,nfreqs))
 allocate(G_ao_itau(ntimes_twice,nBas,nBas))
 allocate(G_ao_iw(nfreqs,nBas,nBas),G_ao_iw_old(nfreqs,nBas,nBas))
 allocate(Sigma_c_w_ao(nfreqs,nBas,nBas),DeltaG_ao_iw(nfreqs,nBas,nBas),G_ao_iw_hf(nfreqs,nBas,nBas))
 allocate(err_current(1))
 allocate(err_currentP(1))
 allocate(G_iw_extrap(1))
 allocate(P_ao_extrap(1))
 allocate(err_diis(1,1))
 allocate(err_diisP(1,1))
 allocate(G_iw_old_diis(1,1))
 allocate(P_ao_old_diis(1,1))
 if(maxDIIS>0) then
  deallocate(err_current)
  deallocate(err_currentP)
  deallocate(G_iw_extrap)
  deallocate(P_ao_extrap)
  deallocate(err_diis)
  deallocate(err_diisP)
  deallocate(G_iw_old_diis)
  deallocate(P_ao_old_diis)
  allocate(err_current(nBasSqnfreqs))
  allocate(err_currentP(nBasSq))
  allocate(G_iw_extrap(nBasSqnfreqs))
  allocate(P_ao_extrap(nBasSq))
  allocate(err_diis(nBasSqnfreqs,maxDIIS))
  allocate(err_diisP(nBasSq,maxDIIS))
  allocate(G_iw_old_diis(nBasSqnfreqs,maxDIIS))
  allocate(P_ao_old_diis(nBasSq,maxDIIS))
 endif

 ! Initialize arrays
 G_ao_itau=czero
 Sigma_c_w_ao=czero
 DeltaG_ao_iw=czero  ! Initialize DeltaG(i w) [ it will be G(i w) - Go(i w) ]  
 G_ao_iw=czero
 err_diis=czero
 G_iw_old_diis=czero
 cHFinv=matmul(transpose(cHF),S)
 P_ao_hf=P_in
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

 ! Read grids 
 call read_scGX_grids(ntimes,nfreqs,tcoord,tweight,wcoord,wweight,sint2w_weight,cost2w_weight, &
                      cosw2t_weight,sinw2t_weight,verbose)

 ! Build Go(i w)
 do ifreq=1,nfreqs
  weval_cpx=im*wcoord(ifreq)
  call G_AO_RHF_w(nBas,nOrb,nO,eta,cHF,eHF,weval_cpx,Mat_ao_tmp)
  G_ao_iw_hf(ifreq,:,:)=Mat_ao_tmp(:,:)
 enddo

 ! If required, read the restart files
 if(restart_scGHF) then
  call read_scGX_restart(nBas,nfreqs,ntimes_twice,chem_pot,P_ao,P_ao_hf,G_ao_iw_hf,G_ao_itau,G_ao_itau,read_SD_chkp)
  P_ao_iter=P_ao
 endif

!------------!
! scGHF loop !
!------------!

 iter=0
 iter_fock=0
 do
  iter=iter+1

  ! Converge with respect to the Fock operator (using only good P_ao matrices -> Tr[P_ao S_ao]=Nelectrons )
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
   call get_1rdm_scGW(nBas,nfreqs,chem_pot,S,F_ao,Sigma_c_w_ao,wcoord,wweight, &
                      Mat_ao_tmp,G_ao_iw_hf,DeltaG_ao_iw,P_ao,P_ao_hf,trace_1_rdm)
!  I. Duchemin's method to compute P_ao from G(i w) 
!   call get_1rdm_scGW_v2(nBas,nOrb,nfreqs,chem_pot,S,X,F_ao,Sigma_c_w_ao,wcoord,wweight, &
!                         Mat_ao_tmp,G_ao_iw_hf,DeltaG_ao_iw,P_ao,trace_1_rdm) 
   if(abs(trace_1_rdm-nElectrons)**2d0>thrs_N .and. chem_pot_scG) &
    call fix_chem_pot_scGW_bisec(iter_fock,nBas,nfreqs,nElectrons,thrs_N,thrs_Ngrad,chem_pot,S,F_ao,Sigma_c_w_ao,wcoord,wweight, &
                                 Mat_ao_tmp,G_ao_iw_hf,DeltaG_ao_iw,P_ao,P_ao_hf,trace_1_rdm,chem_pot_saved,verbose_scGHF)
   ! Check convergence of P_ao
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
      if(abs(err_currentP(idiis_indexP))<1e-12) err_currentP(idiis_indexP)=0d0
      P_ao_extrap(idiis_indexP)=P_ao(ibas,jbas)
      idiis_indexP=idiis_indexP+1
     enddo
    enddo
    call DIIS_extrapolation(rcondP,nBasSq,nBasSq,n_diisP,err_diisP,P_ao_old_diis,err_currentP,P_ao_extrap)
    idiis_indexP=1
    do ibas=1,nBas
     do jbas=1,nBas
      P_ao(ibas,jbas)=P_ao_extrap(idiis_indexP) 
      idiis_indexP=idiis_indexP+1
     enddo
    enddo
   else
    P_ao(:,:)=alpha_mixing*P_ao(:,:)+(1d0-alpha_mixing)*P_ao_old(:,:)
   endif
  
  enddo

  ! Check convergence of P_ao after a scGHF iteration
  diff_Pao=0d0
  do ibas=1,nBas
   do jbas=1,nBas
    diff_Pao=diff_Pao+abs(P_ao(ibas,jbas)-P_ao_iter(ibas,jbas))
   enddo
  enddo
  P_ao_iter=P_ao

  ! Print iter info
  U_mo=-matmul(matmul(cHFinv,P_ao),transpose(cHFinv)) ! Minus to order occ numbers
  call diagonalize_matrix(nOrb,U_mo,Occ)
  Occ=-Occ
  cNO=matmul(cHF,U_mo)
  write(*,*)
  write(*,'(a,f15.8,a,i5,a,i5)') ' Trace scGHF  ',trace_1_rdm,' after ',iter_fock,' Fock iterations at global iter ',iter
  write(*,'(a,f15.8)')        ' Change of P  ',diff_Pao
  write(*,'(a,f15.8)')        ' Chem. Pot.   ',chem_pot
  write(*,'(a,f15.8)')        ' Enuc         ',ENuc
  write(*,'(a,f15.8)')        ' Ehfl         ',Ehfl
  write(*,'(a,f15.8)')        ' Eelec        ',Ehfl
  write(*,'(a,f15.8)')        ' Etot         ',Ehfl+ENuc
  write(*,*)

  if(diff_Pao<=thrs_Pao) exit

  if(iter==maxSCF) exit
 
  G_ao_iw=G_ao_iw_hf+DeltaG_ao_iw 
  ! Do mixing with previous G(i w) to facilitate convergence
  if(maxDIIS>0) then
   n_diis=min(n_diis+1,maxDIIS)
   err_current=czero
   idiis_index=1
   do ifreq=1,nfreqs
    do ibas=1,nBas
     do jbas=1,nBas
      err_current(idiis_index)=G_ao_iw(ifreq,ibas,jbas)-G_ao_iw_old(ifreq,ibas,jbas)
      if(abs(err_current(idiis_index))<1e-12) err_current(idiis_index)=czero
      G_iw_extrap(idiis_index)=G_ao_iw(ifreq,ibas,jbas)
      idiis_index=idiis_index+1
     enddo
    enddo
   enddo
   call complex_DIIS_extrapolation(rcond,nBasSqnfreqs,nBasSqnfreqs,n_diis,err_diis,G_iw_old_diis,err_current,G_iw_extrap)
   idiis_index=1
   do ifreq=1,nfreqs
    do ibas=1,nBas
     do jbas=1,nBas
      G_ao_iw(ifreq,ibas,jbas)=G_iw_extrap(idiis_index) 
      idiis_index=idiis_index+1
     enddo
    enddo
   enddo
  else
   G_ao_iw(:,:,:)=alpha_mixing*G_ao_iw(:,:,:)+(1d0-alpha_mixing)*G_ao_iw_old(:,:,:)
  endif
  G_ao_iw_old(:,:,:)=G_ao_iw(:,:,:)

 enddo
 write(*,*)
 write(*,'(A50)') '---------------------------------------'
 write(*,'(A50)') '     scGHF calculation completed       '
 write(*,'(A50)') '---------------------------------------'
 write(*,*)
 write(*,'(a,f15.8,a,i5,a)') ' Trace scGHF  ',trace_1_rdm,' after ',iter,' global iterations '
 write(*,'(a,f15.8)')        ' Change of P  ',diff_Pao
 write(*,'(a,f15.8)')        ' Chem. Pot.   ',chem_pot
 write(*,'(a,f15.8)')        ' Enuc         ',ENuc
 write(*,'(a,f15.8)')        ' Ehfl         ',Ehfl
 write(*,'(a,f15.8)')        ' Eelec        ',Ehfl
 write(*,'(a,f15.8)')        ' scGHF Energy ',Ehfl+ENuc
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
 call write_scGX_restart(nBas,ntimes,ntimes_twice,nfreqs,chem_pot,P_ao,P_ao_hf,G_ao_itau,G_ao_itau, &
                         G_ao_iw_hf,DeltaG_ao_iw)

 call wall_time(end_scGHFitauiw)
 
 t_scGHFitauiw = end_scGHFitauiw - start_scGHFitauiw
 write(*,'(A65,1X,F9.3,A8)') 'Total wall time for scGHF = ',t_scGHFitauiw,' seconds'
 write(*,*)

 ! Restore values and deallocate dyn arrays
 eHF(:) = eHF(:)+chem_pot_saved
 deallocate(tcoord,tweight) 
 deallocate(sint2w_weight)
 deallocate(cost2w_weight)
 deallocate(cosw2t_weight)
 deallocate(sinw2t_weight)
 deallocate(G_ao_itau)
 deallocate(G_ao_iw,G_ao_iw_old)
 deallocate(Sigma_c_w_ao,DeltaG_ao_iw,G_ao_iw_hf)
 deallocate(P_ao,P_ao_old,P_ao_iter,P_ao_hf,F_ao,U_mo,cHFinv,cNO,Occ) 
 deallocate(Mat_ao_tmp) 
 deallocate(err_current)
 deallocate(err_currentP)
 deallocate(G_iw_extrap)
 deallocate(P_ao_extrap)
 deallocate(err_diis)
 deallocate(err_diisP)
 deallocate(G_iw_old_diis)
 deallocate(P_ao_old_diis)

end subroutine 
