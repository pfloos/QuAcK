! NOTE The usual qsGW recipe might not be good for qsGWB...!! 
! In standard qsGW, the Fock operator incorporates a static contribution from Sigma_c that is otained evaluating Sigma_c
! at the orbital energies; then, projecting it on the basis (i.e., the orbitals). Also, including an ad-hoc Hermitization step.
! The orbital energies can be aligned with the chemical potential to make the orbital energies of the occupied states to be - and
! + for the virtual states. Then, Sigma_c is evaluated using the POSITIVE ENERGIES FOR THE VIRTUAL STATES! 
! However, in HFB the contributions added to the Fock operator from Sigma_c are EVALUATED USING ONLY NEGATIVE ENERGIES energies if we
! follow the usual recipe. For this reason, even when HFB recovers HF, we do not have the proper contribution from Sigma_c to make
! qsGWB -> qsGW

subroutine sigc_MO_basis_RHFB(nOrb,nOrb_twice,offdiag0,eta,shift,Occ_el,U_QP,eqsGWB_state, &
                              vMAT,nfreqs,ntimes,wcoord,wweight,Sigc_mo_he,Sigc_mo_hh)

! Compute Sigma_c matrix in the AO basis

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: offdiag0

  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nOrb_twice
  integer,intent(in)            :: nfreqs
  integer,intent(in)            :: ntimes

  double precision,intent(in)   :: eta
  double precision,intent(in)   :: shift
  double precision,intent(in)   :: wcoord(nfreqs),wweight(nfreqs)
  double precision,intent(in)   :: Occ_el(nOrb)
  double precision,intent(in)   :: U_QP(nOrb_twice,nOrb_twice)
  double precision,intent(in)   :: eqsGWB_state(nOrb_twice)
  double precision,intent(in)   :: vMAT(nOrb*nOrb,nOrb*nOrb)

! Local variables

  logical                       :: doSigc_eh
  logical                       :: doSigc_ee

  integer                       :: iorb,jorb
  integer                       :: nE_eval_global

  double precision,allocatable  :: Sigc_mo_tmp(:,:,:)
  double precision,allocatable  :: Sigc_mo_tmp2(:,:,:)

  complex *16,allocatable       :: Sigc_mo_he_cpx(:,:,:)
  complex *16,allocatable       :: Sigc_mo_hh_cpx(:,:,:)
  complex *16,allocatable       :: Sigc_mo_eh_cpx(:,:,:)
  complex *16,allocatable       :: Sigc_mo_ee_cpx(:,:,:)
  complex *16,allocatable       :: E_eval_global_cpx(:)

! Output variables

  double precision,intent(out)  :: Sigc_mo_he(nOrb,nOrb)
  double precision,intent(out)  :: Sigc_mo_hh(nOrb,nOrb)

! Initialize

  doSigc_eh=.true.
  doSigc_ee=.true.

! Set energies using cluster method or just using two shifts
  if(.false.) then

   allocate(E_eval_global_cpx(1))
   call set_Eeval_cluster(nOrb,nOrb_twice,1,shift,eqsGWB_state,nE_eval_global,&
                          E_eval_global_cpx,0d0)
   deallocate(E_eval_global_cpx)
   allocate(E_eval_global_cpx(nE_eval_global))
   call set_Eeval_cluster(nOrb,nOrb_twice,nE_eval_global,shift,eqsGWB_state,&
                           nE_eval_global,E_eval_global_cpx,0d0)

   ! Run over unique energies
   allocate(Sigc_mo_he_cpx(nE_eval_global,nOrb,nOrb),Sigc_mo_hh_cpx(nE_eval_global,nOrb,nOrb))
   allocate(Sigc_mo_eh_cpx(nE_eval_global,nOrb,nOrb),Sigc_mo_ee_cpx(nE_eval_global,nOrb,nOrb))
   call build_Sigmac_w_RHFB(nOrb,nOrb_twice,nE_eval_global,eta,0,E_eval_global_cpx,eqsGWB_state, &
                           nfreqs,ntimes,wweight,wcoord,vMAT,U_QP,Sigc_mo_he_cpx,Sigc_mo_hh_cpx, &
                           Sigc_mo_eh_cpx,Sigc_mo_ee_cpx,doSigc_eh,doSigc_ee)
   deallocate(E_eval_global_cpx)

   ! TODO
   Sigc_mo_he=0d0
   Sigc_mo_hh=0d0

  else

   nE_eval_global=2*nOrb+1
   allocate(E_eval_global_cpx(nE_eval_global))
   do iorb=1,nOrb
    E_eval_global_cpx(2*iorb-1)=eqsGWB_state(iorb)-shift
    E_eval_global_cpx(2*iorb)  =eqsGWB_state(iorb)+shift
   enddo
   E_eval_global_cpx(nE_eval_global)=czero

   ! Build Sigma_c for all energies
   allocate(Sigc_mo_he_cpx(nE_eval_global,nOrb,nOrb),Sigc_mo_hh_cpx(nE_eval_global,nOrb,nOrb))
   allocate(Sigc_mo_eh_cpx(nE_eval_global,nOrb,nOrb),Sigc_mo_ee_cpx(nE_eval_global,nOrb,nOrb))
   call build_Sigmac_w_RHFB(nOrb,nOrb_twice,nE_eval_global,eta,0,E_eval_global_cpx,eqsGWB_state, &
                           nfreqs,ntimes,wweight,wcoord,vMAT,U_QP,Sigc_mo_he_cpx,Sigc_mo_hh_cpx, &
                           Sigc_mo_eh_cpx,Sigc_mo_ee_cpx,doSigc_eh,doSigc_ee)
   
   ! Interpolate and transform Sigma from MO to AO basis [incl. the usual qsGW recipe] 
   allocate(Sigc_mo_tmp(nOrb,nOrb,nOrb))
   allocate(Sigc_mo_tmp2(nOrb,nOrb,nOrb))
   Sigc_mo_tmp=0d0;Sigc_mo_tmp2=0d0;
   ! Sigma_c_he
   do iorb=1,nOrb ! Interpolation
    Sigc_mo_tmp(iorb,:,:) = 0.5d0*(Real(Sigc_mo_he_cpx(2*iorb-1,:,:))+Real(Sigc_mo_he_cpx(2*iorb,:,:)))
    Sigc_mo_tmp2(iorb,:,:)=-0.5d0*(Real(Sigc_mo_eh_cpx(2*iorb-1,:,:))+Real(Sigc_mo_eh_cpx(2*iorb,:,:)))
   enddo
   do iorb=1,nOrb
    Sigc_mo_he(iorb,:)=Occ_el(iorb)*Sigc_mo_tmp(iorb,iorb,:)+(1d0-Occ_el(iorb))*Sigc_mo_tmp2(iorb,iorb,:)
   enddo
   if(offdiag0) then  ! qsGW version where all the off-diagonal elements are built at the Fermi level 
    do iorb=1,nOrb
     do jorb=1,nOrb
      if(iorb/=jorb) then
       Sigc_mo_he(iorb,jorb) = Real(Sigc_mo_he_cpx(nE_eval_global,iorb,jorb))
      endif
     enddo
    enddo
   endif
   ! Sigma_c_hh
   Sigc_mo_tmp=0d0;Sigc_mo_tmp2=0d0;
   do iorb=1,nOrb ! Interpolation
    Sigc_mo_tmp(iorb,:,:) = 0.5d0*(Real(Sigc_mo_hh_cpx(2*iorb-1,:,:))+Real(Sigc_mo_hh_cpx(2*iorb,:,:)))
    Sigc_mo_tmp2(iorb,:,:)= 0.5d0*(Real(Sigc_mo_ee_cpx(2*iorb-1,:,:))+Real(Sigc_mo_ee_cpx(2*iorb,:,:)))
   enddo
   do iorb=1,nOrb
    Sigc_mo_hh(iorb,:)=Occ_el(iorb)*Sigc_mo_tmp(iorb,iorb,:)+(1d0-Occ_el(iorb))*Sigc_mo_tmp2(iorb,iorb,:)
   enddo
   if(offdiag0) then  ! qsGW version where all the off-diagonal elements are built at the Fermi level 
    do iorb=1,nOrb
     do jorb=1,nOrb
      if(iorb/=jorb) then
       Sigc_mo_hh(iorb,jorb) = Real(Sigc_mo_hh_cpx(nE_eval_global,iorb,jorb)) 
      endif
     enddo
    enddo
   endif
   
   ! Deallocate arrays
   deallocate(Sigc_mo_he_cpx)
   deallocate(Sigc_mo_hh_cpx)
   deallocate(Sigc_mo_eh_cpx)
   deallocate(Sigc_mo_ee_cpx)
   deallocate(E_eval_global_cpx)
   deallocate(Sigc_mo_tmp,Sigc_mo_tmp2)

  endif


end subroutine 

