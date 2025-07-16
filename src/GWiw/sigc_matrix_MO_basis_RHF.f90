subroutine sigc_MO_basis_RHF(nOrb,nO,offdiag0,eta,shift,eqsGW_state,vMAT,nfreqs,ntimes, &
                             wcoord,wweight,Sigc_mo)

! Compute Sigma_c matrix in the AO basis

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: offdiag0

  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nO
  integer,intent(in)            :: nfreqs
  integer,intent(in)            :: ntimes

  double precision,intent(in)   :: eta
  double precision,intent(in)   :: shift
  double precision,intent(in)   :: wcoord(nfreqs),wweight(nfreqs)
  double precision,intent(in)   :: eqsGW_state(nOrb)
  double precision,intent(in)   :: vMAT(nOrb*nOrb,nOrb*nOrb)

! Local variables

  integer                       :: iorb,jorb
  integer                       :: nE_eval_global

  double precision              :: chem_pot
  double precision,allocatable  :: Sigc_mo_tmp(:,:,:)

  complex *16,allocatable       :: Sigc_mo_cpx(:,:,:)
  complex *16,allocatable       :: E_eval_global_cpx(:)

! Output variables

  double precision,intent(out)  :: Sigc_mo(nOrb,nOrb)

! Initialize variables
  
  chem_pot = 0.5d0*(eqsGW_state(nO)+eqsGW_state(nO+1))

! Se energies using cluster method or just using two shifts
  if(.false.) then

   allocate(E_eval_global_cpx(1))
   call set_Eeval_cluster(nOrb,nOrb,1,shift,eqsGW_state,nE_eval_global,&
                          E_eval_global_cpx,chem_pot)
   deallocate(E_eval_global_cpx)
   allocate(E_eval_global_cpx(nE_eval_global))
   call set_Eeval_cluster(nOrb,nOrb,nE_eval_global,shift,eqsGW_state,&
                           nE_eval_global,E_eval_global_cpx,chem_pot)

   !  Run over energies
   allocate(Sigc_mo_cpx(nE_eval_global,nOrb,nOrb))
   call build_Sigmac_w_RHF(nOrb,nO,nE_eval_global,eta,0,E_eval_global_cpx,eqsGW_state,nfreqs,ntimes,&
                           wweight,wcoord,vMAT,Sigc_mo_cpx)

   deallocate(E_eval_global_cpx)

   deallocate(Sigc_mo_cpx)

   Sigc_mo=0d0 ! TODO
   
  else

   nE_eval_global=2*nOrb+1
   allocate(E_eval_global_cpx(nE_eval_global))
   do iorb=1,nOrb
    E_eval_global_cpx(2*iorb-1)=eqsGW_state(iorb)-shift-chem_pot
    E_eval_global_cpx(2*iorb)  =eqsGW_state(iorb)+shift-chem_pot
   enddo
   E_eval_global_cpx(nE_eval_global)=czero

   !  Run over energies
   allocate(Sigc_mo_cpx(nE_eval_global,nOrb,nOrb))
   call build_Sigmac_w_RHF(nOrb,nO,nE_eval_global,eta,0,E_eval_global_cpx,eqsGW_state,nfreqs,ntimes,&
                           wweight,wcoord,vMAT,Sigc_mo_cpx)
   deallocate(E_eval_global_cpx)
   
   ! Interpolate Sigma MO [incl. version 2 for qsGW]
   allocate(Sigc_mo_tmp(nOrb,nOrb,nOrb))
   do iorb=1,nOrb
    Sigc_mo_tmp(iorb,:,:)=0.5d0*(Real(Sigc_mo_cpx(2*iorb-1,:,:))+Real(Sigc_mo_cpx(2*iorb,:,:)))
   enddo
   do iorb=1,nOrb
    Sigc_mo(iorb,:)=Sigc_mo_tmp(iorb,iorb,:)
   enddo
   if(offdiag0) then  ! qsGW version where all the off-diagonal elements are built at the Fermi level 
    do iorb=1,nOrb
     do jorb=1,nOrb
      if(iorb/=jorb) Sigc_mo(iorb,jorb) = Real(Sigc_mo_cpx(nE_eval_global,iorb,jorb))
     enddo
    enddo
   endif

   deallocate(Sigc_mo_cpx)

  endif

end subroutine 

