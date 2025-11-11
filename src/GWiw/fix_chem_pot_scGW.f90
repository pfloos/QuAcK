subroutine fix_chem_pot_scGW(iter_fock,nBas,nfreqs,nElectrons,thrs_N,thrs_Ngrad,chem_pot,S,F_ao,Sigma_c_w_ao,wcoord,wweight, &
                             G_ao,G_ao_iw_hf,DeltaG_ao_iw,P_ao,P_ao_hf,trace_1_rdm,grad_electrons,verbose) 

! Fix the chemical potential for scGW 

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: verbose
  integer,intent(in)            :: iter_fock
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nfreqs
  double precision,intent(in)   :: nElectrons
  double precision,intent(in)   :: thrs_N
  double precision,intent(in)   :: thrs_Ngrad
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: F_ao(nBas,nBas)
  double precision,intent(in)   :: wcoord(nfreqs)
  double precision,intent(in)   :: wweight(nfreqs)
  double precision,intent(in)   :: P_ao_hf(nBas,nBas)
  complex*16,intent(in)         :: Sigma_c_w_ao(nfreqs,nBas,nBas)
  complex*16,intent(in)         :: G_ao_iw_hf(nfreqs,nBas,nBas)

! Local variables

  integer                       :: isteps
  double precision              :: thrs_closer
  double precision              :: delta_chem_pot
  double precision              :: chem_pot_change
  double precision              :: chem_pot_old
  double precision              :: trace_2up
  double precision              :: trace_up
  double precision              :: trace_down
  double precision              :: trace_2down
  double precision              :: trace_old

! Output variables

  double precision,intent(inout):: chem_pot
  double precision,intent(out)  :: trace_1_rdm
  double precision,intent(out)  :: grad_electrons
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

  if(verbose) then
   write(*,*)
   write(*,'(a,i5)') ' Fixing the Tr[1D] at scGW/scGF2 the using gradient at Fock iter ',iter_fock
   write(*,*)
  endif

  ! First approach close the value with an error lower than 1

  trace_old = 1d2
  do while( abs(trace_old-nElectrons) > thrs_closer .and. isteps <= 100 )
   isteps = isteps + 1
   call get_1rdm_scGW(nBas,nfreqs,nElectrons,chem_pot,S,F_ao,Sigma_c_w_ao, &
                      wcoord,wweight,G_ao,G_ao_iw_hf,DeltaG_ao_iw,P_ao,P_ao_hf,trace_old) 
   call get_1rdm_scGW(nBas,nfreqs,nElectrons,chem_pot-delta_chem_pot,S,F_ao,Sigma_c_w_ao, &
                      wcoord,wweight,G_ao,G_ao_iw_hf,DeltaG_ao_iw,P_ao,P_ao_hf,trace_down) 
   call get_1rdm_scGW(nBas,nfreqs,nElectrons,chem_pot+delta_chem_pot,S,F_ao,Sigma_c_w_ao, &
                      wcoord,wweight,G_ao,G_ao_iw_hf,DeltaG_ao_iw,P_ao,P_ao_hf,trace_up) 
   if( abs(trace_up-nElectrons) > abs(trace_old-nElectrons) .and. abs(trace_down-nElectrons) > abs(trace_old-nElectrons) ) then
     delta_chem_pot = 0.5d0*delta_chem_pot
     thrs_closer = 0.5d0*thrs_closer
     if(delta_chem_pot<1d-2) exit
   else
     if( abs(trace_up-nElectrons) < abs(trace_old-nElectrons) ) then
      chem_pot=chem_pot+delta_chem_pot
     else
      if( abs(trace_down-nElectrons) < abs(trace_old-nElectrons) ) then
       chem_pot=chem_pot-delta_chem_pot
      endif
     endif
   endif
  enddo

  ! Do  final search

  isteps = 0
  delta_chem_pot  = 1.0d-3
  do while( abs(trace_1_rdm) > thrs_N .and. abs(grad_electrons) > thrs_Ngrad .and. isteps <= 100 )
   isteps = isteps + 1
   chem_pot = chem_pot + chem_pot_change
   call get_1rdm_scGW(nBas,nfreqs,nElectrons,chem_pot,S,F_ao,Sigma_c_w_ao, &
                      wcoord,wweight,G_ao,G_ao_iw_hf,DeltaG_ao_iw,P_ao,P_ao_hf,trace_1_rdm)
   call get_1rdm_scGW(nBas,nfreqs,nElectrons,chem_pot+2d0*delta_chem_pot,S,F_ao,Sigma_c_w_ao, &
                      wcoord,wweight,G_ao,G_ao_iw_hf,DeltaG_ao_iw,P_ao,P_ao_hf,trace_2up)
   call get_1rdm_scGW(nBas,nfreqs,nElectrons,chem_pot+delta_chem_pot,S,F_ao,Sigma_c_w_ao, &
                      wcoord,wweight,G_ao,G_ao_iw_hf,DeltaG_ao_iw,P_ao,P_ao_hf,trace_up)
   call get_1rdm_scGW(nBas,nfreqs,nElectrons,chem_pot-delta_chem_pot,S,F_ao,Sigma_c_w_ao, &
                      wcoord,wweight,G_ao,G_ao_iw_hf,DeltaG_ao_iw,P_ao,P_ao_hf,trace_down)
   call get_1rdm_scGW(nBas,nfreqs,nElectrons,chem_pot-2d0*delta_chem_pot,S,F_ao,Sigma_c_w_ao, &
                      wcoord,wweight,G_ao,G_ao_iw_hf,DeltaG_ao_iw,P_ao,P_ao_hf,trace_2down)
   trace_1_rdm=(trace_1_rdm-nElectrons)**2d0
   trace_2up  =(trace_2up  -nElectrons)**2d0
   trace_up   =(trace_up   -nElectrons)**2d0
   trace_down =(trace_down -nElectrons)**2d0
   trace_2down=(trace_2down-nElectrons)**2d0
!   grad_electrons = (trace_up-trace_down)/(2d0*delta_chem_pot)
   grad_electrons = (-trace_2up+8d0*trace_up-8d0*trace_down+trace_2down)/(12d0*delta_chem_pot)
   chem_pot_change = -trace_1_rdm/(grad_electrons+1d-10)
   ! Maximum change is bounded within +/- 0.10
   chem_pot_change = max( min( chem_pot_change , 0.1d0 / real(isteps) ), -0.1d0 / real(isteps) )
  enddo
  call get_1rdm_scGW(nBas,nfreqs,nElectrons,chem_pot,S,F_ao,Sigma_c_w_ao, &
                     wcoord,wweight,G_ao,G_ao_iw_hf,DeltaG_ao_iw,P_ao,P_ao_hf,trace_1_rdm) 

end subroutine
