subroutine fix_chem_pot_scGW_bisec(iter_fock,nBas,nfreqs,nElectrons,thrs_N,thrs_Ngrad,chem_pot,S,F_ao,Sigma_c_w_ao,wcoord,wweight, &
                                   G_ao,G_ao_iw_hf,DeltaG_ao_iw,P_ao,P_ao_hf,trace_1_rdm) 

! Fix the chemical potential for scGW 

  implicit none
  include 'parameters.h'

! Input variables

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
  double precision              :: grad_electrons
  double precision              :: chem_pot_up
  double precision              :: chem_pot_down
  double precision              :: trace_up
  double precision              :: trace_down
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
  grad_electrons  = 0d0
  trace_1_rdm      = -1d0

  write(*,*)
  write(*,'(a,i5)') ' Fixing the Tr[1D] at scGW/scGF2 at Fock iter ',iter_fock
  write(*,*)
  call get_1rdm_scGW(nBas,nfreqs,nElectrons,chem_pot,S,F_ao,Sigma_c_w_ao, &
                     wcoord,wweight,G_ao,G_ao_iw_hf,DeltaG_ao_iw,P_ao,P_ao_hf,trace_old) 
  write(*,*)'------------------------------------------------------'
  write(*,'(1X,A1,1X,A15,1X,A1,1X,A15,1X,A1A15,2X,A1)') &
          '|','Tr[1D]','|','Chem. Pot.','|','Grad N','|'
  write(*,*)'------------------------------------------------------'
  write(*,'(1X,A1,F16.10,1X,A1,F16.10,1X,A1F16.10,1X,A1)') &
  '|',trace_old,'|',chem_pot,'|',grad_electrons,'|'
  write(*,*)'------------------------------------------------------'

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
!     write(*,'(1X,A1,F16.10,1X,A1,F16.10,1X,A1F16.10,1X,A1)') &
!     '|',trace_old,'|',chem_pot,'|',grad_electrons,'|'
     delta_chem_pot = 0.5d0*delta_chem_pot
     thrs_closer = 0.5d0*thrs_closer
!     write(*,*) "| contracting ...                                     |"
     if(delta_chem_pot<1d-2) exit
   else
     if( abs(trace_up-nElectrons) < abs(trace_old-nElectrons) ) then
      chem_pot=chem_pot+delta_chem_pot
!      write(*,'(1X,A1,F16.10,1X,A1,F16.10,1X,A1F16.10,1X,A1)') &
!      '|',trace_up,'|',chem_pot,'|',grad_electrons,'|'
     else
      if( abs(trace_down-nElectrons) < abs(trace_old-nElectrons) ) then
       chem_pot=chem_pot-delta_chem_pot
!       write(*,'(1X,A1,F16.10,1X,A1,F16.10,1X,A1F16.10,1X,A1)') &
!       '|',trace_down,'|',chem_pot,'|',grad_electrons,'|'
      endif
     endif
   endif
  enddo

  ! Find bounds for the bisection method
!  write(*,*) "| fiding bisection method bounds ...                  |"
  call get_1rdm_scGW(nBas,nfreqs,nElectrons,chem_pot,S,F_ao,Sigma_c_w_ao, &
                     wcoord,wweight,G_ao,G_ao_iw_hf,DeltaG_ao_iw,P_ao,P_ao_hf,trace_1_rdm) 
  isteps = 1
  delta_chem_pot  = 1d-2
  trace_up   = 0d0
  trace_down = 0d0
  chem_pot_up   = 0d0
  chem_pot_down = 0d0
  do while( abs(trace_1_rdm-nElectrons) > thrs_N .and. isteps <= 100 )
   if(trace_1_rdm>nElectrons) then
    chem_pot_up=chem_pot+isteps*delta_chem_pot
    call get_1rdm_scGW(nBas,nfreqs,nElectrons,chem_pot_up,S,F_ao,Sigma_c_w_ao, &
                       wcoord,wweight,G_ao,G_ao_iw_hf,DeltaG_ao_iw,P_ao,P_ao_hf,trace_up) 
    if(trace_up<nElectrons) exit
   else    
    chem_pot_down=chem_pot-isteps*delta_chem_pot
    call get_1rdm_scGW(nBas,nfreqs,nElectrons,chem_pot_down,S,F_ao,Sigma_c_w_ao, &
                       wcoord,wweight,G_ao,G_ao_iw_hf,DeltaG_ao_iw,P_ao,P_ao_hf,trace_down) 
    if(trace_down>nElectrons) exit
   endif

   isteps = isteps + 1
  enddo
  if(abs(chem_pot_up)>1d-8) then
   chem_pot_down=chem_pot
   trace_down=trace_1_rdm
!   write(*,'(1X,A1,F16.10,1X,A1,F16.10,1X,A1F16.10,1X,A1)') &
!    '|',trace_up,'|',chem_pot_up,'|',grad_electrons,'|'
  else
   chem_pot_up=chem_pot
   trace_up=trace_1_rdm
!   write(*,'(1X,A1,F16.10,1X,A1,F16.10,1X,A1F16.10,1X,A1)') &
!    '|',trace_down,'|',chem_pot_down,'|',grad_electrons,'|'
  endif

!  write(*,*) "| doing bisection method ...                          |"
  ! Do bisection method
  isteps = 0
  do while( abs(trace_1_rdm-nElectrons) > thrs_N .and. isteps <= 100 )
   isteps = isteps + 1
   chem_pot=0.5d0*(chem_pot_up+chem_pot_down)
   call get_1rdm_scGW(nBas,nfreqs,nElectrons,chem_pot,S,F_ao,Sigma_c_w_ao, &
                      wcoord,wweight,G_ao,G_ao_iw_hf,DeltaG_ao_iw,P_ao,P_ao_hf,trace_1_rdm) 
!   write(*,'(1X,A1,F16.10,1X,A1,F16.10,1X,A1F16.10,1X,A1)') &
!    '|',trace_1_rdm,'|',chem_pot,'|',grad_electrons,'|'
   if(trace_1_rdm>nElectrons) then
    chem_pot_down=chem_pot
    trace_down=trace_1_rdm
   else
    chem_pot_up=chem_pot
    trace_up=trace_1_rdm
   endif
  enddo

  ! Print info
!  write(*,*)'------------------------------------------------------'
  call get_1rdm_scGW(nBas,nfreqs,nElectrons,chem_pot,S,F_ao,Sigma_c_w_ao, &
                     wcoord,wweight,G_ao,G_ao_iw_hf,DeltaG_ao_iw,P_ao,P_ao_hf,trace_1_rdm) 
!  write(*,'(1X,A1,1X,A15,1X,A1,1X,A15,1X,A1A15,2X,A1)') &
!          '|','Tr[1D]','|','Chem. Pot.','|','Grad N','|'
!  write(*,*)'------------------------------------------------------'
  write(*,'(1X,A1,F16.10,1X,A1,F16.10,1X,A1F16.10,1X,A1)') &
  '|',trace_1_rdm,'|',chem_pot,'|',grad_electrons,'|'
  write(*,*)'------------------------------------------------------'
  write(*,*)

end subroutine
