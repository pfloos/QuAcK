subroutine fix_chem_pot_scGX_bisec(iter_fock,nBas,nfreqs,nElectrons,thrs_N,thrs_Ngrad,chem_pot,S,F_ao,Sigma_c_w_ao,wcoord,wweight, &
                                   G_ao,G_ao_iw_hf,DeltaG_ao_iw,P_ao,P_ao_hf,trace_1_rdm,chem_pot_hf,verbose) 

! Fix the chemical potential for scGX 

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: verbose
  integer,intent(in)            :: iter_fock
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nfreqs
  double precision,intent(in)   :: chem_pot_hf
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

  logical                       :: normal
  integer                       :: isteps
  double precision              :: thrs_closer
  double precision              :: delta_chem_pot
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

  normal = .true.
  isteps = 0
  delta_chem_pot  = 2d-1
  thrs_closer     = 2d-1
  grad_electrons  = 0d0
  trace_1_rdm     = -1d0

  if(verbose) then
   write(*,*)
   write(*,'(a,i5)') ' Fixing the Tr[1D] at scGX at Fock iter ',iter_fock
   write(*,*)
  endif
  call get_1rdm_scGX(nBas,nfreqs,chem_pot,S,F_ao,Sigma_c_w_ao,wcoord,wweight,&
                     G_ao,G_ao_iw_hf,DeltaG_ao_iw,P_ao,P_ao_hf,trace_old) 
  if(verbose) then
   write(*,*)'------------------------------------------------------'
   write(*,'(1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,A15,2X,A1)') &
           '|','Tr[1D]','|','Chem. Pot.','|','Grad N','|'
   write(*,*)'------------------------------------------------------'
   write(*,'(1X,A1,F16.10,1X,A1,F16.10,1X,A1,F16.10,1X,A1)') &
   '|',trace_old,'|',chem_pot,'|',grad_electrons,'|'
   write(*,*)'------------------------------------------------------'
  endif

  ! First approach close the value with an error lower than 1

  trace_old = 1d2
  do while( abs(trace_old-nElectrons) > thrs_closer .and. isteps <= 100 )
   isteps = isteps + 1
   call get_1rdm_scGX(nBas,nfreqs,chem_pot,S,F_ao,Sigma_c_w_ao, &
                      wcoord,wweight,G_ao,G_ao_iw_hf,DeltaG_ao_iw,P_ao,P_ao_hf,trace_old) 
   call get_1rdm_scGX(nBas,nfreqs,chem_pot-delta_chem_pot,S,F_ao,Sigma_c_w_ao, &
                      wcoord,wweight,G_ao,G_ao_iw_hf,DeltaG_ao_iw,P_ao,P_ao_hf,trace_down) 
   call get_1rdm_scGX(nBas,nfreqs,chem_pot+delta_chem_pot,S,F_ao,Sigma_c_w_ao, &
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
  delta_chem_pot  = 2d-2
  call get_1rdm_scGX(nBas,nfreqs,chem_pot,S,F_ao,Sigma_c_w_ao, &
                     wcoord,wweight,G_ao,G_ao_iw_hf,DeltaG_ao_iw,P_ao,P_ao_hf,trace_1_rdm)
  isteps = 1
  chem_pot_up   = chem_pot + delta_chem_pot
  chem_pot_down = chem_pot - delta_chem_pot
  do
   call get_1rdm_scGX(nBas,nfreqs,chem_pot_up,S,F_ao,Sigma_c_w_ao, &
                      wcoord,wweight,G_ao,G_ao_iw_hf,DeltaG_ao_iw,P_ao,P_ao_hf,trace_up)
   call get_1rdm_scGX(nBas,nfreqs,chem_pot_down,S,F_ao,Sigma_c_w_ao, &
                      wcoord,wweight,G_ao,G_ao_iw_hf,DeltaG_ao_iw,P_ao,P_ao_hf,trace_down)

   if((trace_up>trace_1_rdm .and. trace_down>trace_1_rdm) .or. (trace_up<trace_1_rdm .and. trace_down<trace_1_rdm)) then
    write(*,'(a,f16.10)') ' Chemical potential          ',chem_pot
    write(*,'(a,f16.10)') ' Delta chemical potential    ',delta_chem_pot
    write(*,'(a,f16.10)') ' Trace with chem pot         ',trace_1_rdm
    write(*,'(a,f16.10)') ' Trace with chem pot + delta ',trace_up
    write(*,'(a,f16.10)') ' Trace with chem pot - delta ',trace_down
    chem_pot=chem_pot_hf
    write(*,'(a,f10.5,a)') ' Entering the gradient method to try to optimize the chem. potential [ chem_pot guess ',chem_pot,' ]'
    call fix_chem_pot_scGX(iter_fock,nBas,nfreqs,nElectrons,thrs_N,thrs_Ngrad,chem_pot,S,F_ao,    & 
                           Sigma_c_w_ao,wcoord,wweight,G_ao,G_ao_iw_hf,DeltaG_ao_iw,P_ao,P_ao_hf, &
                           trace_1_rdm,grad_electrons,verbose)
    if((abs(trace_1_rdm-nElectrons))**2d0 > thrs_N) then
     write(*,'(a,f16.10)') ' Chemical potential          ',chem_pot
     write(*,'(a,f16.10)') ' Trace with chem pot         ',trace_1_rdm
     write(*,'(a,f16.10)') ' Grad. Nelectrons            ',grad_electrons
     stop
    else
     if(verbose) then
      write(*,'(1X,A1,F16.10,1X,A1,F16.10,1X,A1,F16.10,1X,A1)') &
      '|',trace_1_rdm,'|',chem_pot,'|',grad_electrons,'|'
      write(*,*)'------------------------------------------------------'
      write(*,*)
     endif
     return
    endif
   endif

   if(trace_1_rdm<nElectrons) then
    if(trace_up>nElectrons .or. trace_down>nElectrons) then
     if(trace_up>nElectrons) then
      chem_pot_down=chem_pot
      normal=.true.     
      exit
     endif
     if(trace_down>nElectrons) then
      chem_pot_up=chem_pot
      normal=.false.     
      exit
     endif
    endif
   else
    if(trace_up<nElectrons .or. trace_down<nElectrons) then
     if(trace_down<nElectrons) then
      chem_pot_up=chem_pot
      normal=.true.     
      exit
     endif
     if(trace_up<nElectrons) then
      chem_pot_down=chem_pot
      normal=.false.     
      exit
     endif
    endif
   endif
   isteps = isteps + 1

   if(isteps==100) exit

   delta_chem_pot=delta_chem_pot + 2d-2
   chem_pot_up   = chem_pot + delta_chem_pot
   chem_pot_down = chem_pot - delta_chem_pot
   
 enddo

!   call get_1rdm_scGX(nBas,nfreqs,chem_pot_up,S,F_ao,Sigma_c_w_ao, &
!                      wcoord,wweight,G_ao,G_ao_iw_hf,DeltaG_ao_iw,P_ao,P_ao_hf,trace_up)
!   call get_1rdm_scGX(nBas,nfreqs,chem_pot_down,S,F_ao,Sigma_c_w_ao, &
!                      wcoord,wweight,G_ao,G_ao_iw_hf,DeltaG_ao_iw,P_ao,P_ao_hf,trace_down)
!   write(*,*) delta_chem_pot,chem_pot_down,chem_pot_up,trace_down,trace_up,normal

!  write(*,*) "| doing bisection method ...                          |"
  ! Do bisection method
  isteps = 0
  do while( abs(trace_1_rdm-nElectrons) > thrs_N .and. isteps <= 100 )
   isteps = isteps + 1
   chem_pot=0.5d0*(chem_pot_up+chem_pot_down)
   call get_1rdm_scGX(nBas,nfreqs,chem_pot,S,F_ao,Sigma_c_w_ao, &
                      wcoord,wweight,G_ao,G_ao_iw_hf,DeltaG_ao_iw,P_ao,P_ao_hf,trace_1_rdm) 
!   write(*,'(1X,A1,F16.10,1X,A1,F16.10,1X,A1F16.10,1X,A1)') &
!    '|',trace_1_rdm,'|',chem_pot,'|',grad_electrons,'|'
   if(normal) then
    if(trace_1_rdm>nElectrons) then
     chem_pot_up=chem_pot
    else
     chem_pot_down=chem_pot
    endif
   else
    if(trace_1_rdm<nElectrons) then
     chem_pot_up=chem_pot
    else
     chem_pot_down=chem_pot
    endif
   endif
  enddo

  ! Print info
!  write(*,*)'------------------------------------------------------'
  call get_1rdm_scGX(nBas,nfreqs,chem_pot,S,F_ao,Sigma_c_w_ao, &
                     wcoord,wweight,G_ao,G_ao_iw_hf,DeltaG_ao_iw,P_ao,P_ao_hf,trace_1_rdm) 
!  write(*,'(1X,A1,1X,A15,1X,A1,1X,A15,1X,A1A15,2X,A1)') &
!          '|','Tr[1D]','|','Chem. Pot.','|','Grad N','|'
!  write(*,*)'------------------------------------------------------'
  if(verbose) then
   write(*,'(1X,A1,F16.10,1X,A1,F16.10,1X,A1,F16.10,1X,A1)') &
   '|',trace_1_rdm,'|',chem_pot,'|',grad_electrons,'|'
   write(*,*)'------------------------------------------------------'
   write(*,*)
  endif

end subroutine
