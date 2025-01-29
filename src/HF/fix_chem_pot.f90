subroutine fix_chem_pot(nO,nOrb,nOrb2,thrs_N,trace_1rdm,chem_pot,H_hfb,cp,R,eHFB_)

! Fix the chemical potential to integrate to the N for 2N electron systems

  implicit none

! Input variables

  integer,intent(in)            :: nO,nOrb,nOrb2
  double precision,intent(in)   :: thrs_N

! Local variables

  logical                       :: is_up
  integer                       :: iorb
  integer                       :: isteps
  double precision              :: trace_curr,trace_down,trace_up
  double precision              :: chem_pot_curr,chem_pot_down,chem_pot_up
  double precision              :: delta_chem_pot
  double precision              :: golden_ratio
  double precision              :: a
  double precision              :: c

! Output variables

  double precision              :: trace_1rdm
  double precision,intent(inout):: chem_pot
  double precision,intent(inout):: cp(nOrb2,nOrb2) 
  double precision,intent(inout):: R(nOrb2,nOrb2) 
  double precision,intent(inout):: eHFB_(nOrb2) 
  double precision,intent(inout):: H_hfb(nOrb2,nOrb2) 

  !  Initialize delta_chem_pot

  is_up=.false.
  isteps=0
  golden_ratio = 1.618033988
  delta_chem_pot = 5d-1
  trace_down = 0d0
  trace_up   = 0d0
  chem_pot_down = 0d0 
  chem_pot_up   = 0d0 

  ! Set H_HFB to its non-chemical potential dependent contribution
 
  do iorb=1,nOrb
   H_hfb(iorb,iorb)=H_hfb(iorb,iorb)+chem_pot
   H_hfb(iorb+nOrb,iorb+nOrb)=H_hfb(iorb+nOrb,iorb+nOrb)-chem_pot
  enddo

  write(*,*)
  write(*,*)'-------------------------------------'
  write(*,'(1X,A1,1X,A15,1X,A1,1X,A15,1X,A1)') &
          '|','Tr[1D]','|','Chem. Pot.','|'
  write(*,*)'-------------------------------------'

  ! Set interval to search (the good chem_pot is in [chem_pot_down,chem_pot_up])
  call diag_H_hfb(nOrb,nOrb2,chem_pot,trace_curr,H_hfb,cp,R,eHFB_)
  write(*,'(1X,A1,F16.10,1X,A1,F16.10,1X,A1)') &
  '|',trace_curr,'|',chem_pot,'|'
  chem_pot_curr=chem_pot
  if(trace_curr<nO) then
   ! Increase chem_pot to occupy more orbs.
   is_up=.true.
   do
    isteps = isteps+1
    chem_pot = chem_pot + delta_chem_pot
    call diag_H_hfb(nOrb,nOrb2,chem_pot,trace_up,H_hfb,cp,R,eHFB_)
    write(*,'(1X,A1,F16.10,1X,A1,F16.10,1X,A1)') &
    '|',trace_up,'|',chem_pot,'|'
    if(trace_up>=nO+1 .or. isteps>100) exit ! max 50 au steps for mu (is a lot)
   enddo
   chem_pot_up=chem_pot
  else
   ! Decrease chem_pot to occupy less orbs.
   do
    isteps = isteps+1
    chem_pot = chem_pot - delta_chem_pot
    call diag_H_hfb(nOrb,nOrb2,chem_pot,trace_down,H_hfb,cp,R,eHFB_)
    write(*,'(1X,A1,F16.10,1X,A1,F16.10,1X,A1)') &
    '|',trace_down,'|',chem_pot,'|'
    if(trace_down<=nO-1 .or. isteps>100) exit ! max 50 au steps for mu (is a lot)
   enddo
   chem_pot_down=chem_pot
  endif
  if(is_up) then
   write(*,'(1X,A1,F16.10,1X,A1,F16.10,1X,A1)') &
    '|',trace_curr,'|',chem_pot_curr,'|'
   write(*,'(1X,A1,F16.10,1X,A1,F16.10,1X,A1)') &
    '|',trace_up  ,'|',chem_pot_up  ,'|'
   trace_down=trace_curr
   chem_pot_down=chem_pot_curr
  else
   write(*,'(1X,A1,F16.10,1X,A1,F16.10,1X,A1)') &
    '|',trace_curr,'|',chem_pot_curr,'|'
   write(*,'(1X,A1,F16.10,1X,A1,F16.10,1X,A1)') &
    '|',trace_down,'|',chem_pot_down,'|'
   trace_up=trace_curr
   chem_pot_up=chem_pot_curr
  endif

  ! Use Golden-section search algorithm to find chem_pot
  isteps=0
  do
    isteps = isteps+1
    a=(chem_pot_up-chem_pot_down)/(1d0+golden_ratio)
    c=a/golden_ratio
    chem_pot_curr=chem_pot_down+a
    chem_pot=chem_pot_curr+c
    call diag_H_hfb(nOrb,nOrb2,chem_pot_curr,trace_curr,H_hfb,cp,R,eHFB_)
    write(*,'(1X,A1,F16.10,1X,A1,F16.10,1X,A1)') &
     '|',trace_curr,'|',chem_pot_curr,'|'
    call diag_H_hfb(nOrb,nOrb2,chem_pot,trace_1rdm,H_hfb,cp,R,eHFB_)
    write(*,'(1X,A1,F16.10,1X,A1,F16.10,1X,A1)') &
     '|',trace_1rdm,'|',chem_pot,'|'
    if(abs(trace_1rdm-nO)<thrs_N .or. isteps>1000) exit ! 1000 steps for finding mu
    if(is_up) then
     if(trace_1rdm>=trace_curr .or. abs(trace_1rdm-trace_curr)<1d-8 ) then
      chem_pot_down=chem_pot_curr
      trace_down=trace_curr
     else
      chem_pot_up=chem_pot
      trace_up=trace_1rdm
     endif
    else
     if(trace_1rdm>=trace_curr .or. abs(trace_1rdm-trace_curr)<1d-8 ) then
      chem_pot_up=chem_pot
      trace_up=trace_1rdm
     else
      chem_pot_down=chem_pot_curr
      trace_down=trace_curr
     endif
    endif
  enddo

  write(*,*)'-------------------------------------'
  write(*,*)
  
  ! Reset H_HFB to its chemical potential version
 
  do iorb=1,nOrb
   H_hfb(iorb,iorb)=H_hfb(iorb,iorb)-chem_pot
   H_hfb(iorb+nOrb,iorb+nOrb)=H_hfb(iorb+nOrb,iorb+nOrb)+chem_pot
  enddo

end subroutine 

subroutine diag_H_hfb(nOrb,nOrb2,chem_pot,trace_1rdm,H_hfb,cp,R,eHFB_)

! Fix the chemical potential to integrate to the N for 2N electron systems

  implicit none

! Input variables

  integer,intent(in)            :: nOrb,nOrb2
  double precision,intent(in)   :: H_hfb(nOrb2,nOrb2) 
  double precision,intent(in)   :: chem_pot

! Local variables

  integer                       :: iorb

! Output variables

  double precision,intent(out)  :: trace_1rdm
  double precision,intent(inout):: cp(nOrb2,nOrb2) 
  double precision,intent(inout):: R(nOrb2,nOrb2) 
  double precision,intent(inout):: eHFB_(nOrb2) 
    
  cp(:,:) = H_hfb(:,:)
  do iorb=1,nOrb
   cp(iorb,iorb) = cp(iorb,iorb) - chem_pot 
   cp(iorb+nOrb,iorb+nOrb) = cp(iorb+nOrb,iorb+nOrb) + chem_pot 
  enddo

  ! Diagonalize H_HFB matrix

  call diagonalize_matrix(nOrb2,cp,eHFB_)
  
  ! Build R and extract P and Panom
    
  trace_1rdm = 0d0 
  R(:,:)     = 0d0
  do iorb=1,nOrb
   R(:,:) = R(:,:) + matmul(cp(:,iorb:iorb),transpose(cp(:,iorb:iorb))) 
  enddo
  do iorb=1,nOrb
   trace_1rdm = trace_1rdm + R(iorb,iorb) 
  enddo

end subroutine 

