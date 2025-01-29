subroutine fix_chem_pot(nO,nOrb,nOrb2,thrs_N,trace_1rdm,chem_pot,H_hfb,cp,R,eHFB_)

! Fix the chemical potential to integrate to the N for 2N electron systems

  implicit none

! Input variables

  integer,intent(in)            :: nO,nOrb,nOrb2
  double precision,intent(in)   :: thrs_N

! Local variables

  integer                       :: iorb
  double precision              :: trace_curr,trace_down,trace_up
  double precision              :: chem_pot_curr,chem_pot_down,chem_pot_up
  double precision              :: delta_chem_pot

! Output variables

  double precision              :: trace_1rdm
  double precision,intent(inout):: chem_pot
  double precision,intent(inout):: cp(nOrb2,nOrb2) 
  double precision,intent(inout):: R(nOrb2,nOrb2) 
  double precision,intent(inout):: eHFB_(nOrb2) 
  double precision,intent(inout):: H_hfb(nOrb2,nOrb2) 

  !  Initialize delta_chem_pot

  delta_chem_pot = 1d0
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

  ! Set interval to search
  call diag_H_hfb(nOrb,nOrb2,chem_pot,trace_curr,H_hfb,cp,R,eHFB_)
  chem_pot_curr=chem_pot
  if(trace_curr<nO) then
   ! Increase chem_pot to occupy more orbs.
   do
    chem_pot = chem_pot + delta_chem_pot
    call diag_H_hfb(nOrb,nOrb2,chem_pot,trace_up,H_hfb,cp,R,eHFB_)
    if(trace_up>=nO+1) exit
   enddo
   chem_pot_up=chem_pot
  else
   ! Decrease chem_pot to occupy less orbs.
   do
    chem_pot = chem_pot - delta_chem_pot
    call diag_H_hfb(nOrb,nOrb2,chem_pot,trace_down,H_hfb,cp,R,eHFB_)
    if(trace_down<=nO-1) exit
   enddo
   chem_pot_down=chem_pot
  endif

  if(abs(chem_pot_up)>1e-4) then
   
   write(*,'(1X,A1,F16.10,1X,A1,F16.10,1X,A1)') &
    '|',trace_curr,'|',chem_pot_curr,'|'
   write(*,'(1X,A1,F16.10,1X,A1,F16.10,1X,A1)') &
    '|',trace_up  ,'|',chem_pot_up  ,'|'


   write(*,*)'-------------------------------------'
   write(*,*)

  endif

  if(abs(chem_pot_down)>1e-4) then
   
   write(*,'(1X,A1,F16.10,1X,A1,F16.10,1X,A1)') &
    '|',trace_curr,'|',chem_pot_curr,'|'
   write(*,'(1X,A1,F16.10,1X,A1,F16.10,1X,A1)') &
    '|',trace_down,'|',chem_pot_down,'|'


   write(*,*)'-------------------------------------'
   write(*,*)

  endif
  
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

