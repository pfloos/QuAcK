subroutine fix_chem_pot(nO,nOrb,nOrb2,thrs_N,trace_1rdm,chem_pot,H_hfb,cp,R,eHFB_)

! Fix the chemical potential to integrate to the N for 2N electron systems

  implicit none

! Input variables

  integer,intent(in)            :: nO,nOrb,nOrb2
  double precision,intent(in)   :: thrs_N

! Local variables

  logical                       :: backward
  integer                       :: iorb
  integer                       :: isteps
  double precision              :: delta_chem_pot
  double precision              :: chem_pot_change
  double precision              :: grad_electrons
  double precision              :: trace_up
  double precision              :: trace_down
  double precision              :: trace_old

! Output variables

  double precision              :: trace_1rdm
  double precision,intent(inout):: chem_pot
  double precision,intent(inout):: cp(nOrb2,nOrb2) 
  double precision,intent(inout):: R(nOrb2,nOrb2) 
  double precision,intent(inout):: eHFB_(nOrb2) 
  double precision,intent(inout):: H_hfb(nOrb2,nOrb2) 

  !  Initialize 

  backward = .false.
  isteps = 0
  delta_chem_pot  = 1.0d-1
  chem_pot_change = 0d0
  grad_electrons  = 1d0
  trace_1rdm      = -1d0

  ! Set H_HFB to its non-chemical potential dependent contribution
 
  do iorb=1,nOrb
   H_hfb(iorb,iorb)=H_hfb(iorb,iorb)+chem_pot
   H_hfb(iorb+nOrb,iorb+nOrb)=H_hfb(iorb+nOrb,iorb+nOrb)-chem_pot
  enddo

  write(*,*)
  write(*,*)'------------------------------------------------------'
  write(*,'(1X,A1,1X,A15,1X,A1,1X,A15,1X,A1A15,2X,A1)') &
          '|','Tr[1D]','|','Chem. Pot.','|','Grad N','|'
  write(*,*)'------------------------------------------------------'

  ! First approach close the value with an error lower than 1

  call diag_H_hfb(nOrb,nOrb2,chem_pot,trace_old,H_hfb,cp,R,eHFB_)
  write(*,'(1X,A1,F16.10,1X,A1,F16.10,1X,A1F16.10,1X,A1)') &
  '|',trace_old,'|',chem_pot,'|',grad_electrons,'|'
  do while( abs(trace_1rdm-nO) > 1.0d0 .and. isteps <= 100 )
   isteps = isteps + 1
   chem_pot = chem_pot + delta_chem_pot
   call diag_H_hfb(nOrb,nOrb2,chem_pot,trace_1rdm,H_hfb,cp,R,eHFB_)
   write(*,'(1X,A1,F16.10,1X,A1,F16.10,1X,A1F16.10,1X,A1)') &
   '|',trace_1rdm,'|',chem_pot,'|',grad_electrons,'|'
   if( abs(trace_1rdm-nO) > abs(trace_old-nO) .and. .not.backward ) then
    backward=.true.
    chem_pot = chem_pot - 2d0*delta_chem_pot
    delta_chem_pot=-delta_chem_pot
   endif
  enddo

  ! Do  final search

  write(*,*)'------------------------------------------------------'
  isteps = 0
  delta_chem_pot  = 1.0d-3
  do while( abs(trace_1rdm-nO) > thrs_N .and. isteps <= 100 )
   isteps = isteps + 1
   chem_pot = chem_pot + chem_pot_change
   call diag_H_hfb(nOrb,nOrb2,chem_pot,trace_1rdm,H_hfb,cp,R,eHFB_)
   call diag_H_hfb(nOrb,nOrb2,chem_pot+delta_chem_pot,trace_up,H_hfb,cp,R,eHFB_)
   call diag_H_hfb(nOrb,nOrb2,chem_pot-delta_chem_pot,trace_down,H_hfb,cp,R,eHFB_)
   grad_electrons = (trace_up-trace_down)/(2.0d0*delta_chem_pot)
   chem_pot_change = -(trace_1rdm-nO)/(grad_electrons+1d-10)
   write(*,'(1X,A1,F16.10,1X,A1,F16.10,1X,A1F16.10,1X,A1)') &
   '|',trace_1rdm,'|',chem_pot,'|',grad_electrons,'|'
  enddo
  write(*,*)'------------------------------------------------------'
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

