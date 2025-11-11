subroutine fix_chem_pot(nO,nOrb,nOrb_twice,nSCF,thrs_N,trace_1rdm,chem_pot,H_hfb,cp,R,eHFB_)

! Fix the chemical potential to integrate to the N for 2N electron systems

  implicit none

! Input variables

  integer,intent(in)            :: nO,nOrb,nOrb_twice,nSCF
  double precision,intent(in)   :: thrs_N

! Local variables

  logical                       :: use_nelectrons
  integer                       :: iorb
  integer                       :: isteps
  double precision              :: nO_
  double precision              :: thrs_Ngrad
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
  double precision,allocatable  :: R_tmp(:,:) 
  double precision,allocatable  :: cp_tmp(:,:) 

! Output variables

  double precision              :: trace_1rdm
  double precision,intent(inout):: chem_pot
  double precision,intent(inout):: cp(nOrb_twice,nOrb_twice) 
  double precision,intent(inout):: R(nOrb_twice,nOrb_twice) 
  double precision,intent(inout):: eHFB_(nOrb_twice) 
  double precision,intent(inout):: H_hfb(nOrb_twice,nOrb_twice) 

  !  Initialize 

  isteps = 0
  delta_chem_pot  = 2d-1
  thrs_closer     = 2d-1
  thrs_Ngrad      = 1d-8
  chem_pot_change = 0d0
  grad_electrons  = 1d0
  trace_1rdm      = -1d0
  nO_             = nO
  inquire(file='Nelectrons_RHFB', exist=use_nelectrons)
  if(use_nelectrons) then
    write(*,*) 'File Nelectrons_RHFB encountered, setting nO = nO_read/2'
    open(unit=314, form='formatted', file='Nelectrons_RHFB', status='old')
    read(314,*) nO_
    close(314)
    nO_=0.5d0*nO_
    write(*,'(a,f10.5,a)') ' Using ',nO_,' electrons in RHFB per spin channel'
  endif
  allocate(R_tmp(nOrb_twice,nOrb_twice))
  allocate(cp_tmp(nOrb_twice,nOrb_twice))

  ! Set H_HFB to its non-chemical potential dependent contribution
 
  do iorb=1,nOrb
   H_hfb(iorb,iorb)=H_hfb(iorb,iorb)+chem_pot
   H_hfb(iorb+nOrb,iorb+nOrb)=H_hfb(iorb+nOrb,iorb+nOrb)-chem_pot
  enddo

!  write(*,*)
!  write(*,'(a,I5)') ' Fixing the Tr[1D] at iteration ',nSCF
!  write(*,*)
!  write(*,*)'------------------------------------------------------'
!  write(*,'(1X,A1,1X,A15,1X,A1,1X,A15,1X,A1A15,2X,A1)') &
!          '|','Tr[1D]','|','Chem. Pot.','|','Grad N','|'
!  write(*,*)'------------------------------------------------------'

  ! First approach close the value with an error lower than 1

  trace_old = 1d2
  do while( abs(trace_old-nO_) > thrs_closer .and. isteps <= 100 )
   isteps = isteps + 1
   call diag_H_hfb(nOrb,nOrb_twice,chem_pot,trace_old,H_hfb,cp,R,eHFB_)
   call diag_H_hfb(nOrb,nOrb_twice,chem_pot-delta_chem_pot,trace_down,H_hfb,cp,R,eHFB_)
   call diag_H_hfb(nOrb,nOrb_twice,chem_pot+delta_chem_pot,trace_up,H_hfb,cp,R,eHFB_)
   if( abs(trace_up-nO_) > abs(trace_old-nO_) .and. abs(trace_down-nO_) > abs(trace_old-nO_) ) then
!     write(*,'(1X,A1,F16.10,1X,A1,F16.10,1X,A1F16.10,1X,A1)') &
!     '|',trace_old,'|',chem_pot,'|',grad_electrons,'|'
     delta_chem_pot = 0.75d0*delta_chem_pot
     thrs_closer = 0.5d0*thrs_closer
!     write(*,*) "| contracting ...                                     |"  
     if(delta_chem_pot<1d-2) exit
   else
     if( abs(trace_up-nO_) < abs(trace_old-nO_) ) then
      chem_pot=chem_pot+delta_chem_pot 
!      write(*,'(1X,A1,F16.10,1X,A1,F16.10,1X,A1F16.10,1X,A1)') &
!      '|',trace_up,'|',chem_pot,'|',grad_electrons,'|'
     else
      if( abs(trace_down-nO_) < abs(trace_old-nO_) ) then
       chem_pot=chem_pot-delta_chem_pot
!       write(*,'(1X,A1,F16.10,1X,A1,F16.10,1X,A1F16.10,1X,A1)') &
!       '|',trace_down,'|',chem_pot,'|',grad_electrons,'|'
      endif
     endif
   endif
  enddo

  ! Do  final search

!  write(*,*)'------------------------------------------------------'
!  write(*,'(1X,A1,1X,A15,1X,A1,1X,A15,1X,A1A15,2X,A1)') &
!          '|','Error Tr[1D]','|','Chem. Pot.','|','Grad N','|'
!  write(*,*)'------------------------------------------------------'
  isteps = 0
  delta_chem_pot  = 1.0d-3
  trace_1rdm=(trace_1rdm-nO_)**2d0
  do while( sqrt(trace_1rdm) > thrs_N .and. abs(grad_electrons) > thrs_Ngrad .and. isteps <= 100 )
   isteps = isteps + 1
   chem_pot = chem_pot + chem_pot_change
   call diag_H_hfb(nOrb,nOrb_twice,chem_pot,trace_1rdm,H_hfb,cp,R,eHFB_)
   call diag_H_hfb(nOrb,nOrb_twice,chem_pot+2d0*delta_chem_pot,trace_2up,H_hfb,cp_tmp,R_tmp,eHFB_)
   call diag_H_hfb(nOrb,nOrb_twice,chem_pot+delta_chem_pot,trace_up,H_hfb,cp_tmp,R_tmp,eHFB_)
   call diag_H_hfb(nOrb,nOrb_twice,chem_pot-delta_chem_pot,trace_down,H_hfb,cp_tmp,R_tmp,eHFB_)
   call diag_H_hfb(nOrb,nOrb_twice,chem_pot-2d0*delta_chem_pot,trace_2down,H_hfb,cp_tmp,R_tmp,eHFB_)
   trace_1rdm =(trace_1rdm -nO_)**2d0
   trace_2up  =(trace_2up  -nO_)**2d0
   trace_up   =(trace_up   -nO_)**2d0
   trace_down =(trace_down -nO_)**2d0
   trace_2down=(trace_2down-nO_)**2d0
   grad_electrons = (-trace_2up+8d0*trace_up-8d0*trace_down+trace_2down)/(12d0*delta_chem_pot)
   chem_pot_change = -trace_1rdm/(grad_electrons+1d-10)
   ! Maximum change is bounded within +/- 0.10
   chem_pot_change = max( min( chem_pot_change , 0.1d0 / real(isteps) ), -0.1d0 / real(isteps) )
!   write(*,'(1X,A1,F16.10,1X,A1,F16.10,1X,A1F16.10,1X,A1)') &
!   '|',trace_1rdm,'|',chem_pot,'|',grad_electrons,'|'
  enddo
  write(*,*)'------------------------------------------------------'
  call diag_H_hfb(nOrb,nOrb_twice,chem_pot,trace_1rdm,H_hfb,cp,R,eHFB_)
  write(*,'(1X,A1,1X,A15,1X,A1,1X,A15,1X,A1A15,2X,A1)') &
          '|','Tr[1D]','|','Chem. Pot.','|','Grad N','|'
  write(*,*)'------------------------------------------------------'
  write(*,'(1X,A1,F16.10,1X,A1,F16.10,1X,A1F16.10,1X,A1)') &
  '|',trace_1rdm,'|',chem_pot,'|',grad_electrons,'|'
  write(*,*)'------------------------------------------------------'
  write(*,*)

  ! Reset H_HFB to its chemical potential version
 
  do iorb=1,nOrb
   H_hfb(iorb,iorb)=H_hfb(iorb,iorb)-chem_pot
   H_hfb(iorb+nOrb,iorb+nOrb)=H_hfb(iorb+nOrb,iorb+nOrb)+chem_pot
  enddo
  
  deallocate(R_tmp,cp_tmp)

end subroutine 

subroutine diag_H_hfb(nOrb,nOrb_twice,chem_pot,trace_1rdm,H_hfb,cp,R,eHFB_)

! Fix the chemical potential to integrate to the N for 2N electron systems

  implicit none

! Input variables

  integer,intent(in)            :: nOrb,nOrb_twice
  double precision,intent(in)   :: H_hfb(nOrb_twice,nOrb_twice) 
  double precision,intent(in)   :: chem_pot

! Local variables

  integer                       :: iorb

! Output variables

  double precision,intent(out)  :: trace_1rdm
  double precision,intent(inout):: cp(nOrb_twice,nOrb_twice) 
  double precision,intent(inout):: R(nOrb_twice,nOrb_twice) 
  double precision,intent(inout):: eHFB_(nOrb_twice) 
    
  cp(:,:) = H_hfb(:,:)
  do iorb=1,nOrb
   cp(iorb,iorb) = cp(iorb,iorb) - chem_pot 
   cp(iorb+nOrb,iorb+nOrb) = cp(iorb+nOrb,iorb+nOrb) + chem_pot 
  enddo

  ! Diagonalize H_HFB matrix

  call diagonalize_matrix(nOrb_twice,cp,eHFB_)
  
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

