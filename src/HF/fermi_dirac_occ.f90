subroutine fermi_dirac_occ(nO,nOrb,thrs_N,temperature,chem_pot,Occ,eHF)

! Use Fermi Dirac distribution to set up fractional Occs numbers and adjust the chemical potential to integrate to the N for 2N electron systems

  implicit none

! Input variables

  integer,intent(in)            :: nO,nOrb
  double precision,intent(in)   :: thrs_N
  double precision,intent(in)   :: temperature
  double precision,intent(in)   :: eHF(nOrb) 

! Local variables

  logical                       :: use_nelectrons
  integer                       :: iorb
  integer                       :: isteps
  double precision              :: nO_,nOp,nOm,nO2p,nO2m
  double precision              :: delta_chem_pot
  double precision              :: chem_pot_change
  double precision              :: grad_electrons
  double precision              :: trace_1rdm
  double precision              :: trace_up,trace_down
  double precision              :: trace_old
  double precision              :: thrs_closer

! Output variables

  double precision,intent(inout):: chem_pot
  double precision,intent(inout):: Occ(nOrb) 

  !  Initialize variables

  isteps = 0
  thrs_closer     = 2d-1
  delta_chem_pot  = 1.0d-1
  chem_pot_change = 0d0
  grad_electrons  = 1d0
  trace_1rdm      = -1d0
  trace_old       = 1d2
  nO_             = nO
  inquire(file='Nelectrons_RHFB', exist=use_nelectrons)
  if(use_nelectrons) then
    write(*,*) 'File Nelectrons_RHFB encountered, setting nO = nO_read/2'
    open(unit=314, form='formatted', file='Nelectrons_RHFB', status='old')
    read(314,*) nO_
    close(314)
    nO_=0.5d0*nO_
  endif

!  write(*,*)
!  write(*,*)' Fermi-Dirac distribution for the occ numbers'
!  write(*,*)
!  write(*,*)'-------------------------------------'
!  write(*,'(1X,A1,1X,A15,1X,A1,1X,A15,1X,A1)') &
!          '|','Tr[1D]','|','Chem. Pot.','|'
!  write(*,*)'-------------------------------------'

  ! First approach close the value with an error lower than 1

  Occ(:) = fermi_dirac(eHF,chem_pot,temperature)
  trace_old=sum(Occ(:))
!  write(*,'(1X,A1,F16.10,1X,A1,F16.10,1X,A1)') &
!  '|',trace_old,'|',chem_pot,'|'
  do while( abs(trace_old-nO_) > 1.0d0 .and. isteps <= 100 )
   isteps = isteps + 1
   trace_old =sum(fermi_dirac(eHF,chem_pot,temperature))
   trace_up  =sum(fermi_dirac(eHF,chem_pot+delta_chem_pot,temperature))
   trace_down=sum(fermi_dirac(eHF,chem_pot-delta_chem_pot,temperature))
!   write(*,'(1X,A1,F16.10,1X,A1,F16.10,1X,A1)') &
!   '|',trace_old,'|',chem_pot,'|'
   if( abs(trace_up-nO_) > abs(trace_old-nO_) .and. abs(trace_down-nO_) > abs(trace_old-nO_) ) then
!     write(*,'(1X,A1,F16.10,1X,A1,F16.10,1X,A1)') &
!     '|',trace_old,'|',chem_pot,'|'
     delta_chem_pot = 0.5d0*delta_chem_pot
     thrs_closer = 0.5d0*thrs_closer
!     write(*,*) "| contracting ...                   |"
     if(delta_chem_pot<1d-2) exit
   else
     if( abs(trace_up-nO_) < abs(trace_old-nO_) ) then
      chem_pot=chem_pot+delta_chem_pot
!      write(*,'(1X,A1,F16.10,1X,A1,F16.10,1X,A1)') &
!      '|',trace_up,'|',chem_pot,'|'
     else
      if( abs(trace_down-nO_) < abs(trace_old-nO_) ) then
       chem_pot=chem_pot-delta_chem_pot
!       write(*,'(1X,A1,F16.10,1X,A1,F16.10,1X,A1)') &
!       '|',trace_down,'|',chem_pot,'|'
      endif
     endif
   endif
  enddo

  ! Do  final search

!  write(*,*)'-------------------------------------'
!  write(*,'(1X,A1,1X,A15,1X,A1,1X,A15,1X,A1)') &
!          '|','Error Tr[1D]','|','Chem. Pot.','|'
!  write(*,*)'-------------------------------------'
  isteps=0
  delta_chem_pot = 1d-3
  trace_1rdm=(trace_1rdm-nO_)**2d0
  do while( sqrt(trace_1rdm) > thrs_N .and. isteps <= 100 )
    isteps = isteps + 1
    chem_pot = chem_pot + chem_pot_change
    trace_1rdm = (sum(fermi_dirac(eHF,chem_pot,temperature)) - nO_)**2d0 
    nOp  = (sum(fermi_dirac(eHF,chem_pot+delta_chem_pot,temperature)) - nO_)**2d0
    nOm  = (sum(fermi_dirac(eHF,chem_pot-delta_chem_pot,temperature)) - nO_)**2d0 
    nO2p = (sum(fermi_dirac(eHF,chem_pot+2d0*delta_chem_pot,temperature)) - nO_)**2d0 
    nO2m = (sum(fermi_dirac(eHF,chem_pot-2d0*delta_chem_pot,temperature)) - nO_)**2d0 
    grad_electrons = (-nO2p+8d0*nOp-8d0*nOm+nO2m)/(12d0*delta_chem_pot)
    chem_pot_change = -trace_1rdm/(grad_electrons+1d-10)
    ! Maximum change is bounded within +/- 0.10
    !chem_pot_change = max( min( chem_pot_change , 0.1d0 / real(isteps) ), -0.1d0 / real(isteps) )
!    write(*,'(1X,A1,F16.10,1X,A1,F16.10,1X,A1)') &
!    '|',trace_1rdm,'|',chem_pot,'|'
  enddo
  write(*,*)'-------------------------------------'
  Occ(:) = fermi_dirac(eHF,chem_pot,temperature)
  trace_1rdm=sum(Occ(:))
  write(*,'(1X,A1,1X,A15,1X,A1,1X,A15,1X,A1)') &
          '|','Tr[1D]','|','Chem. Pot.','|'
  write(*,*)'-------------------------------------'
  write(*,'(1X,A1,F16.10,1X,A1,F16.10,1X,A1)') &
  '|',trace_1rdm,'|',chem_pot,'|'
  write(*,*)'-------------------------------------'
  write(*,*)

  write(*,*)
  write(*,*) ' Initial occ. numbers '
  write(*,*)
  do iorb=1,nOrb
   write(*,'(3X,F16.10)') Occ(iorb)
  enddo

contains 

function fermi_dirac(eHF,chem_pot,temperature)
  implicit none
  double precision,intent(in) :: eHF(nOrb)
  double precision,intent(in) :: chem_pot
  double precision,intent(in) :: temperature
  double precision            :: fermi_dirac(nOrb)

  fermi_dirac(:) = 1d0 / ( 1d0 + exp((eHF(:) - chem_pot ) / temperature ) )

end function fermi_dirac

end subroutine 

