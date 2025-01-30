subroutine fermi_dirac_occ(nO,nOrb,thrs_N,temperature,chem_pot,Occ,eHF)

! Use Fermi Dirac distribution to set up fractional Occs numbers and adjust the chemical potential to integrate to the N for 2N electron systems

  implicit none

! Input variables

  integer,intent(in)            :: nO,nOrb
  double precision,intent(in)   :: thrs_N
  double precision,intent(in)   :: temperature

! Local variables

  integer                       :: iorb
  integer                       :: isteps
  double precision              :: delta_chem_pot
  double precision              :: chem_pot_change
  double precision              :: grad_electrons
  double precision              :: trace_1rdm

! Output variables

  double precision,intent(inout):: chem_pot
  double precision,intent(inout):: eHF(nOrb) 
  double precision,intent(inout):: Occ(nOrb) 

  !  Initialize variables

  isteps=0
  delta_chem_pot = 1.0d-3
  trace_1rdm     = -1.0d0
  chem_pot_change = 0.0d0

  write(*,*)
  write(*,*)' Fermi-Dirac distribution for the occ numbers'
  write(*,*)
  write(*,*)'-------------------------------------'
  write(*,'(1X,A1,1X,A15,1X,A1,1X,A15,1X,A1)') &
          '|','Tr[1D]','|','Chem. Pot.','|'
  write(*,*)'-------------------------------------'


  Occ(:) = fermi_dirac(eHF,chem_pot,temperature)
  trace_1rdm=sum(Occ(:))
  write(*,'(1X,A1,F16.10,1X,A1,F16.10,1X,A1)') &
  '|',trace_1rdm,'|',chem_pot,'|'


  do while( abs(trace_1rdm-nO) > thrs_N .and. isteps <= 100 )
    isteps = isteps + 1
    chem_pot = chem_pot + chem_pot_change
    Occ(:) = fermi_dirac(eHF,chem_pot,temperature)
    trace_1rdm    = sum(Occ(:))
    grad_electrons = ( sum(fermi_dirac(eHF,chem_pot+delta_chem_pot,temperature)) &
                   - sum(fermi_dirac(eHF,chem_pot-delta_chem_pot,temperature)) )/(2.0d0*delta_chem_pot)
    chem_pot_change = -(trace_1rdm-nO)/(grad_electrons+1d-10)
    write(*,'(1X,A1,F16.10,1X,A1,F16.10,1X,A1)') &
    '|',trace_1rdm,'|',chem_pot,'|'
  enddo
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

