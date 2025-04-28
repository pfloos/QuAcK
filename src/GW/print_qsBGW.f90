
! ---

subroutine print_qsBGW(nBas, nOrb, nOrb2, nO, N_anom, Occ, eqsGW_state, ENuc, ET, EV, EJ, EK, EL, EqsGW, chem_pot, dipole, Delta_HL)

! Print one-electron energies and other stuff

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)                 :: nBas, nOrb, nOrb2
  integer,intent(in)                 :: nO
  double precision,intent(in)        :: Occ(nOrb)
  double precision,intent(in)        :: eqsGW_state(nOrb2)
  double precision,intent(in)        :: ENuc
  double precision,intent(in)        :: ET
  double precision,intent(in)        :: EV
  double precision,intent(in)        :: EJ
  double precision,intent(in)        :: EK
  double precision,intent(in)        :: EL
  double precision,intent(in)        :: EqsGW
  double precision,intent(in)        :: chem_pot
  double precision,intent(in)        :: N_anom
  double precision,intent(in)        :: Delta_HL
  double precision,intent(in)        :: dipole(ncart)

! Local variables

  integer                            :: iorb
  integer                            :: ixyz
  double precision                   :: trace_occ

  logical                            :: dump_orb = .false.

! Dump results

  write(*,*)
  write(*,'(A50)')           '---------------------------------------'
  write(*,'(A33)')           ' Summary               '
  write(*,'(A50)')           '---------------------------------------'
  write(*,'(A33,1X,F16.10,A3)') ' One-electron energy = ',ET + EV,' au'
  write(*,'(A33,1X,F16.10,A3)') ' Kinetic      energy = ',ET,' au'
  write(*,'(A33,1X,F16.10,A3)') ' Potential    energy = ',EV,' au'
  write(*,'(A50)')           '---------------------------------------'
  write(*,'(A33,1X,F16.10,A3)') ' Two-electron energy = ',EJ + EK,' au'
  write(*,'(A33,1X,F16.10,A3)') ' Hartree      energy = ',EJ,' au'
  write(*,'(A33,1X,F16.10,A3)') ' Exchange     energy = ',EK,' au'
  write(*,'(A33,1X,F16.10,A3)') ' Anomalous    energy = ',EL,' au'
  write(*,'(A50)')           '---------------------------------------'
  write(*,'(A33,1X,F16.10,A3)') ' Electronic   energy = ',EqsGW,' au'
  write(*,'(A33,1X,F16.10,A3)') ' Nuclear   repulsion = ',ENuc,' au'
  write(*,'(A33,1X,F16.10,A3)') ' qsGW          energy = ',EqsGW + ENuc,' au'
  write(*,'(A50)')           '---------------------------------------'
  write(*,'(A33,1X,F16.10,A3)') ' Chemical potential  = ',chem_pot,' au'
  write(*,'(A33,1X,F16.10,A3)') ' | Anomalous dens |  = ',N_anom,'   '
  write(*,'(A33,1X,F16.10,A3)') ' Delta QP HOMO-LUMO  = ',Delta_HL,' au'
  write(*,'(A50)')           '---------------------------------------'
  write(*,'(A36)')           ' Dipole moment (Debye)    '
  write(*,'(10X,4A10)')      'X','Y','Z','Tot.'
  write(*,'(10X,4F10.4)')    (dipole(ixyz)*auToD,ixyz=1,ncart),norm2(dipole)*auToD
  write(*,'(A50)')           '---------------------------------------'
  write(*,*)

! Print results

  write(*,'(A50)') '---------------------------------------'
  write(*,'(A50)') ' qsGW occupation numbers '
  write(*,'(A50)') '---------------------------------------'
  trace_occ=0d0
  do iorb=1,nOrb
   if(abs(Occ(iorb))>1d-8) then
    write(*,'(I7,10F15.8)') iorb,2d0*Occ(iorb)
   endif
   trace_occ=trace_occ+2d0*Occ(iorb)
  enddo
  write(*,*)
  write(*,'(A33,1X,F16.10,A3)') ' Trace [ 1D ]        = ',trace_occ,'   '
  write(*,*)

  write(*,'(A50)') '---------------------------------------'
  write(*,'(A50)') ' qsGW QP energies '
  write(*,'(A50)') '---------------------------------------'
  do iorb=1,nOrb2
   write(*,'(I7,10F15.8)') iorb,eqsGW_state(iorb)
  enddo
  write(*,*)


end subroutine 
