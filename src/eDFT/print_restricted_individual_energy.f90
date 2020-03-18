subroutine print_restricted_individual_energy(nEns,Ew,EwGOC,ET,EV,EJ,Ex,Ec,Exc,ExDD,EcDD,ExcDD,E, & 
                                              Om,Omx,Omc,Omxc,OmxDD,OmcDD,OmxcDD)

! Print individual energies for eDFT calculation

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)                 :: nEns
  double precision,intent(in)        :: Ew
  double precision,intent(in)        :: EwGOC
  double precision,intent(in)        :: ET(nEns)
  double precision,intent(in)        :: EV(nEns)
  double precision,intent(in)        :: EJ(nEns)
  double precision,intent(in)        :: Ex(nEns),   Ec(nEns),   Exc(nEns)
  double precision,intent(in)        :: ExDD(nEns), EcDD(nEns), ExcDD(nEns)
  double precision,intent(in)        :: Omx(nEns),  Omc(nEns),  Omxc(nEns)
  double precision,intent(in)        :: OmxDD(nEns),OmcDD(nEns),OmxcDD(nEns)
  double precision,intent(in)        :: E(nEns)
  double precision,intent(in)        :: Om(nEns)

! Local variables

  integer                            :: iEns

!------------------------------------------------------------------------
! Ensemble energies
!------------------------------------------------------------------------

  write(*,'(A60)')           '-------------------------------------------------'
  write(*,'(A60)')           ' ENSEMBLE ENERGIES'
  write(*,'(A60)')           '-------------------------------------------------'
  write(*,'(A44,F16.10,A3)') '     Ensemble energy:      ',Ew,   ' au'
  write(*,'(A44,F16.10,A3)') ' GOC Ensemble energy:      ',EwGOC,' au'
  write(*,'(A60)')           '-------------------------------------------------'
  write(*,*)

!------------------------------------------------------------------------
! Kinetic energy
!------------------------------------------------------------------------

  write(*,'(A60)')           '-------------------------------------------------'
  write(*,'(A60)')           ' INDIVIDUAL KINETIC     ENERGIES'
  write(*,'(A60)')           '-------------------------------------------------'
  do iEns=1,nEns
    write(*,'(A40,I2,A2,F16.10,A3)') ' Kinetic     energy state ',iEns,': ',ET(iEns),' au'
  end do
  write(*,'(A60)')           '-------------------------------------------------'
  write(*,*)

!------------------------------------------------------------------------
! Potential energy
!------------------------------------------------------------------------

  write(*,'(A60)')           '-------------------------------------------------'
  write(*,'(A60)')           ' INDIVIDUAL POTENTIAL   ENERGIES'
  write(*,'(A60)')           '-------------------------------------------------'
  do iEns=1,nEns
    write(*,'(A40,I2,A2,F16.10,A3)') ' Potential   energy state ',iEns,': ',EV(iEns),' au'
  end do
  write(*,'(A60)')           '-------------------------------------------------'
  write(*,*)

!------------------------------------------------------------------------
! Hartree energy
!------------------------------------------------------------------------

  write(*,'(A60)')           '-------------------------------------------------'
  write(*,'(A60)')           ' INDIVIDUAL HARTREE     ENERGIES'
  write(*,'(A60)')           '-------------------------------------------------'
  do iEns=1,nEns
    write(*,'(A40,I2,A2,F16.10,A3)') ' Hartree     energy state ',iEns,': ',EJ(iEns),' au'
  end do
  write(*,'(A60)')           '-------------------------------------------------'
  write(*,*)

!------------------------------------------------------------------------
! Exchange energy
!------------------------------------------------------------------------

  write(*,'(A60)')           '-------------------------------------------------'
  write(*,'(A60)')           ' INDIVIDUAL EXCHANGE    ENERGIES'
  write(*,'(A60)')           '-------------------------------------------------'
  do iEns=1,nEns
    write(*,'(A40,I2,A2,F16.10,A3)') ' Exchange    energy state ',iEns,': ',Ex(iEns),' au'
  end do
  write(*,'(A60)')           '-------------------------------------------------'
  write(*,*)

!------------------------------------------------------------------------
! Correlation energy
!------------------------------------------------------------------------

  write(*,'(A60)')           '-------------------------------------------------'
  write(*,'(A60)')           ' INDIVIDUAL CORRELATION ENERGIES'
  write(*,'(A60)')           '-------------------------------------------------'
  do iEns=1,nEns
    write(*,'(A40,I2,A2,F16.10,A3)') ' Correlation energy state ',iEns,': ',Ec(iEns),' au'
  end do
  write(*,'(A60)')           '-------------------------------------------------'
  write(*,*)

!------------------------------------------------------------------------
! Compute derivative discontinuities
!------------------------------------------------------------------------

  write(*,'(A60)')           '-------------------------------------------------'
  write(*,'(A60)')           ' ENSEMBLE DERIVATIVE CONTRIBUTIONS'
  write(*,'(A60)')           '-------------------------------------------------'
  do iEns=1,nEns
    write(*,*)
    write(*,'(A40,I2,A2,F16.10,A3)') '  x ensemble derivative state ',iEns,': ',ExDD(iEns), ' au'
    write(*,'(A40,I2,A2,F16.10,A3)') '  c ensemble derivative state ',iEns,': ',EcDD(iEns), ' au'
    write(*,'(A40,I2,A2,F16.10,A3)') ' xc ensemble derivative state ',iEns,': ',ExcDD(iEns),' au'
  end do
  write(*,'(A60)')           '-------------------------------------------------'
  write(*,*)

!------------------------------------------------------------------------
! Total and Excitation energies
!------------------------------------------------------------------------

  write(*,'(A60)')           '-------------------------------------------------'
  write(*,'(A60)')           ' INDIVIDUAL AND EXCITATION ENERGIES'
  write(*,'(A60)')           '-------------------------------------------------'
  do iEns=1,nEns
    write(*,'(A40,I2,A2,F16.10,A3)') ' Individual energy state ',iEns,': ',E(iEns),' au'
  end do
  write(*,'(A60)')           '-------------------------------------------------'

  do iEns=2,nEns
    write(*,'(A40,I2,A2,F16.10,A3)') ' Excitation energy  1 ->',iEns,': ',Om(iEns),    ' au'
    write(*,*)
    write(*,'(A44,      F16.10,A3)') '     x  energy contribution     : ',Omx(iEns),   ' au'
    write(*,'(A44,      F16.10,A3)') '     c  energy contribution     : ',Omc(iEns),   ' au'
    write(*,'(A44,      F16.10,A3)') '    xc  energy contribution     : ',Omxc(iEns),  ' au'
    write(*,*)
    write(*,'(A44,      F16.10,A3)') '     x  ensemble derivative     : ',OmxDD(iEns), ' au'
    write(*,'(A44,      F16.10,A3)') '     c  ensemble derivative     : ',OmcDD(iEns), ' au'
    write(*,'(A44,      F16.10,A3)') '    xc  ensemble derivative     : ',OmxcDD(iEns),' au'
  end do
  write(*,'(A60)')           '-------------------------------------------------'

  do iEns=2,nEns
    write(*,'(A40,I2,A2,F16.10,A3)') ' Excitation energy  1 ->',iEns,': ',Om(iEns)*HaToeV,    ' eV'
    write(*,*)
    write(*,'(A44,      F16.10,A3)') '     x  energy contribution     : ',Omx(iEns)*HaToeV,   ' eV'
    write(*,'(A44,      F16.10,A3)') '     c  energy contribution     : ',Omc(iEns)*HaToeV,   ' eV'
    write(*,'(A44,      F16.10,A3)') '    xc  energy contribution     : ',Omxc(iEns)*HaToeV,  ' eV'
    write(*,*)
    write(*,'(A44,      F16.10,A3)') '     x  ensemble derivative     : ',OmxDD(iEns)*HaToeV, ' eV'
    write(*,'(A44,      F16.10,A3)') '     c  ensemble derivative     : ',OmcDD(iEns)*HaToeV, ' eV'
    write(*,'(A44,      F16.10,A3)') '    xc  ensemble derivative     : ',OmxcDD(iEns)*HaToeV,' eV'
  end do
  write(*,'(A60)')           '-------------------------------------------------'
  write(*,*)

end subroutine print_restricted_individual_energy
