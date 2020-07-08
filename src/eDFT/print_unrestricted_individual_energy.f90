subroutine print_unrestricted_individual_energy(nEns,ENuc,Ew,ET,EV,EJ,Ex,Ec,Exc,Eaux,ExDD,EcDD,ExcDD,E, & 
                                                Om,Omx,Omc,Omxc,Omaux,OmxDD,OmcDD,OmxcDD)

! Print individual energies for eDFT calculation

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)                 :: nEns
  double precision,intent(in)        :: ENuc
  double precision,intent(in)        :: Ew
  double precision,intent(in)        :: ET(nspin,nEns)
  double precision,intent(in)        :: EV(nspin,nEns)
  double precision,intent(in)        :: EJ(nsp,nEns)
  double precision,intent(in)        :: Ex(nspin,nEns),   Ec(nsp,nEns),   Exc(nEns)
  double precision,intent(in)        :: Eaux(nspin,nEns)
  double precision,intent(in)        :: ExDD(nspin,nEns), EcDD(nsp,nEns), ExcDD(nsp,nEns)
  double precision,intent(in)        :: Omx(nEns),  Omc(nEns),  Omxc(nEns)
  double precision,intent(in)        :: Omaux(nEns)
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
  write(*,'(A44,F16.10,A3)') '     Ensemble energy:      ',Ew    + ENuc,' au'
  write(*,'(A60)')           '-------------------------------------------------'
  write(*,*)

!------------------------------------------------------------------------
! Kinetic energy
!------------------------------------------------------------------------

  write(*,'(A60)')           '-------------------------------------------------'
  write(*,'(A60)')           ' INDIVIDUAL KINETIC     ENERGIES'
  write(*,'(A60)')           '-------------------------------------------------'
  do iEns=1,nEns
    write(*,'(A40,I2,A2,F16.10,A3)') ' Kinetic     energy state ',iEns,': ',sum(ET(:,iEns)),' au'
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
    write(*,'(A40,I2,A2,F16.10,A3)') ' Potential   energy state ',iEns,': ',sum(EV(:,iEns)),' au'
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
    write(*,'(A40,I2,A2,F16.10,A3)') ' Hartree     energy state ',iEns,': ',sum(EJ(:,iEns)),' au'
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
    write(*,'(A40,I2,A2,F16.10,A3)') ' Exchange    energy state ',iEns,': ',sum(Ex(:,iEns)),' au'
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
    write(*,'(A40,I2,A2,F16.10,A3)') ' Correlation energy state ',iEns,': ',sum(Ec(:,iEns)),' au'
  end do
  write(*,'(A60)')           '-------------------------------------------------'
  write(*,*)

!------------------------------------------------------------------------
! Auxiliary energies
!------------------------------------------------------------------------

  write(*,'(A60)')           '-------------------------------------------------'
  write(*,'(A60)')           ' AUXILIARY KS ENERGIES'
  write(*,'(A60)')           '-------------------------------------------------'
  do iEns=1,nEns
    write(*,'(A40,I2,A2,F16.10,A3)') 'Auxiliary KS energy state ',iEns,': ',sum(Eaux(:,iEns)),' au'
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
    write(*,'(A40,I2,A2,F16.10,A3)') '  x ensemble derivative state ',iEns,': ',sum(ExDD(:,iEns)), ' au'
    write(*,'(A40,I2,A2,F16.10,A3)') '  c ensemble derivative state ',iEns,': ',sum(EcDD(:,iEns)), ' au'
    write(*,'(A40,I2,A2,F16.10,A3)') ' xc ensemble derivative state ',iEns,': ',sum(ExcDD(:,iEns)),' au'
  end do
  write(*,'(A60)')           '-------------------------------------------------'
  write(*,*)

!------------------------------------------------------------------------
! Total Energy and IP and EA
!------------------------------------------------------------------------

  write(*,'(A60)') '-------------------------------------------------'
  write(*,'(A60)') ' IP AND EA FROM AUXILIARY ENERGIES '
  write(*,'(A60)') '-------------------------------------------------'

    write(*,'(A43,F16.10,A4)') ' Ionization Potential  1 -> 2:',Omaux(2)+OmxcDD(2),' au'
    write(*,*)
    write(*,'(A44,      F16.10,A3)') ' auxiliary energy contribution  : ',Omaux(2), ' au'
    write(*,'(A44,      F16.10,A3)') '        x  ensemble derivative  : ',OmxDD(2), ' au'
    write(*,'(A44,      F16.10,A3)') '        c  ensemble derivative  : ',OmcDD(2), ' au'
    write(*,'(A44,      F16.10,A3)') '       xc  ensemble derivative  : ',OmxcDD(2),' au'
    write(*,*)
    write(*,'(A43,F16.10,A4)') ' Electronic Affinity  1 -> 3:',Omaux(3)+OmxcDD(3),' au'
    write(*,*)
    write(*,'(A44,      F16.10,A3)') ' auxiliary energy contribution  : ',Omaux(3), ' au'
    write(*,'(A44,      F16.10,A3)') '        x  ensemble derivative  : ',OmxDD(3), ' au'
    write(*,'(A44,      F16.10,A3)') '        c  ensemble derivative  : ',OmcDD(3), ' au'
    write(*,'(A44,      F16.10,A3)') '       xc  ensemble derivative  : ',OmxcDD(3),' au'
    write(*,*)

  write(*,'(A60)') '-------------------------------------------------'
  write(*,*)

    write(*,'(A40,F16.10,A3)') ' Ionization Potential  1 -> 2:',(Omaux(2)+OmxcDD(2))*HaToeV,' eV'
    write(*,*)
    write(*,'(A44,      F16.10,A3)') ' auxiliary energy contribution  : ',Omaux(2)*HaToeV, ' eV'
    write(*,'(A44,      F16.10,A3)') '        x  ensemble derivative  : ',OmxDD(2)*HaToeV, ' eV'
    write(*,'(A44,      F16.10,A3)') '        c  ensemble derivative  : ',OmcDD(2)*HaToeV, ' eV'
    write(*,'(A44,      F16.10,A3)') '       xc  ensemble derivative  : ',OmxcDD(2)*HaToeV,' eV'
    write(*,*)
    write(*,'(A40,F16.10,A3)') ' Electronic Affinity  1 -> 3:',(Omaux(3)+OmxcDD(3))*HaToeV,' eV'
    write(*,*)
    write(*,'(A44,      F16.10,A3)') ' auxiliary energy contribution  : ',Omaux(3)*HaToeV, ' eV'
    write(*,'(A44,      F16.10,A3)') '        x  ensemble derivative  : ',OmxDD(3)*HaToeV, ' eV'
    write(*,'(A44,      F16.10,A3)') '        c  ensemble derivative  : ',OmcDD(3)*HaToeV, ' eV'
    write(*,'(A44,      F16.10,A3)') '       xc  ensemble derivative  : ',OmxcDD(3)*HaToeV,' eV'
    write(*,*)

  write(*,'(A60)') '-------------------------------------------------'
  write(*,*)

  write(*,'(A60)') '-------------------------------------------------'
  write(*,'(A60)') ' IP and EA FROM INDIVIDUAL ENERGIES '
  write(*,'(A60)') '-------------------------------------------------'
  do iEns=1,nEns
    write(*,'(A40,I2,A2,F16.10,A3)') ' Individual energy state ',iEns,': ',E(iEns) + ENuc,' au'
  end do
  write(*,'(A60)') '-------------------------------------------------'

    write(*,'(A43,F16.10,A4)') ' Ionization Potential  1 -> 2:',Om(2),    ' au'
    write(*,*)
    write(*,'(A44,      F16.10,A3)') '     x  energy contribution     : ',Omx(2),   ' au'
    write(*,'(A44,      F16.10,A3)') '     c  energy contribution     : ',Omc(2),   ' au'
    write(*,'(A44,      F16.10,A3)') '    xc  energy contribution     : ',Omxc(2),  ' au'
    write(*,*)
    write(*,'(A44,      F16.10,A3)') '     x  ensemble derivative     : ',OmxDD(2), ' au'
    write(*,'(A44,      F16.10,A3)') '     c  ensemble derivative     : ',OmcDD(2), ' au'
    write(*,'(A44,      F16.10,A3)') '    xc  ensemble derivative     : ',OmxcDD(2),' au'
    write(*,*)
    write(*,'(A43,F16.10,A4)') ' Electronic Affinity  1 -> 3:',Om(3),    ' au'
    write(*,*)
    write(*,'(A44,      F16.10,A3)') '     x  energy contribution     : ',Omx(3),   ' au'
    write(*,'(A44,      F16.10,A3)') '     c  energy contribution     : ',Omc(3),   ' au'
    write(*,'(A44,      F16.10,A3)') '    xc  energy contribution     : ',Omxc(3),  ' au'
    write(*,*)
    write(*,'(A44,      F16.10,A3)') '     x  ensemble derivative     : ',OmxDD(3), ' au'
    write(*,'(A44,      F16.10,A3)') '     c  ensemble derivative     : ',OmcDD(3), ' au'
    write(*,'(A44,      F16.10,A3)') '    xc  ensemble derivative     : ',OmxcDD(3),' au'
    write(*,*)

    write(*,'(A60)') '-------------------------------------------------'

    write(*,'(A43,F16.10,A4)') ' Ionization Potential 1 -> 2:',Om(2)*HaToeV,    ' eV'
    write(*,*)
    write(*,'(A44,      F16.10,A3)') '     x  energy contribution     : ',Omx(2)*HaToeV,   ' eV'
    write(*,'(A44,      F16.10,A3)') '     c  energy contribution     : ',Omc(2)*HaToeV,   ' eV'
    write(*,'(A44,      F16.10,A3)') '    xc  energy contribution     : ',Omxc(2)*HaToeV,  ' eV'
    write(*,*)
    write(*,'(A44,      F16.10,A3)') '     x  ensemble derivative     : ',OmxDD(2)*HaToeV, ' eV'
    write(*,'(A44,      F16.10,A3)') '     c  ensemble derivative     : ',OmcDD(2)*HaToeV, ' eV'
    write(*,'(A44,      F16.10,A3)') '    xc  ensemble derivative     : ',OmxcDD(2)*HaToeV,' eV'
    write(*,*)
    write(*,'(A43,F16.10,A4)') ' Electronic Affinity 1 -> 3:',Om(3)*HaToeV,    ' eV'
    write(*,*)
    write(*,'(A44,      F16.10,A3)') '     x  energy contribution     : ',Omx(3)*HaToeV,   ' eV'
    write(*,'(A44,      F16.10,A3)') '     c  energy contribution     : ',Omc(3)*HaToeV,   ' eV'
    write(*,'(A44,      F16.10,A3)') '    xc  energy contribution     : ',Omxc(3)*HaToeV,  ' eV'
    write(*,*)
    write(*,'(A44,      F16.10,A3)') '     x  ensemble derivative     : ',OmxDD(3)*HaToeV, ' eV'
    write(*,'(A44,      F16.10,A3)') '     c  ensemble derivative     : ',OmcDD(3)*HaToeV, ' eV'
    write(*,'(A44,      F16.10,A3)') '    xc  ensemble derivative     : ',OmxcDD(3)*HaToeV,' eV'
    write(*,*)

  write(*,'(A60)') '-------------------------------------------------'
 
  write(*,*)

end subroutine print_unrestricted_individual_energy
