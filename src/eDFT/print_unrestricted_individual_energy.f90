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
! Individual energies 
!------------------------------------------------------------------------

  write(*,'(A60)')           '-------------------------------------------------'
  write(*,'(A60)')           ' INDIVIDUAL TOTAL       ENERGIES'
  write(*,'(A60)')           '-------------------------------------------------'
  do iEns=1,nEns
    write(*,'(A40,I2,A2,F16.10,A3)') ' Individual energy state ',iEns,': ',E(iEns) + ENuc,' au'
  end do
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

!  write(*,'(A60)') '-------------------------------------------------'
!  write(*,'(A60)') ' IP AND EA FROM AUXILIARY ENERGIES '
!  write(*,'(A60)') '-------------------------------------------------'

!  do iEns=2,nEns
!    write(*,'(A40,I2,A1,F16.10,A3)') ' Energy difference 1 -> ',iEns,':',Omaux(iEns)+OmxcDD(iEns),' au'
!    write(*,*)
!    write(*,'(A44,      F16.10,A3)') ' auxiliary energy contribution  : ',Omaux(iEns), ' au'
!    write(*,'(A44,      F16.10,A3)') '        x  ensemble derivative  : ',OmxDD(iEns), ' au'
!    write(*,'(A44,      F16.10,A3)') '        c  ensemble derivative  : ',OmcDD(iEns), ' au'
!    write(*,'(A44,      F16.10,A3)') '       xc  ensemble derivative  : ',OmxcDD(iEns),' au'
!    write(*,*)

!    write(*,'(A60)') '-------------------------------------------------'
!    write(*,*)

!    write(*,'(A40,I2,A1,F16.10,A3)') ' Energy difference 1 -> ',iEns,':',(Omaux(iEns)+OmxcDD(iEns))*HaToeV,' eV'
!    write(*,*)
!    write(*,'(A44,      F16.10,A3)') ' auxiliary energy contribution  : ',Omaux(iEns)*HaToeV, ' eV'
!    write(*,'(A44,      F16.10,A3)') '        x  ensemble derivative  : ',OmxDD(iEns)*HaToeV, ' eV'
!    write(*,'(A44,      F16.10,A3)') '        c  ensemble derivative  : ',OmcDD(iEns)*HaToeV, ' eV'
!    write(*,'(A44,      F16.10,A3)') '       xc  ensemble derivative  : ',OmxcDD(iEns)*HaToeV,' eV'
!    write(*,*)
!  end do

!  write(*,'(A60)') '-------------------------------------------------'
!  write(*,*)

 write(*,'(A60)') '-------------------------------------------------'
  write(*,'(A60)') ' IP and EA FROM INDIVIDUAL ENERGIES '
  write(*,'(A60)') '-------------------------------------------------'
  do iEns=1,nEns
!    write(*,'(A40,I2,A2,F16.10,A3)') ' Individual energy state ',iEns,': ',E(iEns) + ENuc,' au'
  end do
  write(*,'(A60)') '-------------------------------------------------'

  do iEns=2,nEns
    write(*,'(A40,I2,A1,F16.10,A3)') ' Energy difference 1 -> ',iEns,':',Om(iEns),    ' au'
    write(*,*)
    write(*,'(A44,      F16.10,A3)') '     x  energy contribution     : ',Omx(iEns),   ' au'
    write(*,'(A44,      F16.10,A3)') '     c  energy contribution     : ',Omc(iEns),   ' au'
    write(*,'(A44,      F16.10,A3)') '    xc  energy contribution     : ',Omxc(iEns),  ' au'
    write(*,*)
    write(*,'(A44,      F16.10,A3)') '     x  ensemble derivative     : ',OmxDD(iEns), ' au'
    write(*,'(A44,      F16.10,A3)') '     c  ensemble derivative     : ',OmcDD(iEns), ' au'
    write(*,'(A44,      F16.10,A3)') '    xc  ensemble derivative     : ',OmxcDD(iENs),' au'
    write(*,*)
   
    write(*,'(A60)') '-------------------------------------------------'
   
    write(*,'(A40,I2,A1,F16.10,A3)') ' Energy difference 1 -> ',iEns,':',Om(iEns)*HaToeV,    ' eV'
    write(*,*)
    write(*,'(A44,      F16.10,A3)') '     x  energy contribution     : ',Omx(iEns)*HaToeV,   ' eV'
    write(*,'(A44,      F16.10,A3)') '     c  energy contribution     : ',Omc(iEns)*HaToeV,   ' eV'
    write(*,'(A44,      F16.10,A3)') '    xc  energy contribution     : ',Omxc(iEns)*HaToeV,  ' eV'
    write(*,*)
    write(*,'(A44,      F16.10,A3)') '     x  ensemble derivative     : ',OmxDD(iEns)*HaToeV, ' eV'
    write(*,'(A44,      F16.10,A3)') '     c  ensemble derivative     : ',OmcDD(iEns)*HaToeV, ' eV'
    write(*,'(A44,      F16.10,A3)') '    xc  ensemble derivative     : ',OmxcDD(iEns)*HaToeV,' eV'
   write(*,*)
  end do 
  write(*,'(A60)') '-------------------------------------------------'
  write(*,*)

end subroutine print_unrestricted_individual_energy
