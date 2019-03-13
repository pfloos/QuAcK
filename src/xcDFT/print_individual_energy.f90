subroutine print_individual_energy(nEns,EJ,Ex,Ec,EcLZ,EcDD,E,Om)

! Print individual energies for eDFT calculation

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)                 :: nEns
  double precision,intent(in)        :: EJ(nsp,nEns)
  double precision,intent(in)        :: Ex(nspin,nEns)
  double precision,intent(in)        :: Ec(nsp,nEns)
  double precision,intent(in)        :: EcLZ(nsp)
  double precision,intent(in)        :: EcDD(nsp,nEns)
  double precision,intent(in)        :: E(nEns)
  double precision,intent(in)        :: Om(nEns)

! Local variables

  integer                            :: iEns

!------------------------------------------------------------------------
! Hartree energy
!------------------------------------------------------------------------

  write(*,'(A60)')           '-------------------------------------------------'
  write(*,'(A50)')           ' Individual Hartree     energies'
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
  write(*,'(A50)')           ' Individual exchange    energies'
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
  write(*,'(A50)')           ' Individual correlation energies'
  write(*,'(A60)')           '-------------------------------------------------'
  do iEns=1,nEns
    write(*,'(A40,I2,A2,F16.10,A3)') ' Correlation energy state ',iEns,': ',sum(Ec(:,iEns)),' au'
  end do
  write(*,'(A60)')           '-------------------------------------------------'
  write(*,*)

!------------------------------------------------------------------------
! Compute Levy-Zahariev shift
!------------------------------------------------------------------------

  write(*,'(A60)')              '-------------------------------------------------'
  write(*,'(A40,2X,2X,F16.10,A3)') ' Levy-Zahariev shifts: ',sum(EcLZ(:)),' au'
  write(*,'(A60)')              '-------------------------------------------------'
  write(*,*)

!------------------------------------------------------------------------
! Compute derivative discontinuities
!------------------------------------------------------------------------

  write(*,'(A60)')           '-------------------------------------------------'
  write(*,'(A50)')           ' Derivative discontinuities (DD)    '
  write(*,'(A60)')           '-------------------------------------------------'
  do iEns=1,nEns
    write(*,'(A40,I2,A2,F16.10,A3)') ' Correlation part of DD ',iEns,': ',sum(EcDD(:,iEns)),' au'
  end do
  write(*,'(A60)')           '-------------------------------------------------'
  write(*,*)

!------------------------------------------------------------------------
! Total and Excitation energies
!------------------------------------------------------------------------

  write(*,'(A60)')           '-------------------------------------------------'
  write(*,'(A50)')           ' Individual and excitation energies '
  write(*,'(A60)')           '-------------------------------------------------'
  do iEns=1,nEns
    write(*,'(A40,I2,A2,F16.10,A3)') ' Individual energy state ',iEns,': ',E(iEns),' au'
  end do
  write(*,'(A60)')           '-------------------------------------------------'
  do iEns=2,nEns
    write(*,'(A40,I2,A2,F16.10,A3)') ' Excitation energy  1 ->',iEns,': ',Om(iEns),' au'
  end do
  write(*,'(A60)')           '-------------------------------------------------'
  do iEns=2,nEns
    write(*,'(A40,I2,A2,F16.10,A3)') ' Excitation energy  1 ->',iEns,': ',Om(iEns)*HaToeV,' eV'
  end do
  write(*,'(A60)')           '-------------------------------------------------'
  write(*,*)

end subroutine print_individual_energy
