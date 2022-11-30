subroutine MP2(regularize,nBas,nC,nO,nV,nR,ERI,ENuc,EHF,e,EcMP2)

! Perform second-order Moller-Plesset calculation with and without regularizers

  implicit none

! Input variables

  logical,intent(in)            :: regularize
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: EHF
  double precision,intent(in)   :: e(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: i,j,a,b

  double precision              :: kappa,sigm1,sigm2
  double precision              :: Dijab
  double precision              :: fs,fs2,fk

  double precision              :: E2d,E2ds,E2ds2,E2dk
  double precision              :: E2x,E2xs,E2xs2,E2xk

  double precision              :: EcsMP2,Ecs2MP2,EckMP2

! Output variables

  double precision,intent(out)  :: EcMP2

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'|    Moller-Plesset second-order calculation   |'
  write(*,*)'************************************************'
  write(*,*)

!---------------------------------------------!
! Parameters for regularized MP2 calculations !
!---------------------------------------------!

  kappa = 1.1d0
  sigm1 = 0.7d0
  sigm2 = 0.4d0

!--------------------------------------------------!
! Compute conventinal and regularized MP2 energies !
!--------------------------------------------------!

  E2d   = 0d0
  E2ds  = 0d0
  E2ds2 = 0d0
  E2dk  = 0d0

  E2x   = 0d0
  E2xs  = 0d0
  E2xs2 = 0d0
  E2xk  = 0d0

  do i=nC+1,nO
    do j=nC+1,nO
      do a=nO+1,nBas-nR
        do b=nO+1,nBas-nR

          Dijab = e(a) + e(b) - e(i) - e(j) 

          ! Second-order ring diagram

          fs  = (1d0 - exp(-sigm1*Dijab))/Dijab
          fs2 = (1d0 - exp(-sigm2*Dijab*Dijab))/Dijab
          fk  = (1d0 - exp(-kappa*Dijab))**2/Dijab

          E2d   = E2d   - ERI(i,j,a,b)*ERI(i,j,a,b)/Dijab
          E2ds  = E2ds  - ERI(i,j,a,b)*ERI(i,j,a,b)*fs 
          E2ds2 = E2ds2 - ERI(i,j,a,b)*ERI(i,j,a,b)*fs2
          E2dk  = E2dk  - ERI(i,j,a,b)*ERI(i,j,a,b)*fk

          ! Second-order exchange diagram

          E2x   = E2x   - ERI(i,j,a,b)*ERI(i,j,b,a)/Dijab
          E2xs  = E2xs  - ERI(i,j,a,b)*ERI(i,j,b,a)*fs 
          E2xs2 = E2xs2 - ERI(i,j,a,b)*ERI(i,j,b,a)*fs2
          E2xk  = E2xk  - ERI(i,j,a,b)*ERI(i,j,b,a)*fk

        enddo
      enddo
    enddo
  enddo

  EcMP2   = 2d0*E2d   - E2x
  EcsMP2  = 2d0*E2ds  - E2xs
  Ecs2MP2 = 2d0*E2ds2 - E2xs2
  EckMP2  = 2d0*E2dk  - E2xk

!------------!
! MP2 energy !
!------------!

  write(*,*)
  write(*,'(A32)')           '--------------------------'
  write(*,'(A32)')           ' MP2 calculation          '
  write(*,'(A32)')           '--------------------------'
  write(*,'(A32,1X,F16.10)') ' MP2 correlation energy = ',EcMP2
  write(*,'(A32,1X,F16.10)') ' Direct part            = ',2d0*E2d
  write(*,'(A32,1X,F16.10)') ' Exchange part          = ',-E2x
  write(*,'(A32)')           '--------------------------'
  write(*,'(A32,1X,F16.10)') ' MP2 electronic  energy = ',EHF + EcMP2
  write(*,'(A32,1X,F16.10)') ' MP2 total       energy = ',ENuc + EHF + EcMP2
  write(*,'(A32)')           '--------------------------'
  write(*,*)

  if(regularize) then

!-------------------!
! sigma1-MP2 energy !
!-------------------!

    write(*,*)
    write(*,'(A32)')           '--------------------------'
    write(*,'(A32)')           ' sigma-MP2 calculation    '
    write(*,'(A32)')           '--------------------------'
    write(*,'(A32,1X,F16.10)') ' MP2 correlation energy = ',EcsMP2
    write(*,'(A32,1X,F16.10)') ' Direct part            = ',2d0*E2ds
    write(*,'(A32,1X,F16.10)') ' Exchange part          = ',-E2xs
    write(*,'(A32)')           '--------------------------'
    write(*,'(A32,1X,F16.10)') ' MP2 electronic  energy = ',EHF + EcsMP2
    write(*,'(A32,1X,F16.10)') ' MP2 total       energy = ',ENuc + EHF + EcsMP2
    write(*,'(A32)')           '--------------------------'
    write(*,*)

!--------------------!
! sigma^2-MP2 energy !
!--------------------!

    write(*,*)
    write(*,'(A32)')           '--------------------------'
    write(*,'(A32)')           ' sigma^2-MP2 calculation  '
    write(*,'(A32)')           '--------------------------'
    write(*,'(A32,1X,F16.10)') ' MP2 correlation energy = ',Ecs2MP2
    write(*,'(A32,1X,F16.10)') ' Direct part            = ',2d0*E2ds2
    write(*,'(A32,1X,F16.10)') ' Exchange part          = ',-E2xs2
    write(*,'(A32)')           '--------------------------'
    write(*,'(A32,1X,F16.10)') ' MP2 electronic  energy = ',EHF + Ecs2MP2
    write(*,'(A32,1X,F16.10)') ' MP2 total       energy = ',ENuc + EHF + Ecs2MP2
    write(*,'(A32)')           '--------------------------'
    write(*,*)

!------------------!
! kappa-MP2 energy !
!------------------!

    write(*,*)
    write(*,'(A32)')           '--------------------------'
    write(*,'(A32)')           ' kappa-MP2 calculation    '
    write(*,'(A32)')           '--------------------------'
    write(*,'(A32,1X,F16.10)') ' MP2 correlation energy = ',EckMP2
    write(*,'(A32,1X,F16.10)') ' Direct part            = ',2d0*E2dk
    write(*,'(A32,1X,F16.10)') ' Exchange part          = ',-E2xk
    write(*,'(A32)')           '--------------------------'
    write(*,'(A32,1X,F16.10)') ' MP2 electronic  energy = ',EHF + EckMP2
    write(*,'(A32,1X,F16.10)') ' MP2 total       energy = ',ENuc + EHF + EckMP2
    write(*,'(A32)')           '--------------------------'
    write(*,*)

  end if

end subroutine MP2
