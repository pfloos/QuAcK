subroutine GMP2(dotest,regularize,nBas,nC,nO,nV,nR,ERI,ENuc,EGHF,eHF,EcMP2)

! Perform second-order Moller-Plesset calculation with and without regularizers

  implicit none

! Input variables

  logical,intent(in)            :: dotest

  logical,intent(in)            :: regularize
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: EGHF
  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: i,j,a,b

  double precision              :: kappa,sigm1,sigm2
  double precision              :: num
  double precision              :: Dijab
  double precision              :: fs,fs2,fk

  double precision              :: E2d,E2ds,E2ds2,E2dk
  double precision              :: E2x,E2xs,E2xs2,E2xk

  double precision              :: EcsMP2,Ecs2MP2,EckMP2

! Output variables

  double precision,intent(out)  :: EcMP2

! Hello world

  write(*,*)
  write(*,*)'*******************************'
  write(*,*)'* Generalized MP2 Calculation |'
  write(*,*)'*******************************'
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

          Dijab = eHF(a) + eHF(b) - eHF(i) - eHF(j) 

          ! Second-order ring diagram

          fs  = (1d0 - exp(-sigm1*Dijab))/Dijab
          fs2 = (1d0 - exp(-sigm2*Dijab*Dijab))/Dijab
          fk  = (1d0 - exp(-kappa*Dijab))**2/Dijab

          num   = 0.25d0*(ERI(i,j,a,b)**2 + ERI(i,j,b,a)**2)
          E2d   = E2d   - num/Dijab
          E2ds  = E2ds  - num*fs 
          E2ds2 = E2ds2 - num*fs2
          E2dk  = E2dk  - num*fk

          ! Second-order exchange diagram

          num   = 0.5d0*ERI(i,j,a,b)*ERI(i,j,b,a)
          E2x   = E2x   + num/Dijab
          E2xs  = E2xs  + num*fs 
          E2xs2 = E2xs2 + num*fs2
          E2xk  = E2xk  + num*fk

        end do
      end do
    end do
  end do

  EcMP2   = E2d   + E2x
  EcsMP2  = E2ds  + E2xs
  Ecs2MP2 = E2ds2 + E2xs2
  EckMP2  = E2dk  + E2xk

!------------!
! MP2 energy !
!------------!

  write(*,*)
  write(*,'(A32)')           '---------------------------'
  write(*,'(A32)')           ' GMP2 calculation          '
  write(*,'(A32)')           '---------------------------'
  write(*,'(A32,1X,F16.10)') ' GMP2 correlation energy = ',EcMP2
  write(*,'(A32,1X,F16.10)') ' Direct part             = ',E2d
  write(*,'(A32,1X,F16.10)') ' Exchange part           = ',E2x
  write(*,'(A32)')           '---------------------------'
  write(*,'(A32,1X,F16.10)') ' GMP2 electronic  energy = ',EGHF + EcMP2
  write(*,'(A32,1X,F16.10)') ' GMP2 total       energy = ',ENuc + EGHF + EcMP2
  write(*,'(A32)')           '---------------------------'
  write(*,*)

  if(regularize) then

!-------------------!
! sigma1-MP2 energy !
!-------------------!

    write(*,*)
    write(*,'(A32)')           '---------------------------'
    write(*,'(A32)')           ' sigma-GMP2 calculation    '
    write(*,'(A32)')           '---------------------------'
    write(*,'(A32,1X,F16.10)') ' GMP2 correlation energy = ',EcsMP2
    write(*,'(A32,1X,F16.10)') ' Direct part             = ',E2ds
    write(*,'(A32,1X,F16.10)') ' Exchange part           = ',E2xs
    write(*,'(A32)')           '---------------------------'
    write(*,'(A32,1X,F16.10)') ' GMP2 electronic  energy = ',EGHF + EcsMP2
    write(*,'(A32,1X,F16.10)') ' GMP2 total       energy = ',ENuc + EGHF + EcsMP2
    write(*,'(A32)')           '---------------------------'
    write(*,*)

!--------------------!
! sigma^2-MP2 energy !
!--------------------!

    write(*,*)
    write(*,'(A32)')           '---------------------------'
    write(*,'(A32)')           ' sigma^2-GMP2 calculation  '
    write(*,'(A32)')           '---------------------------'
    write(*,'(A32,1X,F16.10)') ' GMP2 correlation energy = ',Ecs2MP2
    write(*,'(A32,1X,F16.10)') ' Direct part             = ',E2ds2
    write(*,'(A32,1X,F16.10)') ' Exchange part           = ',E2xs2
    write(*,'(A32)')           '---------------------------'
    write(*,'(A32,1X,F16.10)') ' GMP2 electronic  energy = ',EGHF + Ecs2MP2
    write(*,'(A32,1X,F16.10)') ' GMP2 total       energy = ',ENuc + EGHF + Ecs2MP2
    write(*,'(A32)')           '---------------------------'
    write(*,*)

!------------------!
! kappa-MP2 energy !
!------------------!

    write(*,*)
    write(*,'(A32)')           '---------------------------'
    write(*,'(A32)')           ' kappa-GMP2 calculation    '
    write(*,'(A32)')           '---------------------------'
    write(*,'(A32,1X,F16.10)') ' GMP2 correlation energy = ',EckMP2
    write(*,'(A32,1X,F16.10)') ' Direct part             = ',E2dk
    write(*,'(A32,1X,F16.10)') ' Exchange part           = ',-E2xk
    write(*,'(A32)')           '---------------------------'
    write(*,'(A32,1X,F16.10)') ' GMP2 electronic  energy = ',EGHF + EckMP2
    write(*,'(A32,1X,F16.10)') ' GMP2 total       energy = ',ENuc + EGHF + EckMP2
    write(*,'(A32)')           '---------------------------'
    write(*,*)

  end if

  if(dotest) then

    call dump_test_value('G','MP2 correlation energy',EcMP2)

  end if


end subroutine 
