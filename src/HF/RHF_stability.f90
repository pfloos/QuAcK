subroutine RHF_stability(nBas,nC,nO,nV,nR,nS,eHF,ERI)

! Perform a stability analysis of the RHF solution

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS

  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer,parameter             :: maxS = 20
  integer                       :: ia
  integer                       :: ispin

  double precision,allocatable  :: A(:,:)
  double precision,allocatable  :: B(:,:)
  double precision,allocatable  :: AB(:,:)
  double precision,allocatable  :: Om(:)

! Memory allocation

  allocate(A(nS,nS),B(nS,nS),AB(nS,nS),Om(nS))
 
!-------------------------------------------------------------!
! Stability analysis: Real RHF -> Real RHF  
!-------------------------------------------------------------!
 
  ispin = 1

  call phLR_A(ispin,.false.,nBas,nC,nO,nV,nR,nS,1d0,eHF,ERI,A)
  call phLR_B(ispin,.false.,nBas,nC,nO,nV,nR,nS,1d0,ERI,B)

  AB(:,:) = A(:,:) + B(:,:)

  call diagonalize_matrix(nS,AB,Om)
  Om(:) = 0.5d0*Om(:)

  write(*,*)'-------------------------------------------------------------'
  write(*,*)'|       Stability analysis: Real RHF -> Real RHF            |'
  write(*,*)'-------------------------------------------------------------'
  write(*,'(1X,A1,1X,A5,1X,A1,1X,A23,1X,A1,1X,A23,1X,A1,1X)') &
            '|','State','|',' Excitation energy (au) ','|',' Excitation energy (eV) ','|'
  write(*,*)'-------------------------------------------------------------'
  do ia=1,min(nS,maxS)
    write(*,'(1X,A1,1X,I5,1X,A1,1X,F23.6,1X,A1,1X,F23.6,1X,A1,1X)') &
      '|',ia,'|',Om(ia),'|',Om(ia)*HaToeV,'|'
  enddo
  write(*,*)'-------------------------------------------------------------'

  if(minval(Om(:)) < 0d0) then

    write(*,'(1X,A40,1X)')        'Too bad, RHF solution is unstable!'
    write(*,'(1X,A40,1X,F15.10,A3)') 'Largest negative eigenvalue: ',Om(1),' au'

  else

    write(*,'(1X,A40,1X)')        'Well done, RHF solution is stable!'
    write(*,'(1X,A40,1X,F15.10,A3)') 'Smallest eigenvalue: ',Om(1),' au'

  end if
  write(*,*)'-------------------------------------------------------------'
  write(*,*)

!-------------------------------------------------------------!
! Stability analysis: Real RHF -> Complex RHF  
!-------------------------------------------------------------!
 
  AB(:,:) = A(:,:) - B(:,:)

  call diagonalize_matrix(nS,AB,Om)
  Om(:) = 0.5d0*Om(:)

  write(*,*)'-------------------------------------------------------------'
  write(*,*)'|       Stability analysis: Real RHF -> Complex RHF         |'
  write(*,*)'-------------------------------------------------------------'
  write(*,'(1X,A1,1X,A5,1X,A1,1X,A23,1X,A1,1X,A23,1X,A1,1X)') &
            '|','State','|',' Excitation energy (au) ','|',' Excitation energy (eV) ','|'
  write(*,*)'-------------------------------------------------------------'
  do ia=1,min(nS,maxS)
    write(*,'(1X,A1,1X,I5,1X,A1,1X,F23.6,1X,A1,1X,F23.6,1X,A1,1X)') &
      '|',ia,'|',Om(ia),'|',Om(ia)*HaToeV,'|'
  enddo
  write(*,*)'-------------------------------------------------------------'

  if(minval(Om(:)) < 0d0) then

    write(*,'(1X,A40,1X)')        'Too bad, RHF solution is unstable!'
    write(*,'(1X,A40,1X,F15.10,A3)') 'Largest negative eigenvalue: ',Om(1),' au'

  else

    write(*,'(1X,A40,1X)')        'Well done, RHF solution is stable!'
    write(*,'(1X,A40,1X,F15.10,A3)') 'Smallest eigenvalue: ',Om(1),' au'

  end if
  write(*,*)'-------------------------------------------------------------'
  write(*,*)

!-------------------------------------------------------------!
! Stability analysis: Real RHF -> Real UHF  
!-------------------------------------------------------------!
 
  ispin = 2

  call phLR_A(ispin,.false.,nBas,nC,nO,nV,nR,nS,1d0,eHF,ERI,A)
  call phLR_B(ispin,.false.,nBas,nC,nO,nV,nR,nS,1d0,ERI,B)

  AB(:,:) = A(:,:) + B(:,:)

  call diagonalize_matrix(nS,AB,Om)
  Om(:) = 0.5d0*Om(:)

  write(*,*)'-------------------------------------------------------------'
  write(*,*)'|       Stability analysis: Real RHF -> Real UHF            |'
  write(*,*)'-------------------------------------------------------------'
  write(*,'(1X,A1,1X,A5,1X,A1,1X,A23,1X,A1,1X,A23,1X,A1,1X)') &
            '|','State','|',' Excitation energy (au) ','|',' Excitation energy (eV) ','|'
  write(*,*)'-------------------------------------------------------------'
  do ia=1,min(nS,maxS)
    write(*,'(1X,A1,1X,I5,1X,A1,1X,F23.6,1X,A1,1X,F23.6,1X,A1,1X)') &
      '|',ia,'|',Om(ia),'|',Om(ia)*HaToeV,'|'
  enddo
  write(*,*)'-------------------------------------------------------------'

  if(minval(Om(:)) < 0d0) then

    write(*,'(1X,A40,1X)')        'Too bad, RHF solution is unstable!'
    write(*,'(1X,A40,1X,F15.10,A3)') 'Largest negative eigenvalue: ',Om(1),' au'

  else

    write(*,'(1X,A40,1X)')        'Well done, RHF solution is stable!'
    write(*,'(1X,A40,1X,F15.10,A3)') 'Smallest eigenvalue: ',Om(1),' au'

  end if
  write(*,*)'-------------------------------------------------------------'
  write(*,*)
    
end subroutine 
