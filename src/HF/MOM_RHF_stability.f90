subroutine MOM_RHF_stability(nBas,nC,nO,nV,nR,nS,nCVS,FC,eHF,ERI,occupations)

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
  integer,intent(in)            :: nCVS
  integer,intent(in)            :: FC
  integer,intent(in)            :: occupations(nO)

  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer,parameter             :: maxS = 20
  integer                       :: ia,i
  integer                       :: ispin
  integer                       :: nSCVS,nFC
  logical                       :: found

  double precision,allocatable  :: A(:,:)
  double precision,allocatable  :: B(:,:)
  double precision,allocatable  :: AB(:,:)
  double precision,allocatable  :: Om(:)
  
  integer,allocatable          :: virtuals(:)
  integer,allocatable          :: occupations_fc(:)

! Memory allocation

  allocate(virtuals(nV))
  call non_occupied(nO,nBas,occupations,virtuals)
 
! CVS
  
  print *, "No exications to the first", nCVS, "orbital(s) are considered."
  if(nC/=0) then
    print *, "Do not use PyDuck frozen core with CVS !"
    stop
  end if

! Frozen Core

  nFC = MERGE(1,0,FC/=0) 
  allocate(occupations_fc(nO-nFC))
  ! remove FC from occupations
  do ispin=1,nspin
    occupations_fc(1:nO-nFC) = occupations(1:nO - nFC) 
    found = .false.
    do i=1,nO-1
      if(.not. found) then
        if(occupations(i)==FC) then
          found = .true.
          occupations_fc(i) = occupations(i+1) 
        else
          occupations_fc(i) = occupations(i)
        endif
      else
        occupations_fc(i) = occupations(i+1) 
      endif 
    enddo
  enddo
  print *, "Not Frozen orbitals:"
  print *,occupations_fc(1:nO-nFC)
  

  nSCVS = (nV - nCVS)*(nO - nFC)
  allocate(A(nSCVS,nSCVS), B(nSCVS,nSCVS), AB(nSCVS,nSCVS), Om(nSCVS))
  
!-------------------------------------------------------------!
! Stability analysis: Real RHF -> Real RHF  
!-------------------------------------------------------------!
 
  ispin = 1
  
  print *, nSCVS,nFC,nCVS
  call CVS_phRLR_A(ispin,.false.,nBas,nC,nO,nV,nR,nSCVS,nCVS,nFC,occupations_fc,virtuals,1d0,eHF,ERI,A)
  call CVS_phRLR_B(ispin,.false.,nBas,nC,nO,nV,nR,nSCVS,nCVS,nFC,occupations_fc,virtuals,1d0,ERI,B)

  AB(:,:) = A(:,:) + B(:,:)

  call diagonalize_matrix(nSCVS,AB,Om)
  Om(:) = 2d0*Om(:)

  write(*,*)'-------------------------------------------------------------'
  write(*,*)'|       Stability analysis: Real RHF -> Real RHF            |'
  write(*,*)'-------------------------------------------------------------'
  write(*,'(1X,A1,1X,A5,1X,A1,1X,A23,1X,A1,1X,A23,1X,A1,1X)') &
            '|','State','|',' Excitation energy (au) ','|',' Excitation energy (eV) ','|'
  write(*,*)'-------------------------------------------------------------'
  do ia=1,min(nSCVS,maxS)
    write(*,'(1X,A1,1X,I5,1X,A1,1X,F23.6,1X,A1,1X,F23.6,1X,A1,1X)') &
      '|',ia,'|',Om(ia),'|',Om(ia)*HaToeV,'|'
  end do
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

  call diagonalize_matrix(nSCVS,AB,Om)
  Om(:) = 2d0*Om(:)

  write(*,*)'-------------------------------------------------------------'
  write(*,*)'|       Stability analysis: Real RHF -> Complex RHF         |'
  write(*,*)'-------------------------------------------------------------'
  write(*,'(1X,A1,1X,A5,1X,A1,1X,A23,1X,A1,1X,A23,1X,A1,1X)') &
            '|','State','|',' Excitation energy (au) ','|',' Excitation energy (eV) ','|'
  write(*,*)'-------------------------------------------------------------'
  do ia=1,min(nS,maxS)
    write(*,'(1X,A1,1X,I5,1X,A1,1X,F23.6,1X,A1,1X,F23.6,1X,A1,1X)') &
      '|',ia,'|',Om(ia),'|',Om(ia)*HaToeV,'|'
  end do
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

  call CVS_phRLR_A(ispin,.false.,nBas,nC,nO,nV,nR,nSCVS,nCVS,nFC,occupations_fc,virtuals,1d0,eHF,ERI,A)
  call CVS_phRLR_B(ispin,.false.,nBas,nC,nO,nV,nR,nSCVS,nCVS,nFC,occupations_fc,virtuals,1d0,ERI,B)

  AB(:,:) = A(:,:) + B(:,:)

  call diagonalize_matrix(nSCVS,AB,Om)
  Om(:) = 2d0*Om(:)

  write(*,*)'-------------------------------------------------------------'
  write(*,*)'|       Stability analysis: Real RHF -> Real UHF            |'
  write(*,*)'-------------------------------------------------------------'
  write(*,'(1X,A1,1X,A5,1X,A1,1X,A23,1X,A1,1X,A23,1X,A1,1X)') &
            '|','State','|',' Excitation energy (au) ','|',' Excitation energy (eV) ','|'
  write(*,*)'-------------------------------------------------------------'
  do ia=1,min(nSCVS,maxS)
    write(*,'(1X,A1,1X,I5,1X,A1,1X,F23.6,1X,A1,1X,F23.6,1X,A1,1X)') &
      '|',ia,'|',Om(ia),'|',Om(ia)*HaToeV,'|'
  end do
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

  deallocate(A, B, AB, Om)
    
end subroutine 
