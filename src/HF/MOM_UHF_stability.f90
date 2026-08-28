subroutine MOM_UHF_stability(nBas,nC,nO,nV,nR,nS,nCVS,FC,eHF,ERI_aaaa,ERI_aabb,ERI_bbbb,occupations)

! Perform a stability analysis of the UHF solution obtained by MOM with eventually CVS and/or frozen core

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC(nspin)
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nV(nspin)
  integer,intent(in)            :: nR(nspin)
  integer,intent(in)            :: nS(nspin)
  integer,intent(in)            :: nCVS(nspin)
  integer,intent(in)            :: FC(nspin)
  integer,intent(in)            :: occupations(maxval(nO),nspin)

  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: ERI_aaaa(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_aabb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_bbbb(nBas,nBas,nBas,nBas)

! Local variables

  integer,parameter             :: maxS = 20
  integer                       :: ia,i
  integer                       :: ispin
  integer                       :: nSa,nSb,nSt,nFC(nspin)
  logical                       :: found

  double precision,allocatable  :: A(:,:)
  double precision,allocatable  :: B(:,:)
  double precision,allocatable  :: AB(:,:)
  double precision,allocatable  :: Om(:)
  
  integer,allocatable          :: virtuals(:,:)
  integer,allocatable          :: occupations_fc(:,:)

! Memory allocation

  allocate(virtuals(nBas - minval(nO),nspin))
  virtuals = 0
  do ispin=1,nspin
    call non_occupied(nO(ispin),nBas,occupations(1:nO(ispin),ispin),virtuals(1:nBas - nO(ispin),ispin))
  end do


  ! CVS

  print *, "No exications to the first", nCVS(1), "alpha orbital(s) are considered."
  print *, "No exications to the first", nCVS(2), "beta orbital(s) are considered."
  if(any(nC/=0)) then
    print *, "Do not use PyDuck frozen core with CVS !"
    stop
  end if

! Frozen Core

  nFC(1) = MERGE(1,0,FC(1)/=0) 
  nFC(2) = MERGE(1,0,FC(2)/=0)
  allocate(occupations_fc(maxval(nO-nFC),nspin))
  ! remove FC from occupations
  do ispin=1,nspin
    occupations_fc(1:nO(ispin)-nFC(ispin),ispin) = occupations(1:nO(ispin) - nFC(ispin), ispin) 
    found = .false.
    do i=1,nO(ispin)-1
      if(.not. found) then
        if(occupations(i,ispin)==FC(ispin)) then
          found = .true.
          occupations_fc(i,ispin) = occupations(i+1,ispin) 
        else
          occupations_fc(i,ispin) = occupations(i,ispin)
        endif
      else
        occupations_fc(i,ispin) = occupations(i+1,ispin) 
      endif 
    enddo
  enddo
  do ispin=1,nspin
    print *, "Not Frozen orbitals:"
    print *,occupations_fc(1:nO(ispin)-nFC(ispin),ispin)
  end do

  
!-------------------------------------------------------------!
! Stability analysis: Real UHF -> Real UHF  
!-------------------------------------------------------------!
 
  ispin = 1

  ! Memory allocation
  nSa = (nBas - nO(1) - nCVS(1))*(nO(1) - nFC(1))
  nSb = (nBas - nO(2) - nCVS(2))*(nO(2) - nFC(2))
  nSt = nSa + nSb
  allocate(A(nSt,nSt),B(nSt,nSt),AB(nSt,nSt),Om(nSt))

  call MOM_phULR_A(ispin,.false.,nBas,nC,nO,nV,nR,nSa,nSb,nSt,nCVS,nFC,occupations_fc,virtuals,1d0,eHF,ERI_aaaa,ERI_aabb,ERI_bbbb,A)
  call MOM_phULR_B(ispin,.false.,nBas,nC,nO,nV,nR,nSa,nSb,nSt,nCVS,nFC,occupations_fc,virtuals,1d0,ERI_aaaa,ERI_aabb,ERI_bbbb,B)

  AB(:,:) = A(:,:) + B(:,:)

  call diagonalize_matrix(nSt,AB,Om)

  Om(:) = 2d0*Om(:)

  write(*,*)'-------------------------------------------------------------'
  write(*,*)'|       Stability analysis: Real UHF -> Real UHF            |'
  write(*,*)'-------------------------------------------------------------'
  write(*,'(1X,A1,1X,A5,1X,A1,1X,A23,1X,A1,1X,A23,1X,A1,1X)') &
            '|','State','|',' Excitation energy (au) ','|',' Excitation energy (eV) ','|'
  write(*,*)'-------------------------------------------------------------'
  do ia=1,min(nSt,maxS)
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
! Stability analysis: Real UHF -> Complex UHF  
!-------------------------------------------------------------!
 
  AB(:,:) = A(:,:) - B(:,:)

  call diagonalize_matrix(nSt,AB,Om)
  Om(:) = 2d0*Om(:)

  write(*,*)'-------------------------------------------------------------'
  write(*,*)'|       Stability analysis: Real UHF -> Complex UHF         |'
  write(*,*)'-------------------------------------------------------------'
  write(*,'(1X,A1,1X,A5,1X,A1,1X,A23,1X,A1,1X,A23,1X,A1,1X)') &
            '|','State','|',' Excitation energy (au) ','|',' Excitation energy (eV) ','|'
  write(*,*)'-------------------------------------------------------------'
  do ia=1,min(nSt,maxS)
    write(*,'(1X,A1,1X,I5,1X,A1,1X,F23.6,1X,A1,1X,F23.6,1X,A1,1X)') &
      '|',ia,'|',Om(ia),'|',Om(ia)*HaToeV,'|'
  end do
  write(*,*)'-------------------------------------------------------------'

  if(minval(Om(:)) < 0d0) then

    write(*,'(1X,A40,1X)')        'Too bad, UHF solution is unstable!'
    write(*,'(1X,A40,1X,F15.10,A3)') 'Largest negative eigenvalue: ',Om(1),' au'

  else

    write(*,'(1X,A40,1X)')        'Well done, UHF solution is stable!'
    write(*,'(1X,A40,1X,F15.10,A3)') 'Smallest eigenvalue: ',Om(1),' au'

  end if
  write(*,*)'-------------------------------------------------------------'
  write(*,*)

  deallocate(A,B,AB,Om)


!-------------------------------------------------------------!
! Stability analysis: Real UHF -> Real GHF  
!-------------------------------------------------------------!
 
  ispin = 2

  ! Memory allocation
  nSa = (nBas - nO(1) - nCVS(1))*(nO(1) - nFC(1))
  nSb = (nBas - nO(2) - nCVS(2))*(nO(2) - nFC(2))
  nSt = nSa + nSb

  allocate(A(nSt,nSt),B(nSt,nSt),AB(nSt,nSt),Om(nSt))


  call MOM_phULR_A(ispin,.false.,nBas,nC,nO,nV,nR,nSa,nSb,nSt,nCVS,nFC,occupations_fc,virtuals,1d0,eHF,ERI_aaaa,ERI_aabb,ERI_bbbb,A)
  call MOM_phULR_B(ispin,.false.,nBas,nC,nO,nV,nR,nSa,nSb,nSt,nCVS,nFC,occupations_fc,virtuals,1d0,ERI_aaaa,ERI_aabb,ERI_bbbb,B)

  AB(:,:) = A(:,:) + B(:,:)

  call diagonalize_matrix(nSt,AB,Om)
  Om(:) = 2d0*Om(:)

  write(*,*)'-------------------------------------------------------------'
  write(*,*)'|       Stability analysis: Real UHF -> Real GHF            |'
  write(*,*)'-------------------------------------------------------------'
  write(*,'(1X,A1,1X,A5,1X,A1,1X,A23,1X,A1,1X,A23,1X,A1,1X)') &
            '|','State','|',' Excitation energy (au) ','|',' Excitation energy (eV) ','|'
  write(*,*)'-------------------------------------------------------------'
  do ia=1,min(nSt,maxS)
    write(*,'(1X,A1,1X,I5,1X,A1,1X,F23.6,1X,A1,1X,F23.6,1X,A1,1X)') &
      '|',ia,'|',Om(ia),'|',Om(ia)*HaToeV,'|'
  end do
  write(*,*)'-------------------------------------------------------------'

  if(minval(Om(:)) < 0d0) then

    write(*,'(1X,A40,1X)')        'Too bad, UHF solution is unstable!'
    write(*,'(1X,A40,1X,F15.10,A3)') 'Largest negative eigenvalue: ',Om(1),' au'

  else

    write(*,'(1X,A40,1X)')        'Well done, UHF solution is stable!'
    write(*,'(1X,A40,1X,F15.10,A3)') 'Smallest eigenvalue: ',Om(1),' au'

  end if
  write(*,*)'-------------------------------------------------------------'
  write(*,*)

  deallocate(A, B, AB, Om)
    
end subroutine 
