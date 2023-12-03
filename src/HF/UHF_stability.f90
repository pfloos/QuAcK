subroutine UHF_stability(nBas,nC,nO,nV,nR,nS,eHF,ERI_aaaa,ERI_aabb,ERI_bbbb)

! Perform a stability analysis of the UHF solution

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC(nspin)
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nV(nspin)
  integer,intent(in)            :: nR(nspin)
  integer,intent(in)            :: nS(nspin)

  double precision,intent(in)   :: eHF(nBas,nspin)
  double precision,intent(in)   :: ERI_aaaa(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_aabb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_bbbb(nBas,nBas,nBas,nBas)

! Local variables

  integer,parameter             :: maxS = 20
  integer                       :: ia
  integer                       :: ispin

  integer                       :: nS_aa,nS_bb,nS_sc
  double precision,allocatable  :: Om_sc(:)
  double precision,allocatable  :: A_sc(:,:)
  double precision,allocatable  :: B_sc(:,:)
  double precision,allocatable  :: AB_sc(:,:)

  integer                       :: nS_ab,nS_ba,nS_sf
  double precision,allocatable  :: Om_sf(:)
  double precision,allocatable  :: A_sf(:,:)
  double precision,allocatable  :: B_sf(:,:)
  double precision,allocatable  :: AB_sf(:,:)

! Menory allocation

  nS_aa = nS(1)
  nS_bb = nS(2)
  nS_sc = nS_aa + nS_bb

  allocate(Om_sc(nS_sc),A_sc(nS_sc,nS_sc),B_sc(nS_sc,nS_sc),AB_sc(nS_sc,nS_sc))
 
!-------------------------------------------------------------!
! Stability analysis: Real UHF -> Real UHF  
!-------------------------------------------------------------!
 
  ispin = 1

  call phULR_A(ispin,.false.,nBas,nC,nO,nV,nR,nS_aa,nS_bb,nS_sc,1d0,eHF,ERI_aaaa,ERI_aabb,ERI_bbbb,A_sc)
  call phULR_B(ispin,.false.,nBas,nC,nO,nV,nR,nS_aa,nS_bb,nS_sc,1d0,ERI_aaaa,ERI_aabb,ERI_bbbb,B_sc)

  AB_sc(:,:) = A_sc(:,:) + B_sc(:,:)

  call diagonalize_matrix(nS_sc,AB_sc,Om_sc)
  Om_sc(:) = 2d0*Om_sc(:)

  write(*,*)'-------------------------------------------------------------'
  write(*,*)'|       Stability analysis: Real UHF -> Real UHF            |'
  write(*,*)'-------------------------------------------------------------'
  write(*,'(1X,A1,1X,A5,1X,A1,1X,A23,1X,A1,1X,A23,1X,A1,1X)') &
            '|','State','|',' Excitation energy (au) ','|',' Excitation energy (eV) ','|'
  write(*,*)'-------------------------------------------------------------'
  do ia=1,min(nS_sc,maxS)
    write(*,'(1X,A1,1X,I5,1X,A1,1X,F23.6,1X,A1,1X,F23.6,1X,A1,1X)') &
      '|',ia,'|',Om_sc(ia),'|',Om_sc(ia)*HaToeV,'|'
  end do
  write(*,*)'-------------------------------------------------------------'

  if(minval(Om_sc(:)) < 0d0) then

    write(*,'(1X,A40,1X)')        'Too bad, UHF solution is unstable!'
    write(*,'(1X,A40,1X,F15.10,A3)') 'Largest negative eigenvalue: ',Om_sc(1),' au'

  else

    write(*,'(1X,A40,1X)')        'Well done, UHF solution is stable!'
    write(*,'(1X,A40,1X,F15.10,A3)') 'Smallest eigenvalue: ',Om_sc(1),' au'

  end if
  write(*,*)'-------------------------------------------------------------'
  write(*,*)

!-------------------------------------------------------------!
! Stability analysis: Real UHF -> Complex UHF  
!-------------------------------------------------------------!
 
  AB_sc(:,:) = A_sc(:,:) - B_sc(:,:)

  call diagonalize_matrix(nS_sc,AB_sc,Om_sc)
  Om_sc(:) = 2d0*Om_sc(:)

  write(*,*)'-------------------------------------------------------------'
  write(*,*)'|       Stability analysis: Real UHF -> Complex UHF         |'
  write(*,*)'-------------------------------------------------------------'
  write(*,'(1X,A1,1X,A5,1X,A1,1X,A23,1X,A1,1X,A23,1X,A1,1X)') &
            '|','State','|',' Excitation energy (au) ','|',' Excitation energy (eV) ','|'
  write(*,*)'-------------------------------------------------------------'
  do ia=1,min(nS_sc,maxS)
    write(*,'(1X,A1,1X,I5,1X,A1,1X,F23.6,1X,A1,1X,F23.6,1X,A1,1X)') &
      '|',ia,'|',Om_sc(ia),'|',Om_sc(ia)*HaToeV,'|'
  end do
  write(*,*)'-------------------------------------------------------------'

  if(minval(Om_sc(:)) < 0d0) then

    write(*,'(1X,A40,1X)')        'Too bad, UHF solution is unstable!'
    write(*,'(1X,A40,1X,F15.10,A3)') 'Largest negative eigenvalue: ',Om_sc(1),' au'

  else

    write(*,'(1X,A40,1X)')        'Well done, UHF solution is stable!'
    write(*,'(1X,A40,1X,F15.10,A3)') 'Smallest eigenvalue: ',Om_sc(1),' au'

  end if
  write(*,*)'-------------------------------------------------------------'
  write(*,*)

! Menory (de)allocation


  nS_ab = (nO(1) - nC(1))*(nV(2) - nR(2))
  nS_ba = (nO(2) - nC(2))*(nV(1) - nR(1))
  nS_sf = nS_ab + nS_ba

  deallocate(Om_sc,A_sc,B_sc,AB_sc)
  allocate(Om_sf(nS_sf),A_sf(nS_sf,nS_sf),B_sf(nS_sf,nS_sf),AB_sf(nS_sf,nS_sf))

!-------------------------------------------------------------!
! Stability analysis: Real UHF -> Real GHF  
!-------------------------------------------------------------!
 
  ispin = 2

  call phULR_A(ispin,.false.,nBas,nC,nO,nV,nR,nS_ab,nS_ba,nS_sf,1d0,eHF,ERI_aaaa,ERI_aabb,ERI_bbbb,A_sf)
  call phULR_B(ispin,.false.,nBas,nC,nO,nV,nR,nS_ab,nS_ba,nS_sf,1d0,ERI_aaaa,ERI_aabb,ERI_bbbb,B_sf)

  AB_sf(:,:) = A_sf(:,:) + B_sf(:,:)

  call diagonalize_matrix(nS_sf,AB_sf,Om_sf)
  Om_sf(:) = 2d0*Om_sf(:)

  write(*,*)'-------------------------------------------------------------'
  write(*,*)'|       Stability analysis: Real UHF -> Real GHF            |'
  write(*,*)'-------------------------------------------------------------'
  write(*,'(1X,A1,1X,A5,1X,A1,1X,A23,1X,A1,1X,A23,1X,A1,1X)') &
            '|','State','|',' Excitation energy (au) ','|',' Excitation energy (eV) ','|'
  write(*,*)'-------------------------------------------------------------'
  do ia=1,min(nS_sf,maxS)
    write(*,'(1X,A1,1X,I5,1X,A1,1X,F23.6,1X,A1,1X,F23.6,1X,A1,1X)') &
      '|',ia,'|',Om_sf(ia),'|',Om_sf(ia)*HaToeV,'|'
  end do
  write(*,*)'-------------------------------------------------------------'

  if(minval(Om_sf(:)) < 0d0) then

    write(*,'(1X,A40,1X)')        'Too bad, UHF solution is unstable!'
    write(*,'(1X,A40,1X,F15.10,A3)') 'Largest negative eigenvalue: ',Om_sf(1),' au'

  else

    write(*,'(1X,A40,1X)')        'Well done, UHF solution is stable!'
    write(*,'(1X,A40,1X,F15.10,A3)') 'Smallest eigenvalue: ',Om_sf(1),' au'

  end if
  write(*,*)'-------------------------------------------------------------'
  write(*,*)
    
end subroutine 
