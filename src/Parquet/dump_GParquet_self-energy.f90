subroutine dump_GParquet_self_energy(nOrb,nC,nO,nV,nR,Sig_2,Sig_eh,Sig_pp,Sig)

! Print the decomposition of the parquet self-energy

  implicit none

  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  double precision,intent(in)   :: Sig_2(nOrb)
  double precision,intent(in)   :: Sig_eh(nOrb,6)
  double precision,intent(in)   :: Sig_pp(nOrb,6)
  double precision,intent(in)   :: Sig(nOrb)

! Local variables

  integer                       :: i,a

! Output variables

!----------------------------------------!
! Dump parquet self-energy decomposition !
!----------------------------------------!

  write(*,*)'-----------------------------------------------------------'
  write(*,*)'| Parquet self-energy decomposition (in eV)               |'
  write(*,*)'-----------------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X)') &
            '|','#','|','Sig^(2)','|','Sig^eh','|','Sig^pp','|','Total','|'
  write(*,*)'-----------------------------------------------------------'

  ! Occupied states

  do i=nC+1,nO
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X)') &
    '|',i,'|',Sig_2(i)*HaToeV,'|',sum(Sig_eh(i,:))*HaToeV,'|',sum(Sig_pp(i,:))*HaToeV,'|',Sig(i)*HaToeV,'|'
  end do

  ! Fermi level

  write(*,*)'-----------------------------------------------------------'
  write(*,*)'-----------------------------------------------------------'

  ! Vacant states

  do a=nO+1,nOrb-nR
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X)') &
    '|',a,'|',Sig_2(a)*HaToeV,'|',sum(Sig_eh(a,:))*HaToeV,'|',sum(Sig_pp(a,:))*HaToeV,'|',Sig(a)*HaToeV,'|'
  end do
  write(*,*)'-----------------------------------------------------------'
  write(*,*)

!-------------------------------------------!
! Dump parquet eh self-energy decomposition !
!-------------------------------------------!

  write(*,*)'--------------------------------------------------------------------------------------------------'
  write(*,*)'| Parquet eh self-energy decomposition (in eV)                                                   |'
  write(*,*)'--------------------------------------------------------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X)') &
            '|','#','|','1 (3h2p)','|','2 (3h2p)','|','3 (3h2p)','|','4 (3p2h)','|','5 (3p2h)','|','6 (3p2h)','|','Total','|'
  write(*,*)'--------------------------------------------------------------------------------------------------'

  ! Occupied states

  do i=nC+1,nO
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X)') &
    '|',i,'|',Sig_eh(i,1)*HaToeV,'|',Sig_eh(i,2)*HaToeV,'|',Sig_eh(i,3)*HaToeV,'|',Sig_eh(i,4)*HaToeV,'|', & 
              Sig_eh(i,5)*HaToeV,'|',Sig_eh(i,6)*HaToeV,'|',sum(Sig_eh(i,:))*HaToeV,'|'
  end do

  ! Fermi level

  write(*,*)'--------------------------------------------------------------------------------------------------'
  write(*,*)'--------------------------------------------------------------------------------------------------'

  ! Vacant states

  do a=nO+1,nOrb-nR
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X)') &
    '|',a,'|',Sig_eh(a,1)*HaToeV,'|',Sig_eh(a,2)*HaToeV,'|',Sig_eh(a,3)*HaToeV,'|',Sig_eh(a,4)*HaToeV,'|', & 
              Sig_eh(a,5)*HaToeV,'|',Sig_eh(a,6)*HaToeV,'|',sum(Sig_eh(a,:))*HaToeV,'|'
  end do
  write(*,*)'--------------------------------------------------------------------------------------------------'
  write(*,*)  

!-------------------------------------------!
! Dump parquet pp self-energy decomposition !
!-------------------------------------------!

  write(*,*)'--------------------------------------------------------------------------------------------------'
  write(*,*)'| Parquet pp self-energy decomposition (in eV)                                                   |'
  write(*,*)'--------------------------------------------------------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X)') &
            '|','#','|','1 (4h1p)','|','2 (3h2p)','|','3 (3h2p)','|','4 (4p1h)','|','5 (3p2h)','|','6 (3p2h)','|','Total','|'
  write(*,*)'--------------------------------------------------------------------------------------------------'

  ! Occupied states

  do i=nC+1,nO
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X)') &
    '|',i,'|',Sig_pp(i,1)*HaToeV,'|',Sig_pp(i,2)*HaToeV,'|',Sig_pp(i,3)*HaToeV,'|',Sig_pp(i,4)*HaToeV,'|', & 
              Sig_pp(i,5)*HaToeV,'|',Sig_pp(i,6)*HaToeV,'|',sum(Sig_pp(i,:))*HaToeV,'|'
  end do

  ! Fermi level

  write(*,*)'--------------------------------------------------------------------------------------------------'
  write(*,*)'--------------------------------------------------------------------------------------------------'

  ! Vacant states

  do a=nO+1,nOrb-nR
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X)') &
    '|',a,'|',Sig_pp(a,1)*HaToeV,'|',Sig_pp(a,2)*HaToeV,'|',Sig_pp(a,3)*HaToeV,'|',Sig_pp(a,4)*HaToeV,'|', & 
              Sig_pp(a,5)*HaToeV,'|',Sig_pp(a,6)*HaToeV,'|',sum(Sig_pp(a,:))*HaToeV,'|'
  end do
  write(*,*)'--------------------------------------------------------------------------------------------------'
  write(*,*)  

end subroutine
