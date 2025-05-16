subroutine dump_GParquet_self_energy(nOrb,nC,nO,nV,nR,Sig2,Sigeh,Sigpp,Sig)

! Print the decomposition of the parquet self-energy

  implicit none

  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  double precision,intent(in)   :: Sig2(nOrb)
  double precision,intent(in)   :: Sigeh(nOrb)
  double precision,intent(in)   :: Sigpp(nOrb)
  double precision,intent(in)   :: Sig(nOrb)

! Local variables

  integer                       :: i,a

! Output variables

! Initalization

! Dump parquet self-energy decomposition

  write(*,*)'-----------------------------------------------------------'
  write(*,*)'| Parquet self-energy decomposition (in eV)               |'
  write(*,*)'-----------------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X)') &
            '|','#','|','Sig^(2)','|','Sig_eh','|','Sig_pp','|','Total','|'
  write(*,*)'-----------------------------------------------------------'

  ! Occupied states

  do i=nC+1,nO
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X)') &
    '|',i,'|',Sig2(i)*HaToeV,'|',Sigeh(i)*HaToeV,'|',Sigpp(i)*HaToeV,'|',Sig(i)*HaToeV,'|'
  end do

  ! Fermi level

  write(*,*)'-----------------------------------------------------------'
  write(*,*)'-----------------------------------------------------------'

  ! Vacant states

  do a=nO+1,nOrb-nR
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X)') &
    '|',a,'|',Sig2(a)*HaToeV,'|',Sigeh(a)*HaToeV,'|',Sigpp(a)*HaToeV,'|',Sig(a)*HaToeV,'|'
  end do
  write(*,*)'-----------------------------------------------------------'
  write(*,*)

end subroutine
