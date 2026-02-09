subroutine dump_RParquet_self_energy(nOrb,nC,nO,nV,nR,Sig2d,Sig2x,Sig1eh,Sig3eh,Sig1pp,Sig3pp,Sig)

! Print the decomposition of the parquet self-energy

  implicit none

  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  double precision,intent(in)   :: Sig2d(nOrb)
  double precision,intent(in)   :: Sig2x(nOrb)
  double precision,intent(in)   :: Sig1eh(nOrb)
  double precision,intent(in)   :: Sig3eh(nOrb)
  double precision,intent(in)   :: Sig1pp(nOrb)
  double precision,intent(in)   :: Sig3pp(nOrb)
  double precision,intent(in)   :: Sig(nOrb)

! Local variables

  integer                       :: i,a

! Output variables

! Initalization

! Dump parquet self-energy decomposition

  write(*,*)'--------------------------------------------------------------------------------------------------'
  write(*,*)'| Parquet self-energy decomposition (in eV)                                                      |'
  write(*,*)'--------------------------------------------------------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X)') &
            '|','#','|','Sig^(2d)','|','Sig^(2x)','|','^1 Sig_eh','|','^3 Sig_eh','|','^1 Sig_pp','|','^3 Sig_pp','|','Total','|'
  write(*,*)'--------------------------------------------------------------------------------------------------'

  ! Occupied states

  do i=nC+1,nO
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X,&
                 &F10.6,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X)')  &
    '|',i,'|',Sig2d(i)*HaToeV,'|',Sig2x(i)*HaToeV,'|',Sig1eh(i)*HaToeV,'|',      &
    Sig3eh(i)*HaToeV,'|',Sig1pp(i)*HaToeV,'|',Sig3pp(i)*HaToeV,'|',Sig(i)*HaToeV,'|'
  end do

  ! Fermi level

  write(*,*)'--------------------------------------------------------------------------------------------------'
  write(*,*)'--------------------------------------------------------------------------------------------------'

  ! Vacant states

  do a=nO+1,nOrb-nR
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,&
            &A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X)')                       &
    '|',a,'|',Sig2d(a)*HaToeV,'|',Sig2x(a)*HaToeV,'|',Sig1eh(a)*HaToeV,'|',              &
    Sig3eh(a)*HaToeV,'|',Sig1pp(a)*HaToeV,'|',Sig3pp(a)*HaToeV,'|',Sig(a)*HaToeV,'|'
  end do
  write(*,*)'--------------------------------------------------------------------------------------------------'
  write(*,*)

end subroutine
