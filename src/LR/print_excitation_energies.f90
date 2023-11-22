subroutine print_excitation_energies(method,manifold,nS,Om)

! Print excitation energies for a given spin manifold

  implicit none
  include 'parameters.h'

! Input variables

  character(len=*),intent(in)        :: method
  character(len=*),intent(in)        :: manifold
  integer,intent(in)                 :: nS
  double precision,intent(in)        :: Om(nS)

! Local variables

  integer,parameter                  :: maxS = 20
  integer                            :: m

  write(*,*)
  write(*,*)'-------------------------------------------------------------'
  write(*,'(1X,A15,A15,A15,A9)') trim(method),' calculation: ',trim(manifold),' manifold'
  write(*,*)'-------------------------------------------------------------'
  write(*,'(1X,A1,1X,A5,1X,A1,1X,A23,1X,A1,1X,A23,1X,A1,1X)') &
            '|','State','|',' Excitation energy (au) ','|',' Excitation energy (eV) ','|'
  write(*,*)'-------------------------------------------------------------'

  do m=1,min(maxS,nS)
    write(*,'(1X,A1,1X,I5,1X,A1,1X,F23.6,1X,A1,1X,F23.6,1X,A1,1X)') & 
      '|',m,'|',Om(m),'|',Om(m)*HaToeV,'|'
  end do

  write(*,*)'-------------------------------------------------------------'
  write(*,*)

end subroutine 
