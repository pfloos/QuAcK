subroutine complex_print_excitation_energies(method,manifold,nS,Om)

! Print excitation energies for a given spin manifold

  implicit none
  include 'parameters.h'

! Input variables

  character(len=*),intent(in)        :: method
  character(len=*),intent(in)        :: manifold
  integer,intent(in)                 :: nS
  complex*16,intent(in)              :: Om(nS)

! Local variables

  integer,parameter                  :: maxS = 300
  integer                            :: m

  write(*,*)
  write(*,*)'-------------------------------------------------------------'
  write(*,'(1X,A15,A15,A15,A9)') trim(method),' calculation: ',trim(manifold),' manifold'
  write(*,*)'-------------------------------------------------------------'
  write(*,'(1X,A1,1X,A5,1X,A1,1X,A23,1X,A1,1X,A23,1X,A1,1X)') &
            '|','State','|',' Re(Excitation energy)) (au) ','|',' Im(Excitation energy) (au) ','|'
  write(*,*)'-------------------------------------------------------------'

  do m=1,min(maxS,nS)
    write(*,'(1X,A1,1X,I5,1X,A1,1X,F23.6,1X,A1,1X,F23.6,1X,A1,1X)') & 
      '|',m,'|',real(Om(m)),'|',aimag(Om(m)),'|'
  end do

  write(*,*)'-------------------------------------------------------------'
  write(*,*)

end subroutine 
