subroutine print_excitation(method,ispin,nS,Omega)

! Print excitation energies for a given spin manifold

  implicit none
  include 'parameters.h'

! Input variables

  character*6,intent(in)             :: method
  integer,intent(in)                 :: ispin,nS
  double precision,intent(in)        :: Omega(nS)

! Local variables

  character*7                        :: spin_manifold
  integer,parameter                  :: maxS = 32
  integer                            :: ia

  if(ispin == 1) spin_manifold = 'singlet'
  if(ispin == 2) spin_manifold = 'triplet'
  if(ispin == 3) spin_manifold = 'spinorb'

  write(*,*)
  write(*,*)'-------------------------------------------------------------'
  write(*,'(1X,A1,1X,A14,A14,A7,A9,A15)')'|',method,' calculation: ',spin_manifold,' manifold','         |'
  write(*,*)'-------------------------------------------------------------'
  write(*,'(1X,A1,1X,A5,1X,A1,1X,A23,1X,A1,1X,A23,1X,A1,1X)') &
            '|','State','|',' Excitation energy (au) ','|',' Excitation energy (eV) ','|'
  write(*,*)'-------------------------------------------------------------'

  do ia=1,min(nS,maxS)
    write(*,'(1X,A1,1X,I5,1X,A1,1X,F23.6,1X,A1,1X,F23.6,1X,A1,1X)') & 
      '|',ia,'|',Omega(ia),'|',Omega(ia)*HaToeV,'|'
  enddo

  write(*,*)'-------------------------------------------------------------'
  write(*,*)

end subroutine print_excitation


