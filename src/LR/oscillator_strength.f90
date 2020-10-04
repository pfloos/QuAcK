subroutine oscillator_strength(nBas,nC,nO,nV,nR,nS,dipole_int,Omega,XpY,XmY,os)

! Compute linear response

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision              :: dipole_int(nBas,nBas,ncart)
  double precision,intent(in)   :: Omega(nS)
  double precision,intent(in)   :: XpY(nS,nS)
  double precision,intent(in)   :: XmY(nS,nS)
  
! Local variables

  logical                       :: debug = .false.
  integer                       :: ia,jb,i,j,a,b
  integer                       :: ixyz

  double precision,allocatable  :: f(:,:)

! Output variables

  double precision              :: os(nS)

! Memory allocation

  allocate(f(nS,ncart))

! Initialization
   
  f(:,:) = 0d0

! Compute dipole moments and oscillator strengths

  do ia=1,nS
    do ixyz=1,ncart
      jb = 0
      do j=nC+1,nO
        do b=nO+1,nBas-nR
          jb = jb + 1
          f(ia,ixyz) = f(ia,ixyz) + dipole_int(j,b,ixyz)*XpY(ia,jb)
        end do
      end do
    end do
  end do
  f(:,:) = sqrt(2d0)*f(:,:)

  do ia=1,nS
    os(ia) = 2d0/3d0*Omega(ia)*sum(f(ia,:)**2)
  end do
    
  if(debug) then

    write(*,*) '------------------------'
    write(*,*) ' Dipole moments (X Y Z) '
    write(*,*) '------------------------'
    call matout(nS,ncart,f)
    write(*,*)

    write(*,*) '----------------------'
    write(*,*) ' Oscillator strengths '
    write(*,*) '----------------------'
    call matout(nS,1,os)
    write(*,*)

  end if

end subroutine oscillator_strength
