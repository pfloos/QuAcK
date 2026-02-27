subroutine CVS_phLR_oscillator_strength(nOrb,nC,nO,nV,nR,nS,nCVS,nFC,maxS,occupations,virtuals,dipole_int,Om,XpY,XmY,os)

! Compute oscillator strength from a ph linear response calculation

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  integer,intent(in)            :: nCVS,nFC
  integer,intent(in)            :: occupations(nO-nFC),virtuals(nOrb-nO)
  integer,intent(in)            :: maxS
  double precision              :: dipole_int(nOrb,nOrb,ncart)
  double precision,intent(in)   :: Om(nS)
  double precision,intent(in)   :: XpY(nS,nS)
  double precision,intent(in)   :: XmY(nS,nS)
  
! Local variables

  integer                       :: m,jb,i,j,a,b
  integer                       :: ixyz

  double precision,allocatable  :: f(:,:)

! Output variables

  double precision,intent(out)  :: os(nS)

! Memory allocation

  allocate(f(maxS,ncart))

! Initialization
   
  f(:,:) = 0d0

! Compute dipole moments and oscillator strengths

  do m=1,maxS
    do ixyz=1,ncart
      jb = 0
      do j=1,nO-nFC
        do b=nCVS+1,nOrb-nO
          jb = jb + 1
          f(m,ixyz) = f(m,ixyz) + dipole_int(occupations(j),virtuals(b),ixyz)*XpY(m,jb)
        end do
      end do
    end do
  end do
  f(:,:) = sqrt(2d0)*f(:,:)

  do m=1,maxS
    os(m) = 2d0/3d0*Om(m)*sum(f(m,:)**2)
  end do
 
  write(*,*) '---------------------------------------------------------------'
  write(*,*) '                Transition dipole moment (au)                  '
  write(*,*) '---------------------------------------------------------------'
  write(*,'(A3,5A12)') '#','X','Y','Z','dip. str.','osc. str.'
  write(*,*) '---------------------------------------------------------------'
  do m=1,maxS
    write(*,'(I3,5F12.6)') m,(f(m,ixyz),ixyz=1,ncart),sum(f(m,:)**2),os(m)
  end do
  write(*,*) '---------------------------------------------------------------'
  write(*,*)

end subroutine 
