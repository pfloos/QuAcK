subroutine CVS_phULR_oscillator_strength(nBas,nC,nO,nV,nR,nS,nSa,nSb,nSt,nCVS,nFC,maxS,occupations,virtuals,dipole_int_aa,dipole_int_bb,Omega,XpY,XmY,os)

! Compute linear response

  implicit none
  include 'parameters.h'

! Input variables


  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC(nspin)
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nV(nspin)
  integer,intent(in)            :: nR(nspin)
  integer,intent(in)            :: nS(nspin)
  integer,intent(in)            :: nSa
  integer,intent(in)            :: nSb
  integer,intent(in)            :: nSt
  integer,intent(in)            :: maxS
  integer,intent(in)            :: nCVS(nspin),nFC(nspin)
  integer,intent(in)            :: occupations(maxval(nO),nspin)
  integer,intent(in)            :: virtuals(nBas - minval(nO),nspin)
  double precision              :: dipole_int_aa(nBas,nBas,ncart)
  double precision              :: dipole_int_bb(nBas,nBas,ncart)
  double precision,intent(in)   :: Omega(nSt)
  double precision,intent(in)   :: XpY(nSt,nSt)
  double precision,intent(in)   :: XmY(nSt,nSt)

! Local variables

  integer                       :: ia,jb,i,j,a,b
  integer                       :: ixyz

  double precision,allocatable  :: f(:,:)

! Output variables

  double precision              :: os(maxS)

! Memory allocation

  allocate(f(maxS,ncart))

! Initialization
   
  f(:,:) = 0d0

! Compute dipole moments and oscillator strengths

  do ia=1,maxS
    do ixyz=1,ncart

      jb = 0
      do j=1,nO(1) - nFC(1)
        do b=nCVS(1) +1,nBas-nO(1)
          jb = jb + 1
          f(ia,ixyz) = f(ia,ixyz) + dipole_int_aa(occupations(j,1),virtuals(b,1),ixyz)*XpY(ia,jb)
        end do
      end do

      jb = 0
      do j=1,nO(2) - nFC(2)
        do b=nCVS(2)+1,nBas-nO(2)
          jb = jb + 1
          f(ia,ixyz) = f(ia,ixyz) + dipole_int_bb(occupations(j,2),virtuals(b,2),ixyz)*XpY(ia,nSa+jb)
        end do
      end do

    end do
  end do

  do ia=1,maxS
    os(ia) = 2d0/3d0*Omega(ia)*sum(f(ia,:)**2)
  end do
    
  write(*,*) '---------------------------------------------------------------'
  write(*,*) '                Transition dipole moment (au)                  '
  write(*,*) '---------------------------------------------------------------'
  write(*,'(A3,5A12)') '#','X','Y','Z','dip. str.','osc. str.'
  write(*,*) '---------------------------------------------------------------'
  do ia=1,maxS
    write(*,'(I3,5F12.6)') ia,(f(ia,ixyz),ixyz=1,ncart),sum(f(ia,:)**2),os(ia)
  end do
  write(*,*) '---------------------------------------------------------------'
  write(*,*)

end subroutine 
