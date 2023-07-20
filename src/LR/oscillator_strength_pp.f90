subroutine oscillator_strength_pp(nBas,nC,nO,nV,nR,nOO,nVV,maxOO,maxVV,dipole_int,Om1,X1,Y1,Om2,X2,Y2,os1,os2)

! Compute linear response

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nOO
  integer,intent(in)            :: nVV
  integer,intent(in)            :: maxOO
  integer,intent(in)            :: maxVV
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)
  double precision,intent(in)   :: Om1(nVV)
  double precision,intent(in)   :: X1(nVV,nVV)
  double precision,intent(in)   :: Y1(nOO,nVV)
  double precision,intent(in)   :: Om2(nOO)
  double precision,intent(in)   :: X2(nVV,nOO)
  double precision,intent(in)   :: Y2(nOO,nOO)
  
! Local variables

  integer                       :: a,b,c,d,ab,cd
  integer                       :: i,j,k,l,ij,kl
  integer                       :: ixyz

  double precision,allocatable  :: f1(:,:)
  double precision,allocatable  :: f2(:,:)

! Output variables

  double precision,intent(out)  :: os1(nVV)
  double precision,intent(out)  :: os2(nOO)

! Memory allocation

  allocate(f1(maxVV,ncart),f2(maxOO,ncart))

! Initialization
   
  f1(:,:) = 0d0
  f2(:,:) = 0d0

! Compute dipole moments and oscillator strengths

  do ab=1,maxVV
    do ixyz=1,ncart

      cd = 0
      do c=nO+1,nBas-nR
        do d=c,nBas-nR
          cd = cd + 1
          f1(ab,ixyz) = f1(ab,ixyz) + dipole_int(c,d,ixyz)*X1(cd,ab)
        end do
      end do

      kl = 0
      do k=nC+1,nO
        do l=k,nO
          kl = kl + 1
          f1(ab,ixyz) = f1(ab,ixyz) + dipole_int(k,l,ixyz)*Y1(kl,ab)
        end do
      end do

    end do
  end do
  f1(:,:) = sqrt(2d0)*f1(:,:)

  do ab=1,maxVV
    os1(ab) = +2d0/3d0*abs(Om1(ab))*sum(f1(ab,:)**2)
  end do


  do ij=1,maxOO
    do ixyz=1,ncart

      cd = 0
      do c=nO+1,nBas-nR
        do d=c,nBas-nR
          cd = cd + 1
          f2(ij,ixyz) = f2(ij,ixyz) + dipole_int(c,d,ixyz)*X2(cd,ij)
        end do
      end do

      kl = 0
      do k=nC+1,nO
        do l=k,nO
          kl = kl + 1
          f2(ij,ixyz) = f2(ij,ixyz) + dipole_int(k,l,ixyz)*Y2(kl,ij)
        end do
      end do

    end do
  end do
  f2(:,:) = sqrt(2d0)*f2(:,:)

  do ij=1,maxOO
    os2(ij) = 2d0/3d0*abs(Om2(ij))*sum(f2(ij,:)**2)
  end do
 
  write(*,*) '---------------------------------------------------------------'
  write(*,*) '           Transition dipole moment N -> N+2 (au)              '
  write(*,*) '---------------------------------------------------------------'
  write(*,'(A3,5A12)') '#','X','Y','Z','dip. str.','osc. str.'
  write(*,*) '---------------------------------------------------------------'
  do ab=1,maxVV
    write(*,'(I3,5F12.6)') ab,(f1(ab,ixyz),ixyz=1,ncart),sum(f1(ab,:)**2),os1(ab)
  end do
  write(*,*)
  write(*,*) '---------------------------------------------------------------'
  write(*,*) '           Transition dipole moment N -> N-2 (au)              '
  write(*,*) '---------------------------------------------------------------'
  do ij=1,maxOO
    write(*,'(I3,5F12.6)') ij,(f2(ij,ixyz),ixyz=1,ncart),sum(f2(ij,:)**2),os2(ij)
  end do
  write(*,*) '---------------------------------------------------------------'
  write(*,*)

end subroutine 
