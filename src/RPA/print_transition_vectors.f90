subroutine print_transition_vectors(spin_allowed,nBas,nC,nO,nV,nR,nS,dipole_int,Omega,XpY,XmY)

! Print transition vectors for linear response calculation

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: spin_allowed
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
  integer,parameter             :: maxS = 10
  double precision              :: norm
  double precision,parameter    :: thres_vec = 0.1d0
  double precision,allocatable  :: X(:)
  double precision,allocatable  :: Y(:)
  double precision,allocatable  :: f(:,:)
  double precision,allocatable  :: os(:)

! Memory allocation

  allocate(X(nS),Y(nS),f(nS,ncart),os(nS))

! Initialization
   
  f(:,:) = 0d0
  os(:)  = 0d0

! Compute dipole moments and oscillator strengths

  if(spin_allowed) then

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

  end if
  
! Print details about excitations

  do ia=1,min(nS,maxS)

    X(:) = 0.5d0*(XpY(ia,:) + XmY(ia,:))
    Y(:) = 0.5d0*(XpY(ia,:) - XmY(ia,:))

    print*,'---------------------------------------------'
    write(*,'(A15,I3,A2,F10.6,A3,A6,F6.4,A1)') ' Excitation n. ',ia,': ',Omega(ia)*HaToeV,' eV',' (f = ',os(ia),')'
    print*,'---------------------------------------------'

    jb = 0
    do j=nC+1,nO
      do b=nO+1,nBas-nR
        jb = jb + 1
        if(abs(X(jb)) > thres_vec) write(*,'(I3,A4,I3,A3,F10.6)') j,' -> ',b,' = ',X(jb)/sqrt(2d0)
      end do
    end do
 
    jb = 0
    do j=nC+1,nO
      do b=nO+1,nBas-nR
        jb = jb + 1
        if(abs(Y(jb)) > thres_vec) write(*,'(I3,A4,I3,A3,F10.6)') j,' <- ',b,' = ',Y(jb)/sqrt(2d0)
      end do
    end do
   write(*,*)

   print*,'<S**2> = ',2d0*sum(X(:)**2 + Y(:)**2)

  end do

! Thomas-Reiche-Kuhn sum rule

  write(*,'(A30,F10.6)') 'Thomas-Reiche-Kuhn sum rule = ',sum(os(:))
  write(*,*)

end subroutine print_transition_vectors
