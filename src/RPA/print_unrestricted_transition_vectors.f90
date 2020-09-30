subroutine print_unrestricted_transition_vectors(spin_allowed,nBas,nC,nO,nV,nR,nS,nSa,nSb,nSt,dipole_int_aa,dipole_int_bb, & 
                                                 Omega,XpY,XmY)

! Print transition vectors for linear response calculation

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: spin_allowed
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC(nspin)
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nV(nspin)
  integer,intent(in)            :: nR(nspin)
  integer,intent(in)            :: nS(nspin)
  integer,intent(in)            :: nSa
  integer,intent(in)            :: nSb
  integer,intent(in)            :: nSt
  double precision              :: dipole_int_aa(nBas,nBas,ncart)
  double precision              :: dipole_int_bb(nBas,nBas,ncart)
  double precision,intent(in)   :: Omega(nSt)
  double precision,intent(in)   :: XpY(nSt,nSt)
  double precision,intent(in)   :: XmY(nSt,nSt)

! Local variables

  logical                       :: debug = .false.
  integer                       :: ia,jb,i,j,a,b
  integer                       :: ixyz
  integer                       :: ispin
  integer,parameter             :: maxS = 10
  double precision              :: norm
  double precision,parameter    :: thres_vec = 0.1d0
  double precision,allocatable  :: X(:)
  double precision,allocatable  :: Y(:)
  double precision,allocatable  :: f(:,:)
  double precision,allocatable  :: os(:)

! Memory allocation

  allocate(X(nSt),Y(nSt),f(nSt,ncart),os(nSt))

! Initialization
   
  f(:,:) = 0d0
  os(:)  = 0d0

! Compute dipole moments and oscillator strengths


  if(spin_allowed) then

    do ia=1,nSt
      do ixyz=1,ncart

        jb = 0
        do j=nC(1)+1,nO(1)
          do b=nO(1)+1,nBas-nR(1)
            jb = jb + 1
            f(ia,ixyz) = f(ia,ixyz) + dipole_int_aa(j,b,ixyz)*XpY(ia,jb)
          end do
        end do

        jb = 0
        do j=nC(2)+1,nO(2)
          do b=nO(2)+1,nBas-nR(2)
            jb = jb + 1
            f(ia,ixyz) = f(ia,ixyz) + dipole_int_bb(j,b,ixyz)*XpY(ia,nSa+jb)
          end do
        end do

      end do
    end do
 
    do ia=1,nSt
      os(ia) = 2d0/3d0*Omega(ia)*sum(f(ia,:)**2)
    end do
      
    if(debug) then

      write(*,*) '----------------'
      write(*,*) ' Dipole moments '
      write(*,*) '----------------'
      call matout(nSt,ncart,f(:,:))
      write(*,*)
  
      write(*,*) '----------------------'
      write(*,*) ' Oscillator strengths '
      write(*,*) '----------------------'
      call matout(nSt,1,os(:))
      write(*,*)

    end if

  end if
  
! Print details about excitations

  do ia=1,min(nSt,maxS)

    X(:) = 0.5d0*(XpY(ia,:) + XmY(ia,:))
    Y(:) = 0.5d0*(XpY(ia,:) - XmY(ia,:))

    print*,'---------------------------------------------'
    write(*,'(A15,I3,A2,F10.6,A3,A6,F6.4,A1)') ' Excitation n. ',ia,': ',Omega(ia)*HaToeV,' eV',' (f = ',os(ia),')'
    print*,'---------------------------------------------'

    ! Spin-up transitions

    jb = 0
    do j=nC(1)+1,nO(1)
      do b=nO(1)+1,nBas-nR(1)
        jb = jb + 1
        if(abs(X(jb)) > thres_vec) write(*,'(I3,A5,I3,A4,F10.6)') j,'A -> ',b,'A = ',X(jb)/sqrt(2d0)
      end do
    end do
 
    jb = 0
    do j=nC(1)+1,nO(1)
      do b=nO(1)+1,nBas-nR(1)
        jb = jb + 1
        if(abs(Y(jb)) > thres_vec) write(*,'(I3,A5,I3,A4,F10.6)') j,'A <- ',b,'A = ',Y(jb)/sqrt(2d0)
      end do
    end do

    ! Spin-down transitions

    jb = 0
    do j=nC(2)+1,nO(2)
      do b=nO(2)+1,nBas-nR(2)
        jb = jb + 1
        if(abs(X(jb)) > thres_vec) write(*,'(I3,A5,I3,A4,F10.6)') j,'B -> ',b,'B = ',X(jb)/sqrt(2d0)
      end do
    end do
 
    jb = 0
    do j=nC(2)+1,nO(2)
      do b=nO(2)+1,nBas-nR(2)
        jb = jb + 1
        if(abs(Y(jb)) > thres_vec) write(*,'(I3,A5,I3,A4,F10.6)') j,'B <- ',b,'B = ',Y(jb)/sqrt(2d0)
      end do
    end do
   write(*,*)

  end do

  write(*,'(A30,F10.6)') 'Thomas-Reiche-Kuhn sum rule = ',sum(os(:))
  write(*,*)

end subroutine print_unrestricted_transition_vectors
