subroutine print_unrestricted_transition_vectors(spin_allowed,nBas,nC,nO,nV,nR,nS,nSt,dipole_int,Omega,XpY,XmY)

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
  integer,intent(in)            :: nSt
  double precision              :: dipole_int(nBas,nBas,ncart,nspin)
  double precision,intent(in)   :: Omega(nSt)
  double precision,intent(in)   :: XpY(nSt,nSt)
  double precision,intent(in)   :: XmY(nSt,nSt)

! Local variables

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

! Compute dipole moments and oscillator strengths

   
  f(:,:) = 0d0
  if(spin_allowed) then

    do ispin=1,nspin
      do ia=1,nSt
        do ixyz=1,ncart
          jb = 0
          do j=nC(ispin)+1,nO(ispin)
            do b=nO(ispin)+1,nBas-nR(ispin)
              jb = jb + 1
              f(ia,ixyz) = f(ia,ixyz) + dipole_int(j,b,ixyz,ispin)*XpY(ia,jb)
            end do
          end do
        end do
      end do
    end do
 
    write(*,*) '----------------'
    write(*,*) ' Dipole moments '
    write(*,*) '----------------'
    call matout(nSt,ncart,f(:,:))
    write(*,*)
 
    do ia=1,nSt
      os(ia) = 2d0/3d0*Omega(ia)*sum(f(ia,:)**2)
    end do
    
    write(*,*) '----------------------'
    write(*,*) ' Oscillator strengths '
    write(*,*) '----------------------'
    call matout(nSt,1,os(:))
    write(*,*)

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
        if(abs(X(jb)) > thres_vec) write(*,'(I3,A4,I3,A3,F10.6)') j,' -> ',b,' = ',X(jb)/sqrt(2d0)
      end do
    end do
 
    jb = 0
    do j=nC(1)+1,nO(1)
      do b=nO(1)+1,nBas-nR(1)
        jb = jb + 1
        if(abs(Y(jb)) > thres_vec) write(*,'(I3,A4,I3,A3,F10.6)') j,' <- ',b,' = ',Y(jb)/sqrt(2d0)
      end do
    end do
   write(*,*)

    ! Spin-down transitions

    jb = 0
    do j=nC(2)+1,nO(2)
      do b=nO(2)+1,nBas-nR(2)
        jb = jb + 1
        if(abs(X(jb)) > thres_vec) write(*,'(I3,A4,I3,A3,F10.6)') j,' -> ',b,' = ',X(jb)/sqrt(2d0)
      end do
    end do
 
    jb = 0
    do j=nC(2)+1,nO(2)
      do b=nO(2)+1,nBas-nR(2)
        jb = jb + 1
        if(abs(Y(jb)) > thres_vec) write(*,'(I3,A4,I3,A3,F10.6)') j,' <- ',b,' = ',Y(jb)/sqrt(2d0)
      end do
    end do
   write(*,*)

  end do

  write(*,'(A30,F10.6)') 'Thomas-Reiche-Kuhn sum rule = ',sum(os(:))
  write(*,*)

end subroutine print_unrestricted_transition_vectors
