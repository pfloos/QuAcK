subroutine print_transition_vectors(nBas,nC,nO,nV,nR,nS,Omega,XpY,XmY)

! Print transition vectors for linear response calculation

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: Omega(nS)
  double precision,intent(in)   :: XpY(nS,nS)
  double precision,intent(in)   :: XmY(nS,nS)

! Local variables

  integer                       :: ia,jb,i,j,a,b
  integer,parameter             :: maxS = 10
  double precision,parameter    :: thres_vec = 0.1d0
  double precision,allocatable  :: X(:)
  double precision,allocatable  :: Y(:)

! Memory allocation

  allocate(X(nS),Y(nS))

  write(*,*)
  do ia=1,min(nS,maxS)

  X(:) = 0.5d0*(XpY(ia,:) + XmY(ia,:))
  Y(:) = 0.5d0*(XpY(ia,:) - XmY(ia,:))

  print*,'--------------------------------'
  write(*,'(A15,I3,A2,F10.6,A3)') ' Excitation n. ',ia,': ',Omega(ia)*HaToeV,' eV'
  print*,'--------------------------------'
  print*,'Resonant vectors'
  print*,'--------------------------------'

    jb = 0
    do j=nC+1,nO
      do b=nO+1,nBas-nR
        jb = jb + 1
        if(abs(X(jb)) > thres_vec) write(*,'(I3,A5,I3,A3,F10.6)') j,' --> ',b,' = ',X(jb)
      end do
    end do
 
    print*,'--------------------------------'
    print*,'Anti-resonant vectors'
    print*,'--------------------------------'
    jb = 0
    do j=nC+1,nO
      do b=nO+1,nBas-nR
        jb = jb + 1
        if(abs(Y(jb)) > thres_vec) write(*,'(I3,A5,I3,A3,F10.6)') j,' --> ',b,' = ',Y(jb)
      end do
    end do
   write(*,*)

  end do

end subroutine print_transition_vectors
