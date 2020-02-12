subroutine QP_roots(nBas,nC,nO,nV,nR,nS,eta,eHF,Omega,rho)

! Compute all the roots of the QP equation for each orbital

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: eta
  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: Omega(nS)
  double precision,intent(in)   :: rho(nBas,nBas,nS)

! Local variables

  integer                       :: i,j,a,b,x,jb,g
  integer                       :: nRoot
  integer,parameter             :: nGrid = 10000
  double precision,parameter    :: wmin = -10d0
  double precision,parameter    :: wmax = +10d0
  double precision              :: dw
  double precision              :: eps
  double precision              :: left_sign,right_sign
  double precision              :: left_QP,right_QP
  double precision,allocatable  :: w(:)
  double precision,allocatable  :: SigC(:)
  double precision,allocatable  :: dSigC(:)

  integer                       :: nIt
  integer,parameter             :: maxIt = 64
  double precision,parameter    :: thresh = 1d-6
  double precision,external     :: SigmaC,dSigmaC
  double precision              :: f,df
  double precision              :: s,ds
  double precision              :: om

! Construct grid

  allocate(w(nGrid),SigC(nGrid),dSigC(nGrid))

! Minimum and maximum frequency values

  dw = (wmax - wmin)/dble(nGrid)

  do g=1,nGrid
    w(g) = wmin + dble(g)*dw
  enddo

! Main loop over the orbitals

  do x=nC+1,nBas-nR

  SigC(:)  = 0d0
  dSigC(:) = 0d0

    ! Loop over grid points
 
    do g=1,nGrid

      ! Occupied part of the self-energy and spectral weight

      do i=nC+1,nO
        jb = 0
        do j=nC+1,nO
          do b=nO+1,nBas-nR

            jb = jb + 1
            eps = w(g) - eHF(i) + Omega(jb)
            SigC(g)  = SigC(g)  + 2d0*rho(x,i,jb)**2*eps/(eps**2 + eta**2)
            dSigC(g) = dSigC(g) - 2d0*rho(x,i,jb)**2*(eps**2 - eta**2)/(eps**2 + eta**2)**2

          enddo
        enddo
      enddo
 
      ! Virtual part of the self-energy and spectral weight
 
      do a=nO+1,nBas-nR
        jb = 0
        do j=nC+1,nO
          do b=nO+1,nBas-nR

            jb = jb + 1
            eps = w(g) - eHF(a) - Omega(jb)
            SigC(g)  = SigC(g)  + 2d0*rho(x,a,jb)**2*eps/(eps**2 + eta**2)
            dSigC(g) = dSigC(g) - 2d0*rho(x,a,jb)**2*(eps**2 - eta**2)/(eps**2 + eta**2)**2

          enddo                                               
        enddo                                                 
      enddo                                                   

    enddo

! Find the zeros of the QP equation

    nRoot = 0

    write(*,*) '-----------------'
    write(*,'(A10,I3)') 'Orbital ',x
    write(*,*) '-----------------'
    write(*,*) 
    left_QP   = w(1) - eHF(x) - SigC(1)
    left_sign = sign(1d0,left_QP)

    do g=2,nGrid

      right_QP = w(g) - eHF(x) - SigC(g)
      right_sign = sign(1d0,right_QP)

      if(left_sign /= right_sign .and. left_sign == sign(1d0,-1d0)) then

        nRoot = nRoot + 1
        write(*,'(A20,I6,F10.6,F10.6,I6,F10.6,F10.6,F10.6)') & 
          'root right here!',g-1,w(g-1),left_QP,g,w(g),right_QP,1d0/(1d0-dSigC(g))

        ! Run Newton's algorithm to find the root

          om = w(g-1)
          nIt = 0
          f = 1d0
    
          do while (abs(f) > thresh .and. nIt < maxIt)
    
            nIt = nIt + 1

            s  =  SigmaC(x,om,eta,nBas,nC,nO,nV,nR,nS,eHF,Omega,rho)
            ds = dSigmaC(x,om,eta,nBas,nC,nO,nV,nR,nS,eHF,Omega,rho)
            f  = om - eHF(x) - s
            df = 1d0 - ds
    
            write(*,'(A3,I3,A1,1X,4F15.9)') 'It.',nIt,':',om,f,s,1d0/(1d0-ds)
    
            om = om - f/df
    
          end do

      end if

      left_sign = right_sign
      left_QP   = right_QP
   
    end do
    write(*,*)
    write(*,'(A32,I3,A1,I3)') 'Number of roots for orbital ',x,':',nRoot
    write(*,*)

  end do

end subroutine QP_roots
