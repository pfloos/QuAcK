subroutine UCC_lda_exchange_derivative_discontinuity(nEns,wEns,aCC_w1,aCC_w2,nGrid,weight,rhow,Cx_choice,doNcentered,ExDD)

! Compute the unrestricted version of the curvature-corrected exchange ensemble derivative

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  double precision,intent(in)   :: aCC_w1(3)
  double precision,intent(in)   :: aCC_w2(3)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rhow(nGrid)
  integer,intent(in)            :: Cx_choice
  integer,intent(in)            :: doNcentered


! Local variables

  integer                       :: iEns,jEns
  integer                       :: iG
  double precision              :: r,alpha
  double precision,allocatable  :: dExdw(:)
  double precision,external     :: Kronecker_delta

  double precision              :: a1,b1,c1,w1
  double precision              :: a2,b2,c2,w2
  double precision              :: dCxdw1,dCxdw2

  double precision              :: nEli,nElw


! Output variables

  double precision,intent(out)  :: ExDD(nEns)

! External variable

  double precision,external     :: electron_number


! Memory allocation

  allocate(dExdw(nEns))


! Parameters for N -> N-1

  a1 = aCC_w1(1)
  b1 = aCC_w1(2)
  c1 = aCC_w1(3)


! Parameters for N -> N+1

  a2 = aCC_w2(1)
  b2 = aCC_w2(2)
  c2 = aCC_w2(3)
 
  
  nElw = electron_number(nGrid,weight,rhow)


! Cx coefficient for unrestricted Slater LDA exchange

  alpha = -(3d0/2d0)*(3d0/(4d0*pi))**(1d0/3d0)

  w1 = wEns(2)
  w2 = wEns(3)

  select case (Cx_choice)
  case(1)
    dCxdw1 = (0.5d0*b1 + (2d0*a1 + 0.5d0*c1)*(w1 - 0.5d0) - (1d0 - w1)*w1*(3d0*b1 + 4d0*c1*(w1 - 0.5d0)))
    dCxdw2 = 0.d0
  case(2)
    dCxdw1 = 0.d0
    dCxdw2 =(0.5d0*b2 + (2d0*a2 + 0.5d0*c2)*(w2 - 0.5d0) - (1d0 - w2)*w2*(3d0*b2 + 4d0*c2*(w2 - 0.5d0)))
  case(3)
    dCxdw1 = (0.5d0*b1 + (2d0*a1 + 0.5d0*c1)*(w1 - 0.5d0) - (1d0 - w1)*w1*(3d0*b1 + 4d0*c1*(w1 - 0.5d0))) &
           * (1d0 - w2*(1d0 - w2)*(a2 + b2*(w2 - 0.5d0) + c2*(w2 - 0.5d0)**2))

    dCxdw2 = (1d0 - w1*(1d0 - w1)*(a1 + b1*(w1 - 0.5d0) + c1*(w1 - 0.5d0)**2))                            &
           * (0.5d0*b2 + (2d0*a2 + 0.5d0*c2)*(w2 - 0.5d0) - (1d0 - w2)*w2*(3d0*b2 + 4d0*c2*(w2 - 0.5d0)))  
  end select


  dCxdw1 = alpha*dCxdw1
  dCxdw2 = alpha*dCxdw2

  dExdw(:) = 0d0

  do iG=1,nGrid
    
    r = max(0d0,rhow(iG))
    
    if(r > threshold) then
 
      dExdw(1) = 0d0
      dExdw(2) = dExdw(2) + weight(iG)*dCxdw1*r**(4d0/3d0)
      dExdw(3) = dExdw(3) + weight(iG)*dCxdw2*r**(4d0/3d0)

    end if
     
  end do 

  ExDD(:) = 0d0
 
  do iEns=1,nEns
    do jEns=2,nEns

      ExDD(iEns) = ExDD(iEns) + (Kronecker_delta(iEns,jEns) - wEns(jEns))*dExdw(jEns)

    end do
  end do

  if (doNcentered .NE. 0) then
    ExDD(2) = ((nElw-1)/nElw)*ExDD(2)
    ExDD(3) = ((nElw+1)/nElw)*ExDD(3)
  end if

end subroutine UCC_lda_exchange_derivative_discontinuity
