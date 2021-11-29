subroutine UCC_lda_exchange_potential(nEns,wEns,nCC,aCC,nGrid,weight,nBas,AO,rho,Fx,Cx_choice,doNcentered)

! Compute the unrestricted version of the curvature-corrected exchange potential

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nCC
  double precision,intent(in)   :: aCC(nCC,nEns-1)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  integer,intent(in)            :: nBas
  double precision,intent(in)   :: AO(nBas,nGrid)
  double precision,intent(in)   :: rho(nGrid)
  integer,intent(in)            :: Cx_choice
  logical,intent(in)            :: doNcentered

! Local variables

  integer                       :: mu,nu,iG
  double precision              :: r,vAO

  double precision              :: a1,b1,c1,d1,w1
  double precision              :: a2,b2,c2,d2,w2
  double precision              :: Fx1,Fx2,Cx

! Output variables

  double precision,intent(out)  :: Fx(nBas,nBas)


! Defining enhancements factor for weight-dependent functionals

  if(doNcentered) then

! Parameters for first state

    a1 = aCC(1,1)
    b1 = aCC(2,1)
    c1 = aCC(3,1)
    d1 = aCC(4,1)

! Parameters for second state

    a2 = aCC(1,2)
    b2 = aCC(2,2)
    c2 = aCC(3,2)
    d2 = aCC(4,2)

    w1 = wEns(2)
    Fx1 = 1d0 + a1*w1 + b1*w1**2 + c1*w1**3 + d1*w1**4

    w2 = wEns(3)
    Fx2 = 1d0 + a2*w2 + b2*w2**2 + c2*w2**3 + d2*w2**4

  else

! Parameters for first state

     a1 = aCC(1,1)
     b1 = aCC(2,1)
     c1 = aCC(3,1)

! Parameters for second state

     a2 = aCC(1,2)
     b2 = aCC(2,2)
     c2 = aCC(3,2)

     w1 = wEns(2)  
     Fx1 = 1d0 - w1*(1d0 - w1)*(a1 + b1*(w1 - 0.5d0) + c1*(w1 - 0.5d0)**2)           

     w2 = wEns(3)  
     Fx2 = 1d0 - w2*(1d0 - w2)*(a2 + b2*(w2 - 0.5d0) + c2*(w2 - 0.5d0)**2) 
 
  endif

  select case (Cx_choice)

    case(1)
      Cx = CxLSDA*Fx1

    case(2)
      Cx = CxLSDA*Fx2

    case(3)
      Cx = CxLSDA*Fx2*Fx1

    case default
      Cx = CxLSDA

  end select

! Compute LDA exchange matrix in the AO basis

  Fx(:,:) = 0d0

  do mu=1,nBas
    do nu=1,nBas
      do iG=1,nGrid

        r = max(0d0,rho(iG))

        if(r > threshold) then

          vAO = weight(iG)*AO(mu,iG)*AO(nu,iG)
          Fx(mu,nu) = Fx(mu,nu) + vAO*4d0/3d0*Cx*r**(1d0/3d0)

        endif

      enddo
    enddo
  enddo

end subroutine UCC_lda_exchange_potential
