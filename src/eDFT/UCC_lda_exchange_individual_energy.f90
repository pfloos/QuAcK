subroutine UCC_lda_exchange_individual_energy(nEns,wEns,aCC_w1,aCC_w2,nGrid,weight,rhow,rho,Cx_choice,doNcentered,kappa,Ex)

! Compute the unrestricted version of the curvature-corrected exchange functional

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
  double precision,intent(in)   :: rho(nGrid)
  integer,intent(in)            :: Cx_choice
  logical,intent(in)            :: doNcentered
  double precision,intent(in)   :: kappa
 
! Local variables

  integer                       :: iG
  double precision              :: r,rI
  double precision              :: e_p,dedr
  double precision              :: Exrr,ExrI,ExrrI

  double precision              :: a1,b1,c1,w1
  double precision              :: a2,b2,c2,w2
  double precision              :: Fx1,Fx2,Cx

! Output variables

  double precision,intent(out)  :: Ex

! External variable

  double precision,external     :: electron_number


! Parameters for N -> N-1

  a1 = aCC_w1(1)
  b1 = aCC_w1(2)
  c1 = aCC_w1(3)

! Parameters for N -> N+1

  a2 = aCC_w2(1)
  b2 = aCC_w2(2)
  c2 = aCC_w2(3)

  if(doNcentered) then

     w1 = wEns(2)
     Fx1 = 1d0 - w1*(2d0 - w1)*(a1 + b1*(w1 - 1d0) + c1*(w1 - 1d0)**2)

     w2 = wEns(3)
     Fx2 = 1d0 - w2*(2d0 - w2)*(a2 + b2*(w2 - 1d0) + c2*(w2 - 1d0)**2)

  else


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

  Ex = 0d0
  Exrr  = 0d0
  ExrI  = 0d0
  ExrrI = 0d0


  do iG=1,nGrid

    r  = max(0d0,rhow(iG))
    rI = max(0d0,rho(iG))

    if(r > threshold) then

      e_p  =         Cx*r**(1d0/3d0)
      dedr = 1d0/3d0*Cx*r**(-2d0/3d0)
     
      Exrr = Exrr - weight(iG)*dedr*r*r

      if(rI > threshold) then

        ExrI  = ExrI  + weight(iG)*e_p*rI
        ExrrI = ExrrI + weight(iG)*dedr*r*rI

      endif

    endif

  enddo

! De-scaling for N-centered ensemble  

  if(doNcentered) then 
   
    Exrr  = kappa*Exrr
    ExrI  = kappa*ExrI

  endif

  Ex = Exrr + ExrI + ExrrI 


end subroutine UCC_lda_exchange_individual_energy
