subroutine GW_QP_graph(eta,nBas,nC,nO,nV,nR,nS,eHF,Om,rho,eGWlin,eGW,Z)

! Compute the graphical solution of the QP equation

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
  double precision,intent(in)   :: Om(nS)
  double precision,intent(in)   :: rho(nBas,nBas,nS)

  double precision,intent(in)   :: eGWlin(nBas)

! Local variables

  integer                       :: p
  integer                       :: nIt
  integer,parameter             :: maxIt = 64
  double precision,parameter    :: thresh = 1d-6
  double precision,external     :: GW_SigC,GW_dSigC
  double precision              :: SigC,dSigC
  double precision              :: f,df
  double precision              :: w

! Output variables

  double precision,intent(out)  :: eGW(nBas)
  double precision,intent(out)  :: Z(nBas)

! Run Newton's algorithm to find the root
 
  do p=nC+1,nBas-nR

    w = eGWlin(p)
    nIt = 0
    f = 1d0
    
    do while (abs(f) > thresh .and. nIt < maxIt)
    
      nIt = nIt + 1

      SigC  = GW_SigC(p,w,eta,nBas,nC,nO,nV,nR,nS,eGWlin,Om,rho)
      dSigC = GW_dSigC(p,w,eta,nBas,nC,nO,nV,nR,nS,eGWlin,Om,rho)
      f  = w - eHF(p) - SigC
      df = 1d0/(1d0 - dSigC)
    
      w = w - df*f
    
    end do
 
    if(nIt == maxIt) then 

      write(*,*) 'Newton root search has not converged!'
      eGW(p) = eGWlin(p)

    else

      eGW(p) = w
      Z(p)   = df


      write(*,*)'-------------------------------------------------------------------------------'
      write(*,'(A5,1X,A3,1X,A15,1X,A10)') 'Orb.','It.','e_QP (eV)','Z'
      write(*,'(I5,1X,I3,1X,F15.9,1X,F10.6)') p,nIt,eGW(p)*HaToeV,Z(p)
      write(*,*)'-------------------------------------------------------------------------------'
      write(*,*)

    end if
          
  end do

end subroutine 
