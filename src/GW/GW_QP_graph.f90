subroutine GW_QP_graph(eta,nBas,nC,nO,nV,nR,nS,eHF,Om,rho,eGWlin,eGW)

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
  double precision              :: sigC,dsigC
  double precision              :: f,df
  double precision              :: w

! Output variables

  double precision,intent(out)  :: eGW(nBas)

! Run Newton's algorithm to find the root
 
  do p=nC+1,nBas-nR

    write(*,*) '-----------------'
    write(*,'(A10,I3)') 'Orbital ',p
    write(*,*) '-----------------'

    w = eGWlin(p)
    nIt = 0
    f = 1d0
    write(*,'(A3,I3,A1,1X,3F15.9)') 'It.',nIt,':',w*HaToeV,f
    
    do while (abs(f) > thresh .and. nIt < maxIt)
    
      nIt = nIt + 1

      sigC  = GW_SigC(p,w,eta,nBas,nC,nO,nV,nR,nS,eGWlin,Om,rho)
      dsigC = GW_dSigC(p,w,eta,nBas,nC,nO,nV,nR,nS,eGWlin,Om,rho)
      f  = w - eHF(p) - SigC
      df = 1d0 - dsigC
    
      w = w - f/df

      write(*,'(A3,I3,A1,1X,3F15.9)') 'It.',nIt,':',w*HaToeV,f,sigC
    
    
    end do
 
    if(nIt == maxIt) then 

      write(*,*) 'Newton root search has not converged!'

    else

      eGW(p) = w

      write(*,'(A32,F16.10)')   'Quasiparticle energy (eV)   ',eGW(p)*HaToeV
      write(*,*)

    end if
          
  end do

end subroutine 
