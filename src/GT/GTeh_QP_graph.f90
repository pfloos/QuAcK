subroutine GTeh_QP_graph(eta,nBas,nC,nO,nV,nR,nS,eHF,Om,rhoL,rhoR,eGTlin,eGT,Z)

  implicit none
  include 'parameters.h'

! Iput variables
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  
  double precision,intent(in)   :: eta
  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: Om(nS)
  double precision,intent(in)   :: rhoL(nBas,nBas,nS)
  double precision,intent(in)   :: rhoR(nBas,nBas,nS)

  double precision,intent(in)   :: eGTlin(nBas)
  
! Local variables
  integer                       :: p
  integer                       :: nIt
  integer,parameter             :: maxIt = 64
  double precision,parameter    :: thresh = 1d-6
  double precision,external     :: GTeh_SigC,GTeh_dSigC
  double precision              :: sigC,dsigC
  double precision              :: f,df
  double precision              :: w
  
! Output variables

  double precision,intent(out)  :: eGT(nBas)
  double precision,intent(out)  :: Z(nBas)

  sigC = 0d0
  dsigC = 0d0

! Run Newton's algorithm to find the root
  do p=nC+1,nBas-nR

    write(*,*) '-----------------'
    write(*,'(A10,I3)') 'Orbital ',p
    write(*,*) '-----------------'

    w = eGTlin(p)
    nIt = 0
    f = 1d0
    write(*,'(A3,I3,A1,1X,3F15.9)') 'It.',nIt,':',w*HaToeV,f

    do while (abs(f) > thresh .and. nIt < maxIt)

       nIt = nIt + 1

       sigC  = GTeh_SigC(p,w,eta,nBas,nC,nO,nV,nR,nS,eGTlin,Om,rhoL,rhoR)
       dsigC = GTeh_dSigC(p,w,eta,nBas,nC,nO,nV,nR,nS,eGTlin,Om,rhoL,rhoR)
       f  = w - eHF(p)  - sigC 
       df = 1d0/(1d0 - dsigC)
    
       w = w - df*f

       write(*,'(A3,I3,A1,1X,3F15.9)') 'It.',nIt,':',w*HaToeV,df,f

    end do

    if(nIt == maxIt) then 

      write(*,*) 'Newton root search has not converged!'
      eGT(p) = eGTlin(p)

    else

      eGT(p) = w
      Z(p)   = df
      write(*,'(A32,F16.10)')   'Quasiparticle energy (eV)   ',eGT(p)*HaToeV
      write(*,*)

    end if

  end do
  
end subroutine 
