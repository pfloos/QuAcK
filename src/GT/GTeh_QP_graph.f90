subroutine GTeh_QP_graph(eta,nBas,nC,nO,nV,nR,nS,eHF,Om,rhoL,rhoR,eGTlin,eOld,eGT,Z)

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
  double precision,intent(in)   :: rhoL(nBas,nBas,nS)
  double precision,intent(in)   :: rhoR(nBas,nBas,nS)

  double precision,intent(in)   :: eGTlin(nBas)
  double precision,intent(in)   :: eOld(nBas)
  
! Local variables

  integer                       :: p
  integer                       :: nIt
  integer,parameter             :: maxIt = 64
  double precision,parameter    :: thresh = 1d-6
  double precision,external     :: GTeh_SigC,GTeh_dSigC
  double precision              :: SigC,dSigC
  double precision              :: f,df
  double precision              :: w
  
! Output variables

  double precision,intent(out)  :: eGT(nBas)
  double precision,intent(out)  :: Z(nBas)

! Run Newton's algorithm to find the root

  write(*,*)'-----------------------------------------------------'
  write(*,'(A5,1X,A3,1X,A15,1X,A15,1X,A10)') 'Orb.','It.','e_GTehlin (eV)','e_GTehlin (eV)','Z'
  write(*,*)'-----------------------------------------------------'

  do p=nC+1,nBas-nR

    w = eGTlin(p)
    nIt = 0
    f = 1d0

    do while (abs(f) > thresh .and. nIt < maxIt)

      nIt = nIt + 1

      SigC  = GTeh_SigC(p,w,eta,nBas,nC,nO,nV,nR,nS,eOld,Om,rhoL,rhoR)
      dSigC = GTeh_dSigC(p,w,eta,nBas,nC,nO,nV,nR,nS,eOld,Om,rhoL,rhoR)
      f  = w - eHF(p) - SigC 
      df = 1d0/(1d0 - dSigC)
      w = w - df*f

    end do

    if(nIt == maxIt) then 

      eGT(p) = eGTlin(p)
      write(*,'(I5,1X,I3,1X,F15.9,1X,F15.9,1X,F10.6,1X,A12)') p,nIt,eGTlin(p)*HaToeV,eGT(p)*HaToeV,Z(p),'Cvg Failed!'

    else

      eGT(p) = w
      Z(p)   = df

      write(*,'(I5,1X,I3,1X,F15.9,1X,F15.9,1X,F10.6)') p,nIt,eGTlin(p)*HaToeV,eGT(p)*HaToeV,Z(p)

    end if

  end do

  write(*,*)'-----------------------------------------------------'
  write(*,*)
  
end subroutine 
