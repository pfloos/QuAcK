subroutine R_SOSEX_QP_graph(doSRG,eta,flow,nBas,nC,nO,nV,nR,nS,eHF,Om,rhoL,rhoR,eQPlin,eOld,eQP,Z)

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

  logical,intent(in)            :: doSRG
  double precision,intent(in)   :: eta
  double precision,intent(in)   :: flow
  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: Om(nS)
  double precision,intent(in)   :: rhoL(nBas,nBas,nS)
  double precision,intent(in)   :: rhoR(nBas,nBas,nS)

  double precision,intent(in)   :: eQPlin(nBas)
  double precision,intent(in)   :: eOld(nBas)

! Local variables

  integer                       :: p
  integer                       :: nIt
  integer,parameter             :: maxIt = 64
  double precision,parameter    :: thresh = 1d-6
  double precision,external     :: RGW_Re_SigC,RGW_Re_dSigC
  double precision,external     :: RGW_SRG_Re_SigC,RGW_SRG_Re_dSigC
  double precision              :: SigC,dSigC
  double precision              :: f,df
  double precision              :: w

! Output variables

  double precision,intent(out)  :: eQP(nBas)
  double precision,intent(out)  :: Z(nBas)

! Run Newton's algorithm to find the root

  write(*,*)'-----------------------------------------------------'
  write(*,'(A5,1X,A3,1X,A15,1X,A15,1X,A10)') 'Orb.','It.','e_QPlin (eV)','e_QP (eV)','Z'
  write(*,*)'-----------------------------------------------------'

  do p=nC+1,nBas-nR

    w = eQPlin(p)
    nIt = 0
    f = 1d0
    
    do while (abs(f) > thresh .and. nIt < maxIt)
    
      nIt = nIt + 1

      if(doSRG) then

        SigC  = RGW_SRG_Re_SigC(p,w,flow,nBas,nC,nO,nV,nR,nS,eOld,Om,rhoL,rhoR)
        dSigC = RGW_SRG_Re_dSigC(p,w,flow,nBas,nC,nO,nV,nR,nS,eOld,Om,rhoL,rhoR)

        else

        SigC  = RGW_Re_SigC(p,w,eta,nBas,nC,nO,nV,nR,nS,eOld,Om,rhoL,rhoR)
        dSigC = RGW_Re_dSigC(p,w,eta,nBas,nC,nO,nV,nR,nS,eOld,Om,rhoL,rhoR)

      end if

      f  = w - eHF(p) - SigC
      df = 1d0/(1d0 - dSigC)
      w = w - df*f
    
    end do
 
    if(nIt == maxIt) then 

      eQP(p) = eQPlin(p)
      write(*,'(I5,1X,I3,1X,F15.9,1X,F15.9,1X,F10.6,1X,A12)') p,nIt,eQPlin(p)*HaToeV,eQP(p)*HaToeV,Z(p),'Cvg Failed!'

    else

      eQP(p) = w
      Z(p)   = df

      write(*,'(I5,1X,I3,1X,F15.9,1X,F15.9,1X,F10.6)') p,nIt,eQPlin(p)*HaToeV,eQP(p)*HaToeV,Z(p)

    end if
          
  end do

  write(*,*)'-----------------------------------------------------'
  write(*,*)

end subroutine 
