subroutine GGTpp_QP_graph(eta,nBas,nC,nO,nV,nR,nOO,nVV,eHF,Om1,rho1,Om2,rho2,eGTlin,eOld,eGT,Z)

! Compute the graphical solution of the QP equation

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nOO
  integer,intent(in)            :: nVV
  
  double precision,intent(in)   :: eta
  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: Om1(nVV)
  double precision,intent(in)   :: rho1(nBas,nBas,nVV)
  double precision,intent(in)   :: Om2(nOO)
  double precision,intent(in)   :: rho2(nBas,nBas,nOO)

  double precision,intent(in)   :: eGTlin(nBas)
  double precision,intent(in)   :: eOld(nBas)
  
! Local variables

  integer                       :: p
  integer                       :: nIt
  integer,parameter             :: maxIt = 64
  double precision,parameter    :: thresh = 1d-6
  double precision,external     :: GGTpp_SigC,GGTpp_dSigC
  double precision              :: SigC,dSigC
  double precision              :: f,df
  double precision              :: w
  
! Output variables

  double precision,intent(out)  :: eGT(nBas)
  double precision,intent(out)  :: Z(nBas)

! Run Newton's algorithm to find the root

  write(*,*)'-----------------------------------------------------'
  write(*,'(A5,1X,A3,1X,A15,1X,A15,1X,A10)') 'Orb.','It.','e_GTpplin (eV)','e_GTpplin (eV)','Z'
  write(*,*)'-----------------------------------------------------'

  do p=nC+1,nBas-nR

    w = eGTlin(p)
    nIt = 0
    f = 1d0

    do while (abs(f) > thresh .and. nIt < maxIt)

      nIt = nIt + 1

      SigC  = GGTpp_SigC(p,w,eta,nBas,nC,nO,nV,nR,nOO,nVV,eOld,Om1,rho1,Om2,rho2)
      dSigC = GGTpp_dSigC(p,w,eta,nBas,nC,nO,nV,nR,nOO,nVV,eOld,Om1,rho1,Om2,rho2)
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
