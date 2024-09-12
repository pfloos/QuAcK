subroutine RGTpp_QP_graph(eta,nBas,nC,nO,nV,nR,nOOs,nVVs,nOOt,nVVt,eHF,Om1s,rho1s,Om2s,rho2s, & 
                         Om1t,rho1t,Om2t,rho2t,eGTlin,eOld,eGT,Z)

! Compute the graphical solution of the QP equation

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nOOs,nOOt
  integer,intent(in)            :: nVVs,nVVt
  
  double precision,intent(in)   :: eta
  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: Om1s(nVVs),Om1t(nVVt)
  double precision,intent(in)   :: rho1s(nBas,nBas,nVVs),rho1t(nBas,nBas,nVVt)
  double precision,intent(in)   :: Om2s(nOOs),Om2t(nOOt)
  double precision,intent(in)   :: rho2s(nBas,nBas,nOOs),rho2t(nBas,nBas,nOOt)

  double precision,intent(in)   :: eGTlin(nBas)
  double precision,intent(in)   :: eOld(nBas)
  
! Local variables

  integer                       :: p
  integer                       :: nIt
  integer,parameter             :: maxIt = 64
  double precision,parameter    :: thresh = 1d-6
  double precision,external     :: RGTpp_SigC,RGTpp_dSigC
  double precision              :: SigC,dSigC
  double precision              :: f,df
  double precision              :: w
  
! Output variables

  double precision,intent(out)  :: eGT(nBas)
  double precision,intent(out)  :: Z(nBas)

! Run Newton's algorithm to find the root

  write(*,*)'-----------------------------------------------------'
  write(*,'(A5,1X,A3,1X,A15,1X,A15,1X,A10)') 'Orb.','It.','e_GTpplin (eV)','e_GTpp (eV)','Z'
  write(*,*)'-----------------------------------------------------'

  do p=nC+1,nBas-nR

    w = eGTlin(p)
    nIt = 0
    f = 1d0

    do while (abs(f) > thresh .and. nIt < maxIt)

      nIt = nIt + 1

      SigC  = RGTpp_SigC(p,w,eta,nBas,nC,nO,nV,nR,nOOs,nVVs,nOOt,nVVt,eOld,Om1s,rho1s,Om2s,rho2s,Om1t,rho1t,Om2t,rho2t)
      dSigC = RGTpp_dSigC(p,w,eta,nBas,nC,nO,nV,nR,nOOs,nVVs,nOOt,nVVt,eOld,Om1s,rho1s,Om2s,rho2s,Om1t,rho1t,Om2t,rho2t)
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
