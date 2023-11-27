subroutine GW_QP_graph(eta,nBas,nC,nO,nV,nR,nS,eHF,Om,rho,eGWlin,eOld,eGW,Z)

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
  double precision,intent(in)   :: eOld(nBas)

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

  write(*,*)'-----------------------------------------------------'
  write(*,'(A5,1X,A3,1X,A15,1X,A15,1X,A10)') 'Orb.','It.','e_GWlin (eV)','e_GW (eV)','Z'
  write(*,*)'-----------------------------------------------------'

  do p=nC+1,nBas-nR

    w = eGWlin(p)
    nIt = 0
    f = 1d0
    
    do while (abs(f) > thresh .and. nIt < maxIt)
    
      nIt = nIt + 1

      SigC  = GW_SigC(p,w,eta,nBas,nC,nO,nV,nR,nS,eOld,Om,rho)
      dSigC = GW_dSigC(p,w,eta,nBas,nC,nO,nV,nR,nS,eOld,Om,rho)
      f  = w - eHF(p) - SigC
      df = 1d0/(1d0 - dSigC)
    
      w = w - df*f
    
    end do
 
    if(nIt == maxIt) then 

      eGW(p) = eGWlin(p)
      write(*,'(I5,1X,I3,1X,F15.9,1X,F15.9,1X,F10.6,1X,A12)') p,nIt,eGWlin(p)*HaToeV,eGW(p)*HaToeV,Z(p),'Cvg Failed!'

    else

      eGW(p) = w
      Z(p)   = df

      write(*,'(I5,1X,I3,1X,F15.9,1X,F15.9,1X,F10.6)') p,nIt,eGWlin(p)*HaToeV,eGW(p)*HaToeV,Z(p)

    end if
          
  end do
  write(*,*)'-----------------------------------------------------'
  write(*,*)

end subroutine 
