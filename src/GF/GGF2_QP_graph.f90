subroutine GGF2_QP_graph(eta,nBas,nC,nO,nV,nR,eHF,ERI,eGFlin,eOld,eGF,Z)

! Compute the graphical solution of the GF2 QP equation

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: eGFlin(nBas)
  double precision,intent(in)   :: eOld(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: p
  integer                       :: nIt
  integer,parameter             :: maxIt = 64
  double precision,parameter    :: thresh = 1d-6
  double precision,external     :: GGF2_SigC,GGF2_dSigC
  double precision              :: SigC,dSigC
  double precision              :: f,df
  double precision              :: w
  
! Output variables

  double precision,intent(out)  :: eGF(nBas)
  double precision,intent(out)  :: Z(nBas)

! Run Newton's algorithm to find the root
 
  write(*,*)'-----------------------------------------------------'
  write(*,'(A5,1X,A3,1X,A15,1X,A15,1X,A10)') 'Orb.','It.','e_GFlin (eV)','e_GF (eV)','Z'
  write(*,*)'-----------------------------------------------------'

  do p=nC+1,nBas-nR

    w = eGFlin(p)
    nIt = 0
    f = 1d0
    
    do while (abs(f) > thresh .and. nIt < maxIt)
    
      nIt = nIt + 1

      SigC  = GGF2_SigC(p,w,eta,nBas,nC,nO,nV,nR,eOld,ERI)
      dSigC = GGF2_dSigC(p,w,eta,nBas,nC,nO,nV,nR,eOld,ERI)
      f  = w - eHF(p) - SigC
      df = 1d0/(1d0 - dSigC)
    
      w = w - df*f
    
    end do
 
    if(nIt == maxIt) then 

      eGF(p) = eGFlin(p)
      write(*,'(I5,1X,I3,1X,F15.9,1X,F15.9,1X,F10.6,1X,A12)') p,nIt,eGFlin(p)*HaToeV,eGF(p)*HaToeV,Z(p),'Cvg Failed!'

    else

      eGF(p) = w
      Z(p)   = df

      write(*,'(I5,1X,I3,1X,F15.9,1X,F15.9,1X,F10.6)') p,nIt,eGFlin(p)*HaToeV,eGF(p)*HaToeV,Z(p)

    end if

  end do

end subroutine 
