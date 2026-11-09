subroutine RGF2_QP_graph(doSRG,eta,flow,nBas,nC,nO,nV,nR,eHF,ERI,eGFlin,eOld,eGF,Z)

! Compute the graphical solution of the GF2 QP equation

  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: eta
  double precision,intent(in)   :: flow
  logical,intent(in)            :: doSRG
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
  double precision,external     :: RGF2_SigC,RGF2_dSigC,RGF2_SRG_SigC,RGF2_SRG_dSigC
  double precision              :: SigC,dSigC
  double precision              :: f,df
  double precision              :: w
  logical                       :: converged
  
! Output variables

  double precision,intent(out)  :: eGF(nBas)
  double precision,intent(inout):: Z(nBas)

! Run Newton's algorithm to find the root
 
  write(*,*)'-----------------------------------------------------'
  write(*,'(A5,1X,A3,1X,A15,1X,A15,1X,A10)') 'Orb.','It.','e_GFlin (eV)','e_GF (eV)','Z'
  write(*,*)'-----------------------------------------------------'

  eGF(:) = eGFlin(:)
  do p=nC+1,nBas-nR

    w = eGFlin(p)
    converged = .false.

    do nIt=1,maxIt
      if(.not.ieee_is_finite(w)) exit

      if(doSRG) then
        SigC  = RGF2_SRG_SigC(p,w,flow,nBas,nC,nO,nV,nR,eOld,ERI)
        dSigC = RGF2_SRG_dSigC(p,w,flow,nBas,nC,nO,nV,nR,eOld,ERI)
      else
        SigC  = RGF2_SigC(p,w,eta,nBas,nC,nO,nV,nR,eOld,ERI)
        dSigC = RGF2_dSigC(p,w,eta,nBas,nC,nO,nV,nR,eOld,ERI)
      end if

      f  = w - eHF(p) - SigC
      df = 1d0 - dSigC
      if(.not.ieee_is_finite(f) .or. .not.ieee_is_finite(df)) exit
      if(abs(df) <= epsilon(1d0)) exit

      if(abs(f) <= thresh) then
        converged = .true.
        exit
      end if

      w = w - f/df

    end do

    if(.not.converged) then

      nIt = min(nIt,maxIt)
      write(*,'(I5,1X,I3,1X,F15.9,1X,F15.9,1X,F10.6,1X,A12)') p,nIt,eGFlin(p)*HaToeV,eGF(p)*HaToeV,Z(p),'Cvg Failed!'

    else

      eGF(p) = w
      Z(p)   = 1d0/df

      write(*,'(I5,1X,I3,1X,F15.9,1X,F15.9,1X,F10.6)') p,nIt,eGFlin(p)*HaToeV,eGF(p)*HaToeV,Z(p)

    end if

  end do

end subroutine 
