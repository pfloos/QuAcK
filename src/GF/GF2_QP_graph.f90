subroutine GF2_QP_graph(eta,nBas,nC,nO,nV,nR,nS,eHF,ERI,eGF)

! Compute the graphical solution of the GF2 QP equation

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas,nC,nO,nV,nR,nS
  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: p
  integer                       :: nIt
  integer,parameter             :: maxIt = 64
  double precision,parameter    :: thresh = 1d-6
  double precision,external     :: GF2_SigC,GF2_dSigC
  double precision              :: sigC,dsigC
  double precision              :: f,df
  double precision              :: w
  
! Output variables

  double precision,intent(out)  :: eGF(nBas)


! Run Newton's algorithm to find the root
 
  do p=nC+1,nBas-nR

    write(*,*) '-----------------'
    write(*,'(A10,I3)') 'Orbital ',p
    write(*,*) '-----------------'

    w = eHF(p)
    nIt = 0
    f = 1d0
    write(*,'(A3,I3,A1,1X,3F15.9)') 'It.',nIt,':',w*HaToeV,f
    
    do while (abs(f) > thresh .and. nIt < maxIt)
    
      nIt = nIt + 1

      sigC  = GF2_SigC(p,w,eta,nBas,nC,nO,nV,nR,nS,eHF,ERI)
      dsigC = GF2_dSigC(p,w,eta,nBas,nC,nO,nV,nR,nS,eHF,ERI)
      f  = w - eHF(p) - sigC
      df = 1d0 - dsigC
    
      w = w - f/df

      write(*,'(A3,I3,A1,1X,3F15.9)') 'It.',nIt,':',w*HaToeV,f,sigC
    
    
    end do
 
    if(nIt == maxIt) then 

      write(*,*) 'Newton root search has not converged!'

    else

      eGF(p) = w

      write(*,'(A32,F16.10)')   'Quasiparticle energy (eV)   ',eGF(p)*HaToeV
      write(*,*)

   end if

end do

end subroutine 
