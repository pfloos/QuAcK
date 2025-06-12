subroutine complex_cRGF2_QP_graph(doSRG,eta,flow,nBas,nC,nO,nV,nR,Re_eHF,Im_eHF,&
                ERI,Re_eGFlin,Im_eGFlin,Re_eOld,Im_eold,Re_eGF,Im_eGF,Re_Z,Im_Z)

! Compute the graphical solution of the complex GF2 QP equation

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: eta
  double precision,intent(in)   :: flow
  logical, intent(in)           :: doSRG
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  double precision,intent(in)   :: Re_eHF(nBas)
  double precision,intent(in)   :: Im_eHF(nBas)
  double precision,intent(in)   :: Re_eGFlin(nBas)
  double precision,intent(in)   :: Im_eGFlin(nBas)
  double precision,intent(in)   :: Re_eOld(nBas)
  double precision,intent(in)   :: Im_eOld(nBas)
  complex*16,intent(in)         :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: p
  integer                       :: nIt
  integer,parameter             :: maxIt = 64
  double precision,parameter    :: thresh = 1d-6
  double precision              :: Re_SigC,Im_SigC,Re_dSigC,Im_dSigC
  double precision              :: Re_f,Im_f,Re_df,Im_df
  double precision              :: Re_w,Im_w
  
! Output variables

  double precision,intent(out)  :: Re_eGF(nBas),Im_eGF(nBas)
  double precision,intent(out)  :: Re_Z(nBas),Im_Z(nBas)

! Run Newton's algorithm to find the root
 
  write(*,*)'-----------------------------------------------------'
  write(*,'(A5,1X,A3,1X,A15,1X,A15,1X,A10)') 'Orb.','It.','Re(e_GFlin) (eV)','Re(e_GF) (eV)','Re(Z)'
  write(*,'(A5,1X,A3,1X,A15,1X,A15,1X,A10)') 'Orb.','It.','Im(e_GFlin) (eV)','Im(e_GF) (eV)','Im(Z)'
  write(*,*)'-----------------------------------------------------'

  Re_SigC = 0d0
  Im_SigC = 0d0
  Re_dSigC = 0d0
  Im_dSigC = 0d0

  do p=nC+1,nBas-nR

    Re_w = Re_eGFlin(p)
    Im_w = Im_eGFlin(p)
    nIt = 0
    Re_f = 1d0
    Im_f = 1d0
    
    do while (abs(cmplx(Re_f,Im_f,kind=8)) > thresh .and. nIt < maxIt)
    
      nIt = nIt + 1
      if(doSRG) then
        call complex_cRGF2_SRG_SigC_dSigC(flow,p,eta,nBas,nC,nO,nV,nR,Re_w,Im_w,Re_eOld,Im_eOld,ERI,&
              Re_SigC,Im_SigC,Re_dSigC,Im_dSigC)
      else 
        call complex_cRGF2_SigC_dSigC(p,eta,nBas,nC,nO,nV,nR,Re_w,Im_w,Re_eOld,Im_eOld,ERI,&
              Re_SigC,Im_SigC,Re_dSigC,Im_dSigC)
      end if
      Re_f  = Re_w - Re_eHF(p) - Re_SigC 
      Im_f  = Im_w - Im_eHF(p) - Im_SigC
      Re_df = (1d0 - Re_dSigC)/((1d0 - Re_dSigC)**2 + Im_dSigC**2)
      Im_df = Im_dSigC/((1d0 - Re_dSigC)**2 + Im_dSigC**2)
    
      Re_w = Re_w - Re_df*Re_f + Im_df*Im_f
      Im_w = Im_w - Re_f*Im_df - Re_df*Im_f
    
    end do
 
    if(nIt == maxIt) then 

      Re_eGF(p) = Re_eGFlin(p)
      Im_eGF(p) = Im_eGFlin(p)
      write(*,'(I5,1X,I3,1X,F15.9,1X,F15.9,1X,F10.6,1X,A12)') p,nIt,Re_eGFlin(p)*HaToeV,Re_eGF(p)*HaToeV,Re_Z(p),'Cvg Failed!'
      write(*,'(I5,1X,I3,1X,F15.9,1X,F15.9,1X,F10.6,1X,A12)') p,nIt,Im_eGFlin(p)*HaToeV,Im_eGF(p)*HaToeV,Im_Z(p),'Cvg Failed!'

    else

      Re_eGF(p) = Re_w
      Im_eGF(p) = Im_w
      Re_Z(p)   = Re_df
      Im_Z(p)   = Im_df

      write(*,'(I5,1X,I3,1X,F15.9,1X,F15.9,1X,F10.6)') p,nIt,Re_eGFlin(p)*HaToeV,Re_eGF(p)*HaToeV,Re_Z(p)
      write(*,'(I5,1X,I3,1X,F15.9,1X,F15.9,1X,F10.6)') p,nIt,Im_eGFlin(p)*HaToeV,Im_eGF(p)*HaToeV,Im_Z(p)
      !write(*,*) "e_GF - e_HF - \Sigma(e_GF) = ", sqrt(Re_f**2 + Im_f**2)
    write(*,*)'-----------------------------------------------------'
    end if

  end do

end subroutine 
