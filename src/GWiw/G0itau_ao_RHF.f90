
subroutine G0itau_ao_RHF(nBas,nOrb,nO,tau,G0itau,cHF,eHF)

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nO

  double precision,intent(in)   :: tau
  double precision,intent(in)   :: cHF(nBas,nOrb)
  double precision,intent(in)   :: eHF(nOrb)

! Local variables

  integer                       :: iorb

  complex*16,allocatable        :: Gtmp(:,:)

! Output variables
  complex*16,intent(out)        :: G0itau(nBas,nBas)
  
!------------------------------------------------------------------------
! Build G<(i tau) and G>(i tau)
!------------------------------------------------------------------------

  allocate(Gtmp(nOrb,nOrb)) 
  Gtmp=czero
  
  if(tau>0d0) then ! G
   do iorb=1,nO
     Gtmp(iorb,iorb) =  im*Exp(eHF(iorb)*tau)
   enddo
  else             ! G
   do iorb=nO+1,nOrb
     Gtmp(iorb,iorb) = -im*Exp(eHF(iorb)*tau)
   enddo
  endif 

  G0itau=matmul(matmul(cHF,Gtmp),transpose(cHF))
  
  deallocate(Gtmp)

end subroutine

