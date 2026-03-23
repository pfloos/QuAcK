
subroutine G_AO_GHF_itau(nBas2,nO,tau,G0itau,cGHF,e_GHF)

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas2
  integer,intent(in)            :: nO

  double precision,intent(in)   :: tau
  double precision,intent(in)   :: cGHF(nBas2,nBas2)
  double precision,intent(in)   :: e_GHF(nBas2)

! Local variables

  integer                       :: iorb

  complex*16,allocatable        :: Gtmp(:,:)

! Output variables
  complex*16,intent(out)        :: G0itau(nBas2,nBas2)
  
!------------------------------------------------------------------------
! Build G<(i tau) and G>(i tau)
!------------------------------------------------------------------------

  allocate(Gtmp(nBas2,nBas2)) 
  Gtmp=czero
  
  if(tau>0d0) then ! G
   do iorb=1,nO
     Gtmp(iorb,iorb) =  im*Exp(e_GHF(iorb)*tau)
   enddo
  else             ! G
   do iorb=nO+1,nBas2
     Gtmp(iorb,iorb) = -im*Exp(e_GHF(iorb)*tau)
   enddo
  endif 

  G0itau=matmul(matmul(cGHF,Gtmp),transpose(cGHF))
  
  deallocate(Gtmp)

end subroutine

