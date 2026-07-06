
subroutine G_AO_GHFB_itau(nBas2,nOrb2,nOrb4,tau,G0itau,Gen_cHFB,eGHFB,Mat1,Mat2,Mat3,Mat4)

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas2
  integer,intent(in)            :: nOrb2
  integer,intent(in)            :: nOrb4

  double precision,intent(in)   :: tau
  double precision,intent(in)   :: Gen_cHFB(nBas2,nOrb2)
  double precision,intent(in)   :: eGHFB(nOrb4)
  double precision,intent(in)   :: Mat1(nOrb2,nOrb2)
  double precision,intent(in)   :: Mat2(nOrb2,nOrb2)
  double precision,intent(in)   :: Mat3(nOrb2,nOrb2)
  double precision,intent(in)   :: Mat4(nOrb2,nOrb2)


! Local variables

  integer                       :: iorb

  complex*16,allocatable        :: Gtmp(:,:)

! Output variables
  complex*16,intent(out)        :: G0itau(nBas2,nBas2)
  
!------------------------------------------------------------------------
! Build G<(i tau) and G>(i tau)
!------------------------------------------------------------------------

  allocate(Gtmp(nOrb2,nOrb2)) 
  Gtmp=czero
  
  if(tau>0d0) then ! Gocc
   do iorb=1,nOrb2
     Gtmp(:,:) = Gtmp(:,:) + im*Exp(eGHFB(iorb)*tau)                     &
               * matmul(Mat1(:,iorb:iorb),transpose(Mat2(:,iorb:iorb)))
   enddo
  else             ! Ghole 
   do iorb=1,nOrb2
     Gtmp(:,:) = Gtmp(:,:) - im*Exp(-eGHFB(iorb)*tau)                    &
               * matmul(Mat3(:,iorb:iorb),transpose(Mat4(:,iorb:iorb)))
   enddo
  endif 

  G0itau=matmul(matmul(Gen_cHFB,Gtmp),transpose(Gen_cHFB))
  
  deallocate(Gtmp)

end subroutine

