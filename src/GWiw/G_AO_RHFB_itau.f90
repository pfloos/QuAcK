
subroutine G_AO_RHFB_itau(nBas,nOrb,nOrb_twice,tau,G0itau,cHFB,eHFB,Mat1,Mat2,Mat3,Mat4)

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nOrb_twice

  double precision,intent(in)   :: tau
  double precision,intent(in)   :: cHFB(nBas,nOrb)
  double precision,intent(in)   :: eHFB(nOrb_twice)
  double precision,intent(in)   :: Mat1(nOrb,nOrb)
  double precision,intent(in)   :: Mat2(nOrb,nOrb)
  double precision,intent(in)   :: Mat3(nOrb,nOrb)
  double precision,intent(in)   :: Mat4(nOrb,nOrb)


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
  
  if(tau>0d0) then ! Gocc
   do iorb=1,nOrb
     Gtmp(:,:) = Gtmp(:,:) + im*Exp(eHFB(iorb)*tau)                     &
               * matmul(Mat1(:,iorb:iorb),transpose(Mat2(:,iorb:iorb)))
   enddo
  else             ! Ghole 
   do iorb=1,nOrb
     Gtmp(:,:) = Gtmp(:,:) - im*Exp(-eHFB(iorb)*tau)                    &
               * matmul(Mat3(:,iorb:iorb),transpose(Mat4(:,iorb:iorb)))
   enddo
  endif 

  G0itau=matmul(matmul(cHFB,Gtmp),transpose(cHFB))
  
  deallocate(Gtmp)

end subroutine

