subroutine G_AO_RHFB_w(nBas,nOrb,nOrb_twice,eta,cHFB,eHFB,wcoord,Mat1,Mat2,Mat3,Mat4,G_AO)

! G(i w)

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nOrb_twice

  double precision,intent(in)   :: eta
  double precision,intent(in)   :: cHFB(nBas,nOrb)
  double precision,intent(in)   :: Mat1(nOrb,nOrb)
  double precision,intent(in)   :: Mat2(nOrb,nOrb)
  double precision,intent(in)   :: Mat3(nOrb,nOrb)
  double precision,intent(in)   :: Mat4(nOrb,nOrb)


  complex*16,intent(in)         :: wcoord
  complex*16,allocatable        :: Gtmp(:,:)

! Local variables

  integer                       :: iorb

! Output variables
  double precision,intent(inout):: eHFB(nOrb_twice)
  complex*16,intent(out)        :: G_AO(nBas,nBas)
  
!--------------------------
! Build G(i w) in AO basis
!--------------------------

! write(*,*)     
! write(*,*)'**************'
! write(*,*)'* HFB G(i w) *'
! write(*,*)'**************'
! write(*,*)

 allocate(Gtmp(nOrb,nOrb))
 Gtmp(:,:) = czero
  
 do iorb=1,nOrb
  Gtmp(:,:) = Gtmp(:,:) + 1d0/(wcoord-eHFB(iorb)-im*eta)           &
            * matmul(Mat1(:,iorb:iorb),transpose(Mat2(:,iorb:iorb)))
  Gtmp(:,:) = Gtmp(:,:) + 1d0/(wcoord+eHFB(iorb)+im*eta)           &
            * matmul(Mat3(:,iorb:iorb),transpose(Mat4(:,iorb:iorb)))
 enddo

 G_AO=matmul(matmul(cHFB,Gtmp),transpose(cHFB))

 ! deallocate dyn arrays
 deallocate(Gtmp)

end subroutine

