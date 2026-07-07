subroutine G_AO_GHFB_w(nBas2,nOrb2,nOrb4,eta,Gen_cHFB,eGHFB,wcoord,Mat1,Mat2,Mat3,Mat4,G_AO)

! G(i w)

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas2
  integer,intent(in)            :: nOrb2
  integer,intent(in)            :: nOrb4

  double precision,intent(in)   :: eta
  double precision,intent(in)   :: Gen_cHFB(nBas2,nOrb2)
  double precision,intent(in)   :: Mat1(nOrb2,nOrb2)
  double precision,intent(in)   :: Mat2(nOrb2,nOrb2)
  double precision,intent(in)   :: Mat3(nOrb2,nOrb2)
  double precision,intent(in)   :: Mat4(nOrb2,nOrb2)


  complex*16,intent(in)         :: wcoord
  complex*16,allocatable        :: Gtmp(:,:)

! Local variables

  integer                       :: iorb

! Output variables
  double precision,intent(inout):: eGHFB(nOrb4)
  complex*16,intent(out)        :: G_AO(nBas2,nBas2)
  
!--------------------------
! Build G(i w) in AO basis
!--------------------------

! write(*,*)     
! write(*,*)'**************'
! write(*,*)'* HFB G(i w) *'
! write(*,*)'**************'
! write(*,*)

 allocate(Gtmp(nOrb2,nOrb2))
 Gtmp(:,:) = czero
  
 do iorb=1,nOrb2
  Gtmp(:,:) = Gtmp(:,:) + 1d0/(wcoord-eGHFB(iorb)-im*eta)           &
            * matmul(Mat1(:,iorb:iorb),transpose(Mat2(:,iorb:iorb)))
  Gtmp(:,:) = Gtmp(:,:) + 1d0/(wcoord+eGHFB(iorb)+im*eta)           &
            * matmul(Mat3(:,iorb:iorb),transpose(Mat4(:,iorb:iorb)))
 enddo

 G_AO=matmul(matmul(Gen_cHFB,Gtmp),transpose(Gen_cHFB))

 ! deallocate dyn arrays
 deallocate(Gtmp)

end subroutine

