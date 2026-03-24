!  subroutine RG0F3(dotest,renormalization,nBas,nOrb,nC,nO,nV,nR,V,e0)

! ! Perform third-order Green function calculation in diagonal approximation

!   implicit none
!   include 'parameters.h'

! ! Input variables

!   logical,intent(in)            :: dotest
!   integer,intent(in)            :: renormalization
!   integer,intent(in)            :: nBas,nOrb,nC,nO,nV,nR
!   double precision,intent(in)   :: e0(nOrb),V(nOrb,nOrb,nOrb,nOrb)

! ! Local variables

!   double precision              :: eps,eps1,eps2
!   double precision,allocatable  :: Sig2(:),SigInf(:),Sig3(:),eGF3(:),eOld(:)
!   double precision,allocatable  :: App(:,:),Bpp(:,:),Cpp(:,:),Dpp(:,:)
!   double precision,allocatable  :: Z(:),X2h1p(:),X1h2p(:),Sig2h1p(:),Sig1h2p(:)

!   integer                       :: i,j,k,l,a,b,c,d,p

! ! Hello world

!   write(*,*)
!   write(*,*)'************************************************'
!   write(*,*)'|    Third-order Green function calculation    |'
!   write(*,*)'************************************************'
!   write(*,*)

! ! Memory allocation

!   allocate(eGF3(nOrb),Sig2(nOrb),SigInf(nOrb),Sig3(nOrb),   &
!            App(nOrb,6),Bpp(nOrb,2),Cpp(nOrb,6),Dpp(nOrb,6), &
!            Z(nOrb),X2h1p(nOrb),X1h2p(nOrb),Sig2h1p(nOrb),Sig1h2p(nOrb))

! !------------------------------------------------------------------------
! ! Compute third-order frequency-independent contribution
! !------------------------------------------------------------------------

!   App(:,:) = 0d0

!   do p=nC+1,nOrb-nR
!     do i=nC+1,nO
!     do j=nC+1,nO
!     do k=nC+1,nO
!       do a=nO+1,nOrb-nR
!       do b=nO+1,nOrb-nR

!           eps1 = e0(j) + e0(i) - e0(a) - e0(b)
!           eps2 = e0(k) + e0(i) - e0(a) - e0(b)

!           App(p,1) = App(p,1) & 
!                    - (2d0*V(p,k,p,j) - V(p,k,j,p))*(2d0*V(j,i,a,b) - V(j,i,b,a))*V(a,b,k,i)/(eps1*eps2)

!       end do
!       end do
!     end do
!     end do
!     end do
!   end do

!   do p=nC+1,nOrb-nR
!     do i=nC+1,nO
!     do j=nC+1,nO
!       do a=nO+1,nOrb-nR
!       do b=nO+1,nOrb-nR
!       do c=nO+1,nOrb-nR

!           eps1 = e0(j) + e0(i) - e0(a) - e0(b)
!           eps2 = e0(j) + e0(i) - e0(a) - e0(c)

!           App(p,2) = App(p,2) & 
!                    + (2d0*V(p,c,p,b) - V(p,c,b,p))*(2d0*V(j,i,a,b) - V(j,i,b,a))*V(j,i,c,a)/(eps1*eps2)

!       end do
!       end do
!       end do
!     end do
!     end do
!   end do
  
!   do p=nC+1,nOrb-nR
!     do i=nC+1,nO
!     do j=nC+1,nO
!       do a=nO+1,nOrb-nR
!       do b=nO+1,nOrb-nR
!       do c=nO+1,nOrb-nR

!           eps1 = e0(j) + e0(i) - e0(a) - e0(b)
!           eps2 = e0(j)         - e0(c)

!           App(p,3) = App(p,3) & 
!                    + (2d0*V(p,c,p,j) - V(p,c,j,p))*(2d0*V(j,i,a,b) - V(j,i,b,a))*V(a,b,c,i)/(eps1*eps2)

!       end do
!       end do
!       end do
!     end do
!     end do
!   end do

!   App(:,4) = App(:,3)

!   do p=nC+1,nOrb-nR
!     do i=nC+1,nO
!     do j=nC+1,nO
!     do k=nC+1,nO
!       do a=nO+1,nOrb-nR
!       do b=nO+1,nOrb-nR

!           eps1 = e0(j) + e0(i) - e0(a) - e0(b)
!           eps2 = e0(k)         - e0(b)

!           App(p,5) = App(p,5) & 
!                    - (2d0*V(p,b,p,k) - V(p,b,k,p))*(2d0*V(j,i,a,b) - V(j,i,b,a))*V(i,j,k,a)/(eps1*eps2)

!       end do
!       end do
!     end do
!     end do
!     end do
!   end do
  
!   App(:,6) = App(:,5)

! ! Frequency-independent part of the third-order self-energy

!   SigInf(:) = App(:,1) + App(:,2) + App(:,3) + App(:,4) + App(:,5) + App(:,6)

! ! Frequency-dependent second-order contribution

!   Bpp(:,:) = 0d0

!   do p=nC+1,nOrb-nR
!     do i=nC+1,nO
!     do j=nC+1,nO
!       do a=nO+1,nOrb-nR

!         eps = eGF3(p) + e0(a) - e0(i) - e0(j)

!         Bpp(p,1) = Bpp(p,1) & 
!                  + (2d0*V(p,a,i,j) - V(p,a,j,i))*V(p,a,i,j)/eps

!       end do
!     end do
!     end do
!   end do

!   do p=nC+1,nOrb-nR
!     do i=nC+1,nO
!       do a=nO+1,nOrb-nR
!       do b=nO+1,nOrb-nR

!         eps = eGF3(p) + e0(i) - e0(a) - e0(b)

!         Bpp(p,2) = Bpp(p,2) &
!                  + (2d0*V(p,i,a,b) - V(p,i,b,a))*V(p,i,a,b)/eps

!       end do
!       end do
!     end do
!   end do

!   ! Total second-order Green function

!   Sig2(:) = Bpp(:,1) + Bpp(:,2)

!   ! Frequency-dependent third-order contribution: "C" terms

!   Cpp(:,:) = 0d0

!   do p=nC+1,nOrb-nR
!     do i=nC+1,nO
!       do a=nO+1,nOrb-nR
!       do b=nO+1,nOrb-nR
!       do c=nO+1,nOrb-nR
!       do d=nO+1,nOrb-nR

!           eps1 = eGF3(p) + e0(i) - e0(a) - e0(b)
!           eps2 = eGF3(p) + e0(i) - e0(c) - e0(d)

!           Cpp(p,1) = Cpp(p,1) & 
!                    + (2d0*V(p,i,a,b) - V(p,i,b,a))*V(a,b,c,d)*V(p,i,c,d)/(eps1*eps2)

!       end do
!       end do
!       end do
!       end do
!     end do
!   end do

!   do p=nC+1,nOrb-nR
!     do i=nC+1,nO
!     do j=nC+1,nO
!     do k=nC+1,nO
!       do a=nO+1,nOrb-nR
!       do b=nO+1,nOrb-nR

!           eps1 = eGF3(p) + e0(i) - e0(a) - e0(b)
!           eps2 = e0(j)   + e0(k) - e0(a) - e0(b)

!           Cpp(p,2) = Cpp(p,2) & 
!                    + (2d0*V(p,i,a,b) - V(p,i,b,a))*V(a,b,j,k)*V(p,i,j,k)/(eps1*eps2)

!       end do
!       end do
!     end do
!     end do
!     end do
!   end do

!   Cpp(:,3) = Cpp(:,2)

!   do p=nC+1,nOrb-nR
!     do i=nC+1,nO
!     do j=nC+1,nO
!       do a=nO+1,nOrb-nR
!       do b=nO+1,nOrb-nR
!       do c=nO+1,nOrb-nR

!           eps1 = eGF3(p) + e0(a) - e0(i) - e0(j)
!           eps2 = e0(i)   + e0(j) - e0(b) - e0(c)

!           Cpp(p,4) = Cpp(p,4) & 
!                    + (2d0*V(p,a,i,j) - V(p,a,j,i))*V(i,j,b,c)*V(p,a,b,c)/(eps1*eps2)
!       end do
!       end do
!       end do
!     end do
!     end do
!   end do

!   Cpp(:,5) = Cpp(:,4)

!   do p=nC+1,nOrb-nR
!     do i=nC+1,nO
!     do j=nC+1,nO
!     do k=nC+1,nO
!     do l=nC+1,nO
!       do a=nO+1,nOrb-nR

!           eps1 = eGF3(p) + e0(a) - e0(i) - e0(j)
!           eps2 = eGF3(p) + e0(a) - e0(k) - e0(l)

!           Cpp(p,6) = Cpp(p,6) & 
!                    - (2d0*V(p,a,k,l) - V(p,a,l,k))*V(k,l,i,j)*V(p,a,i,j)/(eps1*eps2)
!       end do
!     end do
!     end do
!     end do
!     end do
!   end do

!   ! Frequency-dependent third-order contribution: "D" terms

!   Dpp(:,:) = 0d0

!   do p=nC+1,nOrb-nR
!     do i=nC+1,nO
!     do j=nC+1,nO
!       do a=nO+1,nOrb-nR
!       do b=nO+1,nOrb-nR
!       do c=nO+1,nOrb-nR

!           eps1 = eGF3(p) + e0(i) - e0(a) - e0(b)
!           eps2 = eGF3(p) + e0(j) - e0(b) - e0(c)

!           Dpp(p,1) = Dpp(p,1) &
!                    + V(p,i,a,b)*(V(a,j,i,c)*(    V(p,j,c,b) - 2d0*V(p,j,b,c)) &
!                                + V(a,j,c,i)*(    V(p,j,b,c) - 2d0*V(p,j,c,b)))/(eps1*eps2)

!           Dpp(p,1) = Dpp(p,1) &
!                    + V(p,i,b,a)*(V(a,j,i,c)*(4d0*V(p,j,b,c) - 2d0*V(p,j,c,b)) &
!                                + V(a,j,c,i)*(    V(p,j,c,b) - 2d0*V(p,j,b,c)))/(eps1*eps2)

!       end do
!       end do
!       end do
!     end do
!     end do
!   end do

!   do p=nC+1,nOrb-nR
!     do i=nC+1,nO
!     do j=nC+1,nO
!       do a=nO+1,nOrb-nR
!       do b=nO+1,nOrb-nR
!       do c=nO+1,nOrb-nR

!           eps1 = eGF3(p) + e0(i) - e0(a) - e0(c)
!           eps2 = e0(i)   + e0(j) - e0(a) - e0(b)

!           Dpp(p,2) = Dpp(p,2) &
!                    + V(p,i,c,a)*(V(a,b,i,j)*(4d0*V(p,b,c,j) - 2d0*V(p,b,j,c)) &
!                                + V(a,b,j,i)*(    V(p,b,j,c) - 2d0*V(p,b,c,j)))/(eps1*eps2)

!           Dpp(p,2) = Dpp(p,2) &
!                    + V(p,i,a,c)*(V(a,b,i,j)*(    V(p,b,j,c) - 2d0*V(p,b,c,j)) &
!                                + V(a,b,j,i)*(    V(p,b,c,j) - 2d0*V(p,b,j,c)))/(eps1*eps2)

!       end do
!       end do
!       end do
!     end do
!     end do
!   end do

!   Dpp(:,3) = Dpp(:,2)

!   do p=nC+1,nOrb-nR
!     do i=nC+1,nO
!     do j=nC+1,nO
!     do k=nC+1,nO
!       do a=nO+1,nOrb-nR
!       do b=nO+1,nOrb-nR

!           eps1 = eGF3(p) + e0(a) - e0(j) - e0(k)
!           eps2 = e0(i)   + e0(j) - e0(a) - e0(b)

!           Dpp(p,4) = Dpp(p,4) &
!                    + V(p,a,k,j)*(V(j,i,a,b)*(4d0*V(p,i,k,b) - 2d0*V(p,i,b,k)) &
!                                + V(j,i,b,a)*(    V(p,i,b,k) - 2d0*V(p,i,k,b)))/(eps1*eps2)

!           Dpp(p,4) = Dpp(p,4) &
!                    + V(p,a,j,k)*(V(j,i,a,b)*(    V(p,i,b,k) - 2d0*V(p,i,k,b)) &
!                                + V(j,i,b,a)*(    V(p,i,k,b) - 2d0*V(p,i,b,k)))/(eps1*eps2)

!       end do
!       end do
!     end do
!     end do
!     end do
!   end do

!   Dpp(:,5) = Dpp(:,4)

!   do p=nC+1,nOrb-nR
!     do i=nC+1,nO
!     do j=nC+1,nO
!     do k=nC+1,nO
!       do a=nO+1,nOrb-nR
!       do b=nO+1,nOrb-nR

!           eps1 = eGF3(p) + e0(a) - e0(i) - e0(k)
!           eps2 = eGF3(p) + e0(b) - e0(j) - e0(k)

!           Dpp(p,6) = Dpp(p,6) &
!                    - V(p,a,k,i)*(V(i,b,a,j)*(4d0*V(p,b,k,j) - 2d0*V(p,b,j,k)) &
!                                + V(i,b,j,a)*(    V(p,b,j,k) - 2d0*V(p,b,k,j)))/(eps1*eps2)

!           Dpp(p,6) = Dpp(p,6) &
!                    - V(p,a,i,k)*(V(i,b,a,j)*(    V(p,b,j,k) - 2d0*V(p,b,k,j)) &
!                                + V(i,b,j,a)*(    V(p,b,k,j) - 2d0*V(p,b,j,k)))/(eps1*eps2)

!       end do
!       end do
!     end do
!     end do
!     end do
!   end do

!   ! Compute renormalization factor (if required)

!   Z(:) = 1d0

!   if(renormalization == 0) then

!     Sig3(:) = SigInf(:) &
!             + Cpp(:,1) + Cpp(:,2) + Cpp(:,3) + Cpp(:,4) + Cpp(:,5) + Cpp(:,6) &
!             + Dpp(:,1) + Dpp(:,2) + Dpp(:,3) + Dpp(:,4) + Dpp(:,5) + Dpp(:,6)

!   elseif(renormalization == 1) then

!     Sig3(:) = SigInf(:) &
!             + Cpp(:,1) + Cpp(:,2) + Cpp(:,3) + Cpp(:,4) + Cpp(:,5) + Cpp(:,6) &
!             + Dpp(:,1) + Dpp(:,2) + Dpp(:,3) + Dpp(:,4) + Dpp(:,5) + Dpp(:,6)

!     Z(:) = Cpp(:,2) + Cpp(:,3) + Cpp(:,4) + Cpp(:,5) & 
!          + Dpp(:,2) + Dpp(:,3) + Dpp(:,4) + Dpp(:,5)

!     Z(nC+1:nOrb-nR) = Z(nC+1:nOrb-nR)/Sig2(nC+1:nOrb-nR)
!     Z(:) = 1d0/(1d0 - Z(:))
  
!     Sig3(:) = Z(:)*Sig3(:)

!   elseif(renormalization == 2) then

!     Sig2h1p(:) = Cpp(:,4) + Cpp(:,5) + Cpp(:,6) + Dpp(:,4) + Dpp(:,5) + Dpp(:,6)
!     Sig1h2p(:) = Cpp(:,1) + Cpp(:,2) + Cpp(:,3) + Dpp(:,1) + Dpp(:,2) + Dpp(:,3)

!     X2h1p(:) = Cpp(:,4) + Cpp(:,5) + Dpp(:,4) + Dpp(:,5)
!     X1h2p(:) = Cpp(:,2) + Cpp(:,3) + Dpp(:,2) + Dpp(:,3)

!     X2h1p(nC+1:nOrb-nR) = X2h1p(nC+1:nOrb-nR)/Bpp(nC+1:nOrb-nR,1)
!     X1h2p(nC+1:nOrb-nR) = X1h2p(nC+1:nOrb-nR)/Bpp(nC+1:nOrb-nR,2)

!     Sig3(:) = SigInf(:) +                     &
!             + 1d0/(1d0 - X2h1p(:))*Sig2h1p(:) &
!             + 1d0/(1d0 - X1h2p(:))*Sig1h2p(:)

!   elseif(renormalization == 3) then

!     Sig3(:) = SigInf(:) &
!             + Cpp(:,1) + Cpp(:,2) + Cpp(:,3) + Cpp(:,4) + Cpp(:,5) + Cpp(:,6) &
!             + Dpp(:,1) + Dpp(:,2) + Dpp(:,3) + Dpp(:,4) + Dpp(:,5) + Dpp(:,6)

!     Sig2h1p(:) = Cpp(:,4) + Cpp(:,5) + Cpp(:,6) + Dpp(:,4) + Dpp(:,5) + Dpp(:,6)
!     Sig1h2p(:) = Cpp(:,1) + Cpp(:,2) + Cpp(:,3) + Dpp(:,1) + Dpp(:,2) + Dpp(:,3)

!     X2h1p(:) = Cpp(:,4) + Cpp(:,5) + Dpp(:,4) + Dpp(:,5)
!     X1h2p(:) = Cpp(:,2) + Cpp(:,3) + Dpp(:,2) + Dpp(:,3)

!     X2h1p(nC+1:nOrb-nR) = X2h1p(nC+1:nOrb-nR)/Bpp(nC+1:nOrb-nR,1)
!     X1h2p(nC+1:nOrb-nR) = X1h2p(nC+1:nOrb-nR)/Bpp(nC+1:nOrb-nR,2)

!     Z(:) = X2h1p(:)*Sig2h1p(:) + X1h2p(:)*Sig1h2p(:)
!     Z(nC+1:nOrb-nR) = Z(nC+1:nOrb-nR)/(Sig3(nC+1:nOrb-nR) - SigInf(nC+1:nOrb-nR))
!     Z(:) = 1d0/(1d0 - Z(:))

!     Sig3(:) = Z(:)*Sig3(:)

!   end if

!   ! Total third-order Green function

!    eGF3(:) = e0(:) + Sig2(:) + Sig3(:)

!   ! Print results

!   call print_G0F3(nOrb,nO,e0,Z,eGF3)

! end subroutine RG0F3

subroutine RG0F3(dotest,linearize,eta,doSRG,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,eHF)

! Perform a one-shot third-order Green function calculation

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dotest

  logical,intent(in)            :: linearize
  double precision,intent(in)   :: eta
  logical,intent(in)            :: doSRG

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: ERHF
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: eHF(nOrb)

! Local variables

  double precision              :: Ec
  double precision              :: flow
  double precision,allocatable  :: eGFlin(:)
  double precision,allocatable  :: eGF(:)
  double precision,allocatable  :: SigC(:)
  double precision,allocatable  :: Z(:)

! Hello world

  write(*,*)
  write(*,*)'*******************************'
  write(*,*)'* Restricted G0F3 Calculation *'
  write(*,*)'*******************************'
  write(*,*)

! SRG regularization

  flow = 500d0

  if(doSRG) then

    write(*,*) '*** SRG regularized G0F3 scheme ***'
    write(*,*)

  end if
  
! Memory allocation

  allocate(SigC(nOrb), Z(nOrb), eGFlin(nOrb), eGF(nOrb))
  SigC(:) = 0d0
  Z(:) = 0d0
  eGFlin(:) = 0d0
  eGF(:) = 0d0

! Frequency-dependent third-order contribution

  if(doSRG) then 

    ! TODO call RGF3_SRG_self_energy_diag(flow,nOrb,nC,nO,nV,nR,eHF,ERI,Ec,SigC,Z)

  else

    ! TODO call RGF3_self_energy_diag(eta,nOrb,nC,nO,nV,nR,eHF,ERI,Ec,SigC,Z)

  end if
  
  eGFlin(:) = eHF(:) + Z(:)*SigC(:)

  if(linearize) then

    write(*,*) '*** Quasiparticle energies obtained by linearization ***'

    eGF(:) = eGFlin(:)

  else

    write(*,*) ' *** Quasiparticle energies obtained by root search *** '
    write(*,*)

    ! call RGF3_QP_graph(doSRG,eta,flow,nOrb,nC,nO,nV,nR,eHF,ERI,eGFlin,eHF,eGF,Z)

  end if

  ! Print results

  ! TODO call print_RG0F3(nOrb,nC,nO,nV,nR,eHF,SigC,eGF,Z,ENuc,ERHF,Ec)

! Testing zone

  if(dotest) then

    call dump_test_value('R','G0F2 correlation energy',Ec)
    call dump_test_value('R','G0F2 HOMO energy',eGF(nO))
    call dump_test_value('R','G0F2 LUMO energy',eGF(nO+1))

  end if

  deallocate(SigC, Z, eGFlin, eGF)
  
end subroutine RG0F3
