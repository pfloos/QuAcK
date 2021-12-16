subroutine unrestricted_linear_response_pp(ispin,TDA,nBas,nC,nO,nV,nR,nPaa,nPab,nPbb,nPt, &
nHaa,nHab,nHbb,nHt,nS_sc,lambda,e,ERI_aaaa,ERI_aabb,ERI_bbbb,Omega1,X1,Y1,Omega2,X2,Y2,&
EcRPA)

! Compute linear response for unrestricted formalism

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: ispin 
  logical,intent(in)            :: TDA  
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC(nspin)
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nV(nspin)
  integer,intent(in)            :: nR(nspin)
  integer,intent(in)            :: nPaa
  integer,intent(in)            :: nPab
  integer,intent(in)            :: nPbb
  integer,intent(in)            :: nPt
  integer,intent(in)            :: nHaa
  integer,intent(in)            :: nHab
  integer,intent(in)            :: nHbb
  integer,intent(in)            :: nHt
  integer,intent(in)            :: nS_sc
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: e(nBas,nspin)
  double precision,intent(in)   :: ERI_aaaa(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_aabb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_bbbb(nBas,nBas,nBas,nBas)
 
! Local variables

  double precision,external     :: trace_matrix
  double precision              :: EcRPA1
  double precision              :: EcRPA2
  double precision,allocatable  :: C(:,:)
  double precision,allocatable  :: B(:,:)
  double precision,allocatable  :: D(:,:)
  double precision,allocatable  :: M(:,:)
  double precision,allocatable  :: Z(:,:)
  double precision,allocatable  :: Omega(:)

! Output variables

  double precision,intent(out)  :: Omega1(nPt)
  double precision,intent(out)  :: X1(nPt,nPt)
  double precision,intent(out)  :: Y1(nHt,nPt)
  double precision,intent(out)  :: Omega2(nHt)
  double precision,intent(out)  :: X2(nPt,nHt)
  double precision,intent(out)  :: Y2(nHt,nHt)
  double precision,intent(out)  :: EcRPA
  

! Memory allocation

 allocate(C(nPt,nPt),B(nPt,nHt),D(nHt,nHt),M(nPt+nHt,nPt+nHt),Z(nPt+nHt,nPt+nHt),&          
Omega(nPt+nHt))

! Build C, B and D matrices for the pp channel

  call unrestricted_linear_response_C_pp(ispin,nBas,nC,nO,nV,nR,nPaa,nPab,nPbb,nPt,lambda,&
e,ERI_aaaa,ERI_aabb,ERI_bbbb,C)

    call unrestricted_linear_response_B_pp(ispin,nBas,nC,nO,nV,nR,nPaa,nPab,nPbb,nPt,nHaa,&
nHab,nHbb,nHt,lambda,ERI_aaaa,ERI_aabb,ERI_bbbb,B)

call unrestricted_linear_response_D_pp(ispin,nBas,nC,nO,nV,nR,nHaa,nHab,nHbb,nHt,lambda,&
ERI_aaaa,ERI_aabb,ERI_bbbb,D)

! Diagonal blocks 

    M(    1:nPt    ,    1:nPt)     = + C(1:nPt,1:nPt)
    M(nPt+1:nPt+nHt,nPt+1:nPt+nHt) = - D(1:nHt,1:nHt)

  ! Off-diagonal blocks

    M(    1:nPt    ,nPt+1:nHt+nPt) = -           B(1:nPt,1:nHt)
    M(nPt+1:nHt+nPt,    1:nPt)     = + transpose(B(1:nPt,1:nHt))

! Diagonalize the p-h matrix

    if(nHt+nPt > 0) call diagonalize_general_matrix(nHt+nPt,M,Omega,Z)

  ! Split the various quantities in p-p and h-h parts

    call sort_ppRPA(nHt,nPt,Omega(:),Z(:,:),Omega1(:),X1(:,:),Y1(:,:),Omega2(:),X2(:,:),&
Y2(:,:))

  ! end if Pourquoi ne faut-il pas de end if ici ? 

! Compute the RPA correlation energy

  EcRPA = 0.5d0*( sum(Omega1(:)) - sum(Omega2(:)) - trace_matrix(nPt,C(:,:)) - trace_matrix(nHt,D(:,:)) )
  EcRPA1 = +sum(Omega1(:)) - trace_matrix(nPt,C(:,:))
  EcRPA2 = -sum(Omega2(:)) - trace_matrix(nHt,D(:,:))
  if(abs(EcRPA - EcRPA1) > 1d-6 .or. abs(EcRPA - EcRPA2) > 1d-6) &
    print*,'!!! Issue in pp-RPA linear reponse calculation RPA1 != RPA2 !!!'
 
end subroutine unrestricted_linear_response_pp