subroutine huckel_guess(nBas,nO,S,Hc,ERI,J,K,X,cp,F,Fp,e,c,P)

!  Hickel guess of the molecular orbitals for HF calculation

  implicit none

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nO
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: Hc(nBas,nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(inout):: J(nBas,nBas)
  double precision,intent(inout):: K(nBas,nBas)
  double precision,intent(in)   :: X(nBas,nBas)
  double precision,intent(inout):: cp(nBas,nBas)
  double precision,intent(inout):: F(nBas,nBas)
  double precision,intent(inout):: Fp(nBas,nBas)
  double precision,intent(inout):: e(nBas)
  double precision,intent(inout):: P(nBas,nBas)

! Local variables

  integer                       :: mu,nu
  double precision              :: a

! Output variables

  double precision,intent(out)  :: c(nBas,nBas)

  a = 1.75d0

  Fp(:,:) = matmul(transpose(X(:,:)),matmul(Hc(:,:),X(:,:)))
  cp(:,:) = Fp(:,:)
  call diagonalize_matrix(nBas,cp,e)
  c(:,:) = matmul(X(:,:),cp(:,:))

  call Coulomb_matrix_AO_basis(nBas,P,ERI,J)
  call exchange_matrix_AO_basis(nBas,P,ERI,K)

  F(:,:) = Hc(:,:) + J(:,:) + 0.5d0*K(:,:)

  do mu=1,nBas
    do nu=mu+1,nBas

      F(mu,nu) = 0.5d0*a*S(mu,nu)*(Hc(mu,mu) + Hc(nu,nu))
      F(nu,mu) = F(mu,nu)

    enddo
  enddo
  
  Fp(:,:) = matmul(transpose(X(:,:)),matmul(F(:,:),X(:,:)))
  cp(:,:) = Fp(:,:)
  call diagonalize_matrix(nBas,cp,e)
  c(:,:) = matmul(X(:,:),cp(:,:))

end subroutine
