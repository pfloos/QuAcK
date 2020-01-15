subroutine huckel_guess(nBas,nO,S,Hc,ERI,J,K,X,cp,Fp,e,c,P)

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
  double precision,intent(inout):: Fp(nBas,nBas)
  double precision,intent(inout):: e(nBas)
  double precision,intent(inout):: P(nBas,nBas)

! Local variables

  integer                       :: mu,nu
  double precision              :: a

! Output variables

  double precision,intent(out)  :: c(nBas,nBas)

  a = 1.75d0

  Fp(:,:) = Hc(:,:)

  do mu=1,nBas
    do nu=mu+1,nBas

      Fp(mu,nu) = 0.5d0*a*S(mu,nu)*(Hc(mu,mu) + Hc(nu,nu))
      Fp(nu,mu) = Fp(mu,nu)

    enddo
  enddo
  
  Fp(:,:) = matmul(transpose(X(:,:)),matmul(Fp(:,:),X(:,:)))
  cp(:,:) = Fp(:,:)
  call diagonalize_matrix(nBas,cp,e)
  c(:,:) = matmul(X(:,:),cp(:,:))
  P(:,:) = 2d0*matmul(c(:,1:nO),transpose(c(:,1:nO)))

! call Coulomb_matrix_AO_basis(nBas,P,ERI,J)
! call exchange_matrix_AO_basis(nBas,P,ERI,K)

! do mu=1,nBas

!   Fp(mu,mu) = Hc(mu,mu) + J(mu,mu) + 0.5d0*K(mu,mu)

!   do nu=mu+1,nBas

!     Fp(mu,nu) = 0.5d0*a*S(mu,nu)*(Hc(mu,mu) + Hc(nu,nu))
!     Fp(nu,mu) = Fp(mu,nu)

!   enddo

! enddo

!   Fp = matmul(transpose(X),matmul(Fp,X))

!   cp(:,:) = Fp(:,:)
!   call diagonalize_matrix(nBas,cp,e)
!   c = matmul(X,cp)

end subroutine
