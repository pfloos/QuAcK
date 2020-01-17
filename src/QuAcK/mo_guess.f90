subroutine mo_guess(nBas,nO,guess_type,S,Hc,ERI,J,K,X,cp,F,Fp,e,c,P)

!  Guess of the molecular orbitals for HF calculation

  implicit none

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nO
  integer,intent(in)            :: guess_type
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

  integer                       :: nSCF

! Output variables

  double precision,intent(out)  :: c(nBas,nBas)

  if(guess_type == 1) then

    Fp = matmul(transpose(X),matmul(Hc,X))
    cp(:,:) = Fp(:,:)
    call diagonalize_matrix(nBas,cp,e)
    c = matmul(X,cp)

  elseif(guess_type == 2) then

    call huckel_guess(nBas,nO,S,Hc,ERI,J,K,X,cp,F,Fp,e,c,P)

  elseif(guess_type == 3) then

    call random_number(c)

  else

    print*,'Wrong guess option'
    stop

  endif

  P(:,:) = 2d0*matmul(c(:,1:nO),transpose(c(:,1:nO)))

end subroutine
