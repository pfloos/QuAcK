subroutine eNcusp(nEl,nBas,S,T,V,G,X,ENuc,EHF,c,e,P,F)

! Perform restricted Hartree-Fock calculation

  implicit none

! Input variables

  integer,intent(in)            :: nEl,nBas
  double precision,intent(in)   :: ENuc,EHF
  double precision,intent(in)   :: S(nBas,nBas),T(nBas,nBas),V(nBas,nBas),G(nBas,nBas,nBas,nBas),X(nBas,nBas)
  double precision,intent(out)  :: c(nBas,nBas),e(nBas),P(nBas,nBas),F(nBas,nBas)

! Local variables

  integer,parameter             :: maxSCF = 128
  double precision,parameter    :: thresh = 1d-6
  integer                       :: nO,nSCF,lwork,info
  double precision              :: ET,EV,Conv,Gap
  double precision,allocatable  :: Hc(:,:),cp(:,:),cO(:,:),Fp(:,:),work(:)

  integer                       :: mu,nu,lambda,sigma,i

! Output variables

! Number of occupied orbitals
  if(mod(nEl,2) /= 0) then
    write(*,*) 'closed-shell system required!'
    stop
  endif
  nO = nEl/2

! Memory allocation
  allocate(Hc(nBas,nBas),cp(nBas,nBas),cO(nBas,nO),Fp(nBas,nBas))
  lwork = 3*nBas
  allocate(work(lwork))

! Core Hamiltonian
  Hc = T + V


end subroutine eNcusp
