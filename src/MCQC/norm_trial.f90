subroutine norm_trial(nBas,nO,c,P,Norm,NormSq)

! Initialize weight function

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas,nO
  double precision,intent(inout):: c(nBas,nO),P(nBas,nBas)

! Local variables

  double precision,allocatable  :: S(:,:),T(:,:),V(:,:),Hc(:,:),G(:,:,:,:)

  integer                       :: mu,nu,la,si

! Output variables

  double precision,intent(inout):: Norm,NormSq

! Memory allocation for one- and two-electron integrals

  allocate(S(nBas,nBas),T(nBas,nBas),V(nBas,nBas),Hc(nBas,nBas),G(nBas,nBas,nBas,nBas))

! Read integrals

  call read_integrals(nBas,S,T,V,Hc,G)

! Compute normalization factor

  P = 2d0*matmul(c,transpose(c))

  Norm = 0d0
  do mu=1,nBas
    do nu=1,nBas
      do la=1,nBas
        do si=1,nBas
          Norm = Norm + P(mu,nu)*P(la,si)*G(mu,la,nu,si)
        enddo
      enddo
    enddo
  enddo

  Norm = Norm*Norm
  NormSq = Norm*Norm
  
  write(*,*)
  write(*,*) 'Normalization of trial wave function: ',Norm
  write(*,*)

end subroutine norm_trial
