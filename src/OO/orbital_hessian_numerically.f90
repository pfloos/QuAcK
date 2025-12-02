subroutine orbital_hessian_numerically(O,V,N,nS,Nsq,Hc,c,ERI_AO,delta,hess)

  implicit none
  include 'parameters.h'

  integer,intent(in)            :: O,V,N,Nsq,nS
  double precision,intent(in)   :: Hc(N,N), c(N,N)
  double precision,intent(in)   :: ERI_AO(N,N,N,N)
  double precision,intent(in)   :: delta
  double precision,intent(out)  :: hess(Nsq,Nsq)

  double precision,allocatable  :: Kap(:,:), c_local(:,:)
  double precision,allocatable  :: grad_plus(:), grad_minus(:)

  integer                       :: r,s,rs,sr

  allocate(Kap(N,N), c_local(N,N))
  allocate(grad_plus(Nsq), grad_minus(Nsq))

  hess(:,:) = 0d0

  rs = 0
  do r = 1, N-1
    do s = r+1, N

      rs = s + (r-1)*N
      sr = r + (s-1)*N

      ! ----- Gradient(+delta) -----
      c_local(:,:) = c(:,:)
      Kap = 0d0
      Kap(r,s) =  delta
      Kap(s,r) = -delta
      call rotate_orbitals(N,Kap,c_local)

      call orbital_gradient_numerically(O,V,N,nS,Nsq,Hc,c_local,ERI_AO,delta,grad_plus)

      ! ----- Gradient(-delta) -----
      c_local(:,:) = c(:,:)
      Kap = 0d0
      Kap(r,s) = -delta
      Kap(s,r) =  delta
      call rotate_orbitals(N,Kap,c_local)

      call orbital_gradient_numerically(O,V,N,nS,Nsq,Hc,c_local,ERI_AO,delta,grad_minus)

      hess(:,rs) = (grad_plus(:) - grad_minus(:)) / (2d0*delta)

      ! Use antisymmetry of kappa
      hess(:,sr) = -hess(:,rs)

    enddo
  enddo

  deallocate(Kap, c_local, grad_plus, grad_minus)

end subroutine
