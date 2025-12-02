subroutine orbital_gradient_numerically(O,V,N,nS,Nsq,Hc,c,ERI_AO,delta,grad)
      
! Compute the orbital gradient given numerically for RPA energy

  implicit none
  include 'parameters.h'
      
! Input variables

  integer,intent(in)            :: O,V,N,Nsq,nS
  double precision,intent(in)   :: Hc(N,N),c(N,N)
  double precision,intent(in)   :: ERI_AO(N,N,N,N)
  double precision,intent(in)   :: delta

! Local variables
  double precision              :: Eplus,Eminus
  double precision,allocatable  :: Kap(:,:)
  double precision,allocatable  :: c_local(:,:)
  integer                       :: p,q,r,s,pq,qp

! Output variables

  double precision,intent(out)  :: grad(Nsq)

  allocate(Kap(N,N),c_local(N,N))
  
  Kap(:,:)  = 0d0
  grad(:)   = 0d0
  
  pq = 0
  do p=1,N-1
    do q=p+1,N
      
      ! ERPA(+delta)
      c_local(:,:)  = c(:,:)
      Kap      = 0d0
      Kap(p,q) = delta
      Kap(q,p) = - delta
      call rotate_orbitals(N,Kap,c_local)
      call get_rpa_energy(O,V,N,nS,Hc,c_local,ERI_AO,Eplus)
  
      ! ERPA(-delta)
      c_local(:,:)  = c(:,:)
      Kap      = 0d0
      Kap(p,q) = - delta
      Kap(q,p) = delta
      call rotate_orbitals(N,Kap,c_local)
      call get_rpa_energy(O,V,N,nS,Hc,c_local,ERI_AO,Eminus)
      pq = q + (p-1)*N
      qp = p + (q-1)*N
      grad(pq) = (Eplus - Eminus)/(2*delta)
      grad(qp) = - grad(pq)
    enddo
  enddo

deallocate(Kap,c_local)
end subroutine
