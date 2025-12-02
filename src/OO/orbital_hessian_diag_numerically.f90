subroutine orbital_hessian_diag_numerically(O,V,N,nS,Nsq,Hc,c,ERI_AO,delta,hess)

  implicit none
  include 'parameters.h'

  ! Input
  integer,intent(in)            :: O,V,N,Nsq,nS
  double precision,intent(in)   :: Hc(N,N), c(N,N)
  double precision,intent(in)   :: ERI_AO(N,N,N,N)
  double precision,intent(in)   :: delta

  ! Output
  double precision,intent(out)  :: hess(Nsq,1)

  ! Local
  double precision              :: E0, Eplus, Eminus
  double precision,allocatable  :: Kap(:,:), c_local(:,:)
  integer                       :: r,s,rs

  allocate(Kap(N,N), c_local(N,N))

  hess(:,:) = 0d0

  ! --- Reference energy E(0) ---
  call get_rpa_energy(O,V,N,nS,Hc,c,ERI_AO,E0)

  do r = 1, N-1
    do s = r+1, N

      rs = s + (r-1)*N

      ! ---- E(+delta) ----
      c_local(:,:) = c(:,:)
      Kap(:,:) = 0d0
      Kap(r,s) =  delta
      Kap(s,r) = -delta
      call rotate_orbitals(N,Kap,c_local)
      call get_rpa_energy(O,V,N,nS,Hc,c_local,ERI_AO,Eplus)

      ! ---- E(-delta) ----
      c_local(:,:) = c(:,:)
      Kap(:,:) = 0d0
      Kap(r,s) = -delta
      Kap(s,r) =  delta
      call rotate_orbitals(N,Kap,c_local)
      call get_rpa_energy(O,V,N,nS,Hc,c_local,ERI_AO,Eminus)

      ! ---- Diagonal Hessian entry ----
      hess(rs,1) = (Eplus - 2d0*E0 + Eminus) / (delta*delta)

    enddo
  enddo

  deallocate(Kap, c_local)

end subroutine

