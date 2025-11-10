subroutine G_optimize_orbitals(nBas,nBas2,nV,nR,nC,nO,N,Nsq,O,V,ERI_AO,ERI_MO,h,F,rdm1_hf,rdm1_c,rdm2_hf,rdm2_c,c,OOConv)

! Calculate gradient and Hessian and rotate orbitals accordingly. Returns ERI and h in new MOs.

  implicit none
  include 'parameters.h'
  include 'quadrature.h'

! Input variables

  integer,intent(in)            :: nBas,nBas2,nV,nR,nC,nO,N,V,O,Nsq
  double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: h(nBas2,nBas2)
  double precision,intent(in)   :: F(nBas2,nBas2)
  double precision,intent(in)   :: rdm1_hf(N,N)
  double precision,intent(in)   :: rdm1_c(N,N)
  double precision,intent(in)   :: rdm2_hf(N,N,N,N)
  double precision,intent(in)   :: rdm2_c(N,N,N,N)

! Local variables
  
  integer                       :: p,q,pq,r,s,rs
  double precision,allocatable  :: hess(:,:), grad(:), grad_tmp(:), hessInv(:,:), hess_tmp(:,:)
  double precision,allocatable  :: Kap(:,:), ExpKap(:,:)

! Output variables

  double precision,intent(inout)   :: ERI_MO(nBas2,nBas2,nBas2,nBas2)
  double precision,intent(inout)   :: c(nBas2,nBas2)
  double precision,intent(out)     :: OOConv

  !--------------------------!
  ! Compute orbital gradient !
  !--------------------------!
 
  allocate(grad(Nsq),grad_tmp(Nsq))
  call orbital_gradient(O,V,N,Nsq,h,ERI_MO,rdm1_hf,rdm2_hf,grad_tmp)
  grad = grad_tmp
  call orbital_gradient(O,V,N,Nsq,F,ERI_MO,rdm1_c,rdm2_c,grad_tmp)
  grad  = grad + grad_tmp
  deallocate(grad_tmp)

  ! Check convergence of orbital optimization
  OOConv = maxval(grad)
  
  !-------------------------!
  ! Compute orbital Hessian !
  !-------------------------!
 
  allocate(hess(Nsq,Nsq),hess_tmp(Nsq,Nsq))
  hess_tmp(:,:) = 0d0 
  call orbital_hessian(O,V,N,Nsq,h,ERI_MO,rdm1_hf,rdm2_hf,hess_tmp)
  hess = hess_tmp
  call orbital_hessian(O,V,N,Nsq,F,ERI_MO,rdm1_c,rdm2_c,hess_tmp)
  hess = hess_tmp + hess
  
  deallocate(hess_tmp)
  allocate(hessInv(Nsq,Nsq))
 
  call pseudo_inverse_matrix(Nsq,hess,hessInv)
  
  deallocate(hess)
 
  allocate(Kap(N,N))
 
  Kap(:,:) = 0d0
 
  pq = 0
  do p=1,nBas2
    do q=1,nBas2
 
      pq = pq + 1
 
      rs = 0
      do r=1,nBas2
        do s=1,nBas2
 
          rs = rs + 1
 
            Kap(p,q) = Kap(p,q) - hessInv(pq,rs)*grad(rs)
 
        end do
      end do
 
    end do
  end do
 
  deallocate(hessInv,grad)
 
  allocate(ExpKap(N,N))
  call matrix_exp(N,Kap,ExpKap)
  deallocate(Kap)
 
  c = matmul(c,ExpKap)
  deallocate(ExpKap)

end subroutine
