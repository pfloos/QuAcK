subroutine R_optimize_orbitals(diagHess,dRPA,nBas,nOrb,nV,nR,nC,nO,N,Nsq,O,V,&
                ERI_AO,ERI_AO_AS,ERI_MO,ERI_MO_AS,                           &
                Hc,h,F,rdm1_hf,rdm1_c,rdm2_hf,rdm2_c,c,OOConv)

! Calculate gradient and Hessian and rotate orbitals accordingly. Returns ERI and h in new MOs.

  implicit none
  include 'parameters.h'
  include 'quadrature.h'

! Input variables

  integer,intent(in)               :: nBas,nOrb,nV,nR,nC,nO,N,V,O,Nsq
  double precision,intent(in)      :: ERI_AO(nBas,nBas,nBas,nBas)
  double precision,intent(in)      :: ERI_AO_AS(nBas,nBas,nBas,nBas)
  double precision,intent(in)      :: Hc(nBas,nBas)
  double precision,intent(in)      :: F(nOrb,nOrb)
  double precision,intent(in)      :: rdm1_hf(N,N)
  double precision,intent(in)      :: rdm1_c(N,N)
  double precision,intent(in)      :: rdm2_hf(N,N,N,N)
  double precision,intent(in)      :: rdm2_c(N,N,N,N)
  logical,intent(in)               :: diagHess
  logical,intent(in)               :: dRPA

! Local variables
  
  integer                          :: p,q,pq,r,s,rs
  integer                          :: nhess,mhess
  double precision,allocatable     :: hess(:,:), grad(:),grad_tmp(:),hess_tmp(:,:), hessInv(:,:)
  double precision,allocatable     :: Kap(:,:), ExpKap(:,:)

! Output variables

  double precision,intent(inout)   :: ERI_MO(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(inout)   :: ERI_MO_AS(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(inout)   :: c(nBas,nOrb)
  double precision,intent(out)     :: OOConv
  double precision,intent(in)      :: h(nOrb,nOrb)

  !--------------------------!
  ! Compute orbital gradient !
  !--------------------------!
 
  allocate(grad(Nsq),grad_tmp(Nsq))
  grad(:) = 0d0
  nhess = Nsq
  mhess = Nsq
  call orbital_gradient(O,V,N,Nsq,h,ERI_MO,rdm1_hf,rdm2_hf,grad_tmp)
  grad = grad + grad_tmp
  if(dRPA) then
    call orbital_gradient(O,V,N,Nsq,F,ERI_MO,rdm1_c,rdm2_c,grad_tmp)
  else
    call orbital_gradient(O,V,N,Nsq,F,ERI_MO_AS,rdm1_c,rdm2_c,grad_tmp)
  endif
  grad = grad + grad_tmp 
  deallocate(grad_tmp)

  ! Check convergence of orbital optimization
  OOConv = maxval(grad)
  
  !-------------------------!
  ! Compute orbital Hessian !
  !-------------------------!

  if(diagHess) then
    mhess = 1
    allocate(hess(nhess,mhess),hess_tmp(nhess,mhess))
    hess(:,:) = 0d0
    call orbital_hessian_diag(O,V,N,Nsq,h,ERI_MO,rdm1_hf,rdm2_hf,hess_tmp)
    hess = hess + hess_tmp
    if(dRPA) then
      call orbital_hessian_diag(O,V,N,Nsq,F,ERI_MO,rdm1_c,rdm2_c,hess_tmp)
    else
      call orbital_hessian_diag(O,V,N,Nsq,F,ERI_MO_AS,rdm1_c,rdm2_c,hess_tmp)
    endif
    hess = hess + hess_tmp
    deallocate(hess_tmp)
    allocate(hessInv(nhess,mhess))
    do pq=1,Nsq
      if(abs(hess(pq,1))>1e-15) then
        hessInv(pq,1) = 1/hess(pq,1)
      endif
    enddo
  else
    allocate(hess(nhess,mhess),hess_tmp(nhess,mhess))
    hess(:,:) = 0d0
    call orbital_hessian(O,V,N,Nsq,h,ERI_MO,rdm1_hf,rdm2_hf,hess_tmp)
    hess = hess + hess_tmp 
    if(dRPA) then
      call orbital_hessian(O,V,N,Nsq,F,ERI_MO,rdm1_c,rdm2_c,hess_tmp)
    else
      call orbital_hessian(O,V,N,Nsq,F,ERI_MO_AS,rdm1_c,rdm2_c,hess_tmp)
    endif
    hess = hess + hess_tmp 
    deallocate(hess_tmp)
    allocate(hessInv(nhess,mhess))
    call pseudo_inverse_matrix(Nsq,hess,hessInv)
    !call inverse_matrix(Nsq,hess,hessInv)
  endif
  
  deallocate(hess)
 
  allocate(Kap(N,N))
 
  Kap(:,:) = 0d0
 
  pq = 0
  do p=1,nOrb
    do q=1,nOrb
 
      pq = pq + 1
 
      rs = 0
      if(diagHess) then
        Kap(p,q) = Kap(p,q) - hessInv(pq,1)*grad(pq)
      else
        do r=1,nOrb
          do s=1,nOrb 
            rs = rs + 1
            Kap(p,q) = Kap(p,q) - hessInv(pq,rs)*grad(rs) 
          end do
        end do
      endif
    end do
  end do
 
  deallocate(hessInv,grad)
  
!  write(*,*) 'kappa'
!  call matout(nOrb,nOrb,Kap)
!  write(*,*)
 
  allocate(ExpKap(N,N))
  call matrix_exp(N,Kap,ExpKap)
  deallocate(Kap)
 
!  write(*,*) 'e^kappa'
!  call matout(N,N,ExpKap)
!  write(*,*)
 
!  write(*,*) 'Old orbitals'
!  call matout(nBas,nOrb,c)
!  write(*,*)
 
  c = matmul(c,ExpKap)
  deallocate(ExpKap)
 
!  write(*,*) 'Rotated orbitals'
!  call matout(nBas,nOrb,c)
!  write(*,*)

  ! Transform integrals and Hc
  call AOtoMO(nBas,nOrb,c,Hc,h)
  call AOtoMO_ERI_RHF(nBas,nOrb,c,ERI_AO,ERI_MO)
  if(.not. dRPA) call AOtoMO_ERI_RHF(nBas,nOrb,c,ERI_AO_AS,ERI_MO_AS)

end subroutine
