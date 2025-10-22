subroutine R_optimize_orbitals(diagHess,nBas,nOrb,nV,nR,nC,nO,N,Nsq,O,V,ERI_AO,ERI_MO,Hc,h,rdm1,rdm2,c,OOConv)

! Calculate gradient and Hessian and rotate orbitals accordingly. Returns ERI and h in new MOs.

  implicit none
  include 'parameters.h'
  include 'quadrature.h'

! Input variables

  integer,intent(in)               :: nBas,nOrb,nV,nR,nC,nO,N,V,O,Nsq
  double precision,intent(in)      :: ERI_AO(nBas,nBas,nBas,nBas)
  double precision,intent(in)      :: Hc(nBas,nBas)
  double precision,intent(in)      :: rdm1(N,N)
  double precision,intent(in)      :: rdm2(N,N,N,N)
  logical,intent(in)               :: diagHess

! Local variables
  
  integer                          :: p,q,pq,r,s,rs
  integer                          :: nhess,mhess
  double precision,allocatable     :: hess(:,:), grad(:), hessInv(:,:)
  double precision,allocatable     :: Kap(:,:), ExpKap(:,:)
  double precision                 :: reg = 1e-14 ! regularisation of Hessian H <- H + reg*I

! Output variables

  double precision,intent(inout)   :: ERI_MO(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(inout)   :: c(nBas,nOrb)
  double precision,intent(out)     :: OOConv
  double precision,intent(in)      :: h(nOrb,nOrb)

  !--------------------------!
  ! Compute orbital gradient !
  !--------------------------!
 
  allocate(grad(Nsq))
  grad(:) = 0d0
  nhess = Nsq
  mhess = Nsq
  call orbital_gradient(O,V,N,Nsq,h,ERI_MO,rdm1,rdm2,grad)
 
  ! Check convergence of orbital optimization
  OOConv = maxval(grad)
  
  !-------------------------!
  ! Compute orbital Hessian !
  !-------------------------!

  if(diagHess) then
    mhess = 1
    allocate(hess(nhess,mhess),hessInv(nhess,mhess))
    hess(:,:) = 0d0
    call orbital_hessian_diag(O,V,N,Nsq,h,ERI_MO,rdm1,rdm2,hess)
    do pq=1,Nsq
      hess(pq,1) = hess(pq,1) + reg
      hessInv(pq,1) = 1/hess(pq,1)
    enddo
  else
    allocate(hess(nhess,mhess),hessInv(nhess,mhess))
    hess(:,:) = 0d0 
    call orbital_hessian(O,V,N,Nsq,h,ERI_MO,rdm1,rdm2,hess)
    call pseudo_inverse_matrix(Nsq,hess,hessInv)
  endif
  
!  write(*,*) "Hessian"
!  call matout(nhess,mhess,hess)
  deallocate(hess)
!  write(*,*) "Inv Hessian"
!  call matout(nhess,mhess,hessInv)
 
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
  call matrix_exponential(N,Kap,ExpKap)
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

end subroutine
