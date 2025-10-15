subroutine R_optimize_orbitals(nBas,nOrb,nV,nR,nC,nO,N,Nsq,O,V,ERI_AO,ERI_MO,h,rdm1,rdm2,c,OOConv)

! Calculate gradient and Hessian and rotate orbitals accordingly. Returns ERI and h in new MOs.

  implicit none
  include 'parameters.h'
  include 'quadrature.h'

! Input variables

  integer,intent(in)            :: nBas,nOrb,nV,nR,nC,nO,N,V,O,Nsq
  double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: h(nOrb,nOrb)
  double precision,intent(in)   :: rdm1(N,N)
  double precision,intent(in)   :: rdm2(N,N,N,N)

! Local variables
  
  integer                       :: p,q,pq,r,s,rs
  double precision,allocatable  :: hess(:,:), grad(:), hessInv(:,:)
  double precision,allocatable  :: Kap(:,:), ExpKap(:,:)

! Output variables

  double precision,intent(inout)   :: ERI_MO(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(inout)   :: c(nBas,nOrb)
  double precision,intent(out)     :: OOConv

  !--------------------------!
  ! Compute orbital gradient !
  !--------------------------!
 
  allocate(grad(Nsq))
  grad(:) = 0d0
  call orbital_gradient(O,V,N,Nsq,h,ERI_MO,rdm1,rdm2,grad)
 
  ! Check convergence of orbital optimization
  OOConv = maxval(grad)
  
  !-------------------------!
  ! Compute orbital Hessian !
  !-------------------------!
 
  allocate(hess(Nsq,Nsq))
  hess(:,:) = 0d0 
  call orbital_hessian(O,V,N,Nsq,h,ERI_MO,rdm1,rdm2,hess)
  
  write(*,*) "Hessian"
  call matout(Nsq,Nsq,hess)
 
  allocate(hessInv(Nsq,Nsq))
 
  call inverse_matrix(Nsq,hess,hessInv)
  
  write(*,*) "Inv Hessian"
  call matout(Nsq,Nsq,hessInv)
  
  deallocate(hess)
 
  allocate(Kap(N,N))
 
  Kap(:,:) = 0d0
 
  pq = 0
  do p=1,nOrb
    do q=1,nOrb
 
      pq = pq + 1
 
      rs = 0
      do r=1,nOrb
        do s=1,nOrb
 
          rs = rs + 1
 
            Kap(p,q) = Kap(p,q) - hessInv(pq,rs)*grad(rs)
 
        end do
      end do
 
    end do
  end do
 
  deallocate(hessInv,grad)

  write(*,*) 'kappa'
  call matout(nOrb,nOrb,Kap)
  write(*,*)
 
  allocate(ExpKap(N,N))
  call matrix_exponential(N,Kap,ExpKap)
  deallocate(Kap)
 
  write(*,*) 'e^kappa'
  call matout(N,N,ExpKap)
  write(*,*)
 
  write(*,*) 'Old orbitals'
  call matout(nBas,nOrb,c)
  write(*,*)
 
  c = matmul(c,ExpKap)
  deallocate(ExpKap)
 
  write(*,*) 'Rotated orbitals'
  call matout(nBas,nOrb,c)
  write(*,*)

end subroutine
