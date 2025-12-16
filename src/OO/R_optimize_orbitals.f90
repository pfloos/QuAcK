subroutine R_optimize_orbitals(diagHess,OVRotOnly,dRPA,nBas,nOrb,nV,nR,nC,nO,N,Nsq,O,V,&
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
  logical,intent(in)               :: OVRotOnly

! Local variables
  
  integer                          :: p,q,pq,r,s,rs,ia,jb,i,j,a,b
  integer                          :: nhess,mhess
  double precision,allocatable     :: hess(:,:), grad(:),grad_tmp(:),hess_tmp(:,:), hessInv(:,:)
  double precision,allocatable     :: Kap(:,:), ExpKap(:,:)
  double precision,allocatable     :: rdm1(:,:), rdm2(:,:,:,:)
  double precision,allocatable     :: VL(:,:), VR(:,:)
  double precision,allocatable     :: WR(:)
  double precision                 :: tol
  double precision                 :: tstart,tend,tdiff

! Output variables

  double precision,intent(inout)   :: ERI_MO(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(inout)   :: ERI_MO_AS(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(inout)   :: c(nBas,nOrb)
  double precision,intent(out)     :: OOConv
  double precision,intent(in)      :: h(nOrb,nOrb)

  !--------------------------!
  ! Compute orbital gradient !
  !--------------------------!
 
  allocate(grad(Nsq),grad_tmp(Nsq),rdm1(N,N),rdm2(N,N,N,N))
  rdm1(:,:) = 0d0
  rdm2(:,:,:,:) = 0d0
  grad(:) = 0d0
  nhess = Nsq
  mhess = Nsq
  rdm1 = rdm1_hf + rdm1_c
  rdm2 = rdm2_hf + rdm2_c
  call orbital_gradient(O,V,N,Nsq,h,ERI_MO,rdm1_hf,rdm2_hf,grad_tmp)
  grad = grad + grad_tmp
  write(*,*) "grad_hf from rdm"
  call matout(N,N,grad)
  if(dRPA) then
    call orbital_gradient(O,V,N,Nsq,h,ERI_MO,rdm1_c,rdm2_c,grad_tmp)
  else
    call orbital_gradient(O,V,N,Nsq,F,ERI_MO_AS,rdm1_c,rdm2_c,grad_tmp)
  endif
  grad = grad + grad_tmp
  write(*,*) "rdm1_crpa"
  call matout(N,N,rdm1_c)
  write(*,*) "rdm2_crpa"
  call matout(Nsq,Nsq,rdm2_c)
  write(*,*) "rdm2_hf"
  call matout(Nsq,Nsq,rdm2_hf)
  write(*,*) "grad_hf + grad_crpa from rdm"
  call matout(N,N,grad)
  write(*,*) "grad from rdm"
  grad = 0d0
  call orbital_gradient(O,V,N,Nsq,h,ERI_MO,rdm1,rdm2,grad)
  call matout(N,N,grad)

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
   ! if(dRPA) then
   !   call orbital_hessian_diag(O,V,N,Nsq,h,ERI_MO,rdm1_c,rdm2_c,hess_tmp)
   ! else
   !   call orbital_hessian_diag(O,V,N,Nsq,F,ERI_MO_AS,rdm1_c,rdm2_c,hess_tmp)
   ! endif
   ! hess = hess + hess_tmp
    deallocate(hess_tmp)
    allocate(hessInv(nhess,mhess))
    tol = 1d-12 * maxval(abs(hess(:,1)))
    do pq=1,Nsq
      if(abs(hess(pq,1))>tol) then
        hessInv(pq,1) = 1/hess(pq,1)
      endif
    enddo
  elseif(.not. OVRotOnly) then
    allocate(hess(nhess,mhess),hess_tmp(nhess,mhess))
    hess(:,:) = 0d0
    !call orbital_hessian(O,V,N,Nsq,h,ERI_MO,rdm1_hf,rdm2_hf,hess_tmp)
    call orbital_hessian(O,V,N,Nsq,h,ERI_MO,rdm1,rdm2,hess_tmp)
    hess = hess + hess_tmp 
   ! if(dRPA) then
   !   call orbital_hessian(O,V,N,Nsq,h,ERI_MO,rdm1_c,rdm2_c,hess_tmp)
   ! else
   !   call orbital_hessian(O,V,N,Nsq,F,ERI_MO_AS,rdm1_c,rdm2_c,hess_tmp)
   ! endif
   ! hess = hess + hess_tmp
   ! write(*,*) "Hessian from rdms"
   ! call matout(Nsq,Nsq,hess)
    deallocate(hess_tmp)
    allocate(hessInv(nhess,mhess))
    call pseudo_inverse_matrix(Nsq,hess,hessInv)
   ! allocate(VL(Nsq,Nsq),VR(Nsq,Nsq),WR(Nsq))
   ! call diagonalize_general_matrix_LR(Nsq,hess,WR,VL,VR)
   ! call vecout(Nsq,WR)
   ! deallocate(VL,VR,WR)
    deallocate(hess)
  else
    nhess = O*V
    mhess = nhess
    allocate(hess(nhess,mhess),hess_tmp(nhess,mhess))
    hess(:,:) = 0d0
    call wall_time(tstart)
    call orbital_hessian_ov(O,V,N,Nsq,h,ERI_MO,rdm1_hf,rdm2_hf,hess_tmp)
    !call orbital_hessian_ov(O,V,N,Nsq,h,ERI_MO,rdm1,rdm2,hess_tmp)
    hess = hess + hess_tmp 
    if(dRPA) then
      call orbital_hessian_ov(O,V,N,Nsq,h,ERI_MO,rdm1_c,rdm2_c,hess_tmp)
    else
      call orbital_hessian_ov(O,V,N,Nsq,F,ERI_MO_AS,rdm1_c,rdm2_c,hess_tmp)
    endif
    call wall_time(tend)
    tdiff = tend - tstart
    write(*,*) "Building Hessian took ", tdiff, "seconds."
    hess = hess + hess_tmp 
    deallocate(hess_tmp)
    allocate(hessInv(nhess,mhess))
    call wall_time(tstart)
    call pseudo_inverse_general_matrix(O*V,hess,hessInv)
    call wall_time(tend)
    tdiff = tend - tstart
    write(*,*) "Inverting Hessian took ", tdiff, "seconds."
   ! allocate(VL(O*V,O*V),VR(O*V,O*V),WR(O*V))
   ! call diagonalize_general_matrix_LR(O*V,hess,WR,VL,VR)
   ! call vecout(O*V,WR)
   ! deallocate(VL,VR,WR)
    deallocate(hess)
  endif
  if((.not. OVRotOnly) .or. diagHess) then
    
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
  
  else
    allocate(Kap(N,N))
    
    Kap(:,:) = 0d0 
    
    do i=1,O
      do a=O+1,N
   
        ia = a - O + (i-1)*V
        do j=1,O
          do b=O+1,N 
            jb = b - O + (j-1)*V
            Kap(i,a) = Kap(i,a) - hessInv(ia,jb)*grad(b + (j-1)*N)
          end do
        end do
        Kap(a,i) = - Kap(i,a)
      end do
    end do
   
  end if
 
  call rotate_orbitals(N,Kap,c)

  ! Transform integrals and Hc
  call AOtoMO(nBas,nOrb,c,Hc,h)
  call AOtoMO_ERI_RHF(nBas,nOrb,c,ERI_AO,ERI_MO)
  if(.not. dRPA) call AOtoMO_ERI_RHF(nBas,nOrb,c,ERI_AO_AS,ERI_MO_AS)
deallocate(rdm1,rdm2)
end subroutine
