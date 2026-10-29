program gw_parquet_checks

! Algebraic regressions for frozen external frequencies, derivatives,
! GW-kernel consistency, pair normalization, and active-space indexing.

  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer,parameter :: no=6,nc=1,nocc=3,nv=3,nr=1,ns=4,nos=3,nvs=3,not=1,nvt=1
  integer :: p,q,r,s,n,spin,mode,i,j,a,b,ia,jb,nh,np
  double precision :: eri(no,no,no,no),w0(no,no,no,no),phi(no,no,no,no),e(no),omega(no)
  double precision :: es(ns),et(ns),ps(nvs),hs(nos),pt(nvt),ht(not)
  double precision :: res(no,no,ns),ret(no,no,ns),rps(no,no,nvs),rhs(no,no,nos)
  double precision :: rpt(no,no,nvt),rht(no,no,not),rw(no,no,ns)
  double precision :: sig(no),z(no),sp(no),sm(no),ztmp(no),sig_old(no),z_old(no),ec,fd,expected
  double precision :: xp(ns,ns),xm(ns,ns),rho(no,no,ns),ref(no,no,ns)
  double precision :: aa(ns,ns),bb(ns,ns),ka(ns,ns),kb(ns,ns)
  double precision,allocatable :: c0(:,:),d0(:,:),b0(:,:),c1(:,:),d1(:,:),b1(:,:)
  double precision,allocatable :: kc(:,:),kd(:,:),kpp(:,:),x1(:,:),y1(:,:),x2(:,:),y2(:,:)
  double precision,allocatable :: rr1(:,:,:),rr2(:,:,:),rf1(:,:,:),rf2(:,:,:)
  logical :: reg
  double precision,parameter :: eta=0.7d0,h=1d-6

  e=[-2d0,-1d0,-0.6d0,0.4d0,0.8d0,1.4d0]
  es=[0.8d0,1.1d0,1.4d0,1.7d0]
  et=es+0.13d0
  ps=[0.9d0,1.2d0,1.5d0]
  hs=[-1.2d0,-1.5d0,-1.9d0]
  pt=1.3d0
  ht=-1.4d0
  do s=1,no
    do r=1,no
      do q=1,no
        do p=1,no
          eri(p,q,r,s)=0.001d0*(sin(dble(p+r))*sin(dble(q+s)) &
                                     +cos(dble(p+r))*cos(dble(q+s)))
        end do
      end do
    end do
  end do
  do n=1,ns
    do q=1,no
      do p=1,no
        res(p,q,n)=0.03d0*sin(dble(p+2*q+n))
        ret(p,q,n)=0.02d0*cos(dble(2*p+q+n))
        rw(p,q,n)=0.01d0*sin(dble(p+q+n))
      end do
    end do
  end do
  rps=res(:,:,1:nvs)
  rhs=ret(:,:,1:nos)
  rpt=res(:,:,1:nvt)
  rht=ret(:,:,1:not)

! The old entry point and the generalized entry preserve legacy results.
  omega=e
  call R_Parquet_self_energy_alt_diag(eta,no,nc,nocc,nv,nr,ns,nos,nvs,not,nvt,e,eri, &
        res,es,ret,et,rps,ps,rpt,pt,rhs,hs,rht,ht,ec,sig_old,z_old)
  call R_Parquet_self_energy_alt_omega(.true.,.false.,0,eta,omega, &
        no,nc,nocc,nv,nr,ns,nos,nvs,not,nvt,e,eri, &
        res,es,ret,et,rps,ps,rpt,pt,rhs,hs,rht,ht,ec,sig,z)
  call check(maxval(abs(sig-sig_old)) < 1d-13,'legacy self-energy entry point')
  call check(maxval(abs(z-z_old)) < 1d-13,'legacy derivative convention')

! Check exact derivatives of the frequency-dependent expression, including SRG.
  do mode=1,2
    reg=(mode == 1)
    omega=e+0.127d0
    call evaluate(0,sig,z)
    do p=nc+1,no-nr
      omega(p)=omega(p)+h
      call evaluate(p,sp,ztmp)
      omega(p)=omega(p)-2d0*h
      call evaluate(p,sm,ztmp)
      omega(p)=omega(p)+h
      fd=(sp(p)-sm(p))/(2d0*h)
      call check(abs(fd-(1d0-1d0/z(p))) < 1d-8,'external-frequency derivative')
      call check(all(abs(sm(1:nc)) < 1d-15),'inactive orbitals remain inactive')
    end do
    call evaluate(nc+1,sp,ztmp)
    call check(abs(sig(nc+1)-sp(nc+1)) < 1d-13,'single-orbital evaluator')
    call check(all(abs(sp(nc+2:no)) < 1d-15),'single-orbital selection')
  end do
  reg=.true.
  omega=e
  omega(2)=e(2)+e(3)-e(4)
  call evaluate(2,sig,z)
  call check(all(ieee_is_finite(sig)) .and. all(ieee_is_finite(z)),'finite SRG zero-denominator limit')

! Bare-kernel limit agrees with existing parquet half-vertex contractions,
! including backward amplitudes and singlet unique-pair normalization.
  phi=0d0
  do j=1,ns
    do i=1,ns
      xp(i,j)=cos(dble(i+2*j))
      xm(i,j)=sin(dble(2*i+j))
    end do
  end do
  do spin=1,2
    call R_GW_Parquet_eh_integrals(spin,no,nc,nocc,nr,ns,eri,eri,xp,xm,rho)
    if(spin == 1) then
      call R_eh_singlet_screened_integral(no,nc,nocc,nr,ns,eri,phi,phi,phi,phi,xp,xm,ref)
    else
      call R_eh_triplet_screened_integral(no,nc,nocc,nr,ns,eri,phi,phi,phi,phi,xp,xm,ref)
    end if
    call check(maxval(abs(rho-ref)) < 1d-13,'bare eh vertex normalization and orientation')
  end do

! Build a nonzero W0 and compare its pair projection with the existing
! GW ppBSE kernels, independently of the new external coupling routines.
  w0=eri
  do s=1,no
    do r=1,no
      do q=1,no
        do p=1,no
          do n=1,ns
            w0(p,q,r,s)=w0(p,q,r,s)-4d0*rw(p,r,n)*rw(q,s,n)/es(n)
          end do
        end do
      end do
    end do
  end do
  do spin=1,2
    nh=nos
    np=nvs
    if(spin == 2) then
      nh=not
      np=nvt
    end if
    allocate(c0(np,np),c1(np,np),d0(nh,nh),d1(nh,nh),b0(np,nh),b1(np,nh))
    allocate(kc(np,np),kd(nh,nh),kpp(np,nh),x1(np,np),y1(nh,np),x2(np,nh),y2(nh,nh))
    allocate(rr1(no,no,np),rr2(no,no,nh),rf1(no,no,np),rf2(no,no,nh))
    call ppRLR_C(spin,no,nc,nocc,nv,nr,np,1d0,e,eri,c0)
    call ppRLR_D(spin,no,nc,nocc,nv,nr,nh,1d0,e,eri,d0)
    call ppRLR_B(spin,no,nc,nocc,nv,nr,nh,np,1d0,eri,b0)
    call ppRLR_C(spin,no,nc,nocc,nv,nr,np,1d0,e,w0,c1)
    call ppRLR_D(spin,no,nc,nocc,nv,nr,nh,1d0,e,w0,d1)
    call ppRLR_B(spin,no,nc,nocc,nv,nr,nh,np,1d0,w0,b1)
    call RGW_ppBSE_static_kernel_C(spin,0d0,no,nc,nocc,nv,nr,ns,np,1d0,eri,es,rw,kc)
    call RGW_ppBSE_static_kernel_D(spin,0d0,no,nc,nocc,nv,nr,ns,nh,1d0,eri,es,rw,kd)
    call RGW_ppBSE_static_kernel_B(spin,0d0,no,nc,nocc,nv,nr,ns,nh,np,1d0,eri,es,rw,kpp)
    call check(maxval(abs(c1-c0-kc)) < 1d-13,'GW ppBSE C kernel')
    call check(maxval(abs(d1-d0-kd)) < 1d-13,'GW ppBSE D kernel')
    call check(maxval(abs(b1-b0-kpp)) < 1d-13,'GW ppBSE B kernel')
    x1=0.2d0
    y1=0.03d0
    x2=-0.04d0
    y2=0.3d0
    call R_GW_Parquet_pp_integrals(spin,no,nc,nocc,nr,nh,np,eri,x1,y1,x2,y2,rr1,rr2)
    if(spin == 1) then
      call R_pp_singlet_screened_integral(no,nc,nocc,nr,nh,np,eri,phi,phi,x1,y1,rf1,x2,y2,rf2)
    else
      call R_pp_triplet_screened_integral(no,nc,nocc,nr,nh,np,eri,phi,phi,x1,y1,rf1,x2,y2,rf2)
    end if
    call check(maxval(abs(rr1-rf1)) < 1d-13,'bare ee vertex normalization')
    call check(maxval(abs(rr2-rf2)) < 1d-13,'bare hh vertex normalization')
    deallocate(c0,c1,d0,d1,b0,b1,kc,kd,kpp,x1,y1,x2,y2,rr1,rr2,rf1,rf2)
  end do

  do spin=1,2
    call phRLR_A(spin,.false.,no,nc,nocc,nv,nr,ns,1d0,e,eri,aa)
    call phRLR_B(spin,.false.,no,nc,nocc,nv,nr,ns,1d0,eri,bb)
    call RGW_phBSE_static_kernel_A(0d0,no,nc,nocc,nv,nr,ns,1d0,eri,es,rw,ka)
    call RGW_phBSE_static_kernel_B(0d0,no,nc,nocc,nv,nr,ns,1d0,eri,es,rw,kb)
    xp=0d0
    xm=0d0
    do i=1,ns
      xp(i,i)=1d0
      xm(i,i)=-1d0
    end do
    call R_GW_Parquet_eh_integrals(spin,no,nc,nocc,nr,ns,eri,w0,xp,xm,rho)
    xm=xp
    call R_GW_Parquet_eh_integrals(spin,no,nc,nocc,nr,ns,eri,w0,xp,xm,ref)
    ia=0
    do i=nc+1,nocc
      do a=nocc+1,no-nr
        ia=ia+1
        jb=0
        do j=nc+1,nocc
          do b=nocc+1,no-nr
            jb=jb+1
            call check(abs(ka(ia,jb)-eri(b,j,a,i)+w0(b,j,a,i)) < 1d-13,'GW ehBSE A kernel')
            call check(abs(kb(ia,jb)-eri(b,j,i,a)+w0(b,j,i,a)) < 1d-13,'GW ehBSE B kernel')
            expected=aa(ia,jb)+ka(ia,jb)
            if(i == j .and. a == b) expected=expected-e(a)+e(i)
            call check(abs(rho(b,j,ia)-expected) < 1d-13,'eh half vertex matches GW A block')
            call check(abs(ref(b,j,ia)-bb(ia,jb)-kb(ia,jb)) < 1d-13,'eh half vertex matches GW B block')
          end do
        end do
      end do
    end do
  end do
  print *, 'GW-PARQUET ALGEBRAIC CHECKS PASSED'

contains
  subroutine evaluate(orb,sv,zv)
    integer,intent(in) :: orb
    double precision,intent(out) :: sv(no),zv(no)
    call R_Parquet_self_energy_alt_omega(reg,.true.,orb,eta,omega,no,nc,nocc,nv,nr,ns,nos,nvs,not,nvt,e,eri, &
          res,es,ret,et,rps,ps,rpt,pt,rhs,hs,rht,ht,ec,sv,zv)
  end subroutine evaluate
  subroutine check(ok,message)
    logical,intent(in) :: ok
    character(len=*),intent(in) :: message
    if(.not.ok) then
      print *, 'FAILED: ',message
      error stop 1
    end if
  end subroutine check
end program gw_parquet_checks
