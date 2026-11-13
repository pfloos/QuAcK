subroutine RGW_iw_Z(flow,nOrb,nC,nO,nV,nR,eHF,ERI,nfreqs,Z)

! Quasiparticle renormalization factor  Z_p = 1/(1 - Re dSigma_c,pp/dw |_{w=e_p})  from imaginary-axis data only
! (no poles of W). Sigma_c,pp(z) is evaluated by contour deformation exactly as in RGW_SRG_iw_self_energy
! (energy-specific least-squares weights on the greenX grid + residues of the continued W) and differentiated
! by a central finite difference along the real direction, at z = e_p +- h + i dZ (defaults h = 0.002, dZ = 0.01 Ha;
! GWIW_ZH / GWIW_ZDELTA override). A large dZ (e.g. 1/sqrt(s)) smears Sigma_c and biases Z (H2: 0.64 vs 0.51 for orbital 7).
!
! Intended to be called once, at the end of the SCF cycles. Frozen core is not implemented (nC = nR = 0).

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: flow
  integer,intent(in)            :: nOrb,nC,nO,nV,nR,nfreqs
  double precision,intent(in)   :: eHF(nOrb)
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)

! Output variables

  double precision,intent(out)  :: Z(nOrb)

! Local variables

  integer,parameter             :: nTest = 600
  integer,parameter             :: nblkmax = 400

  integer                       :: n2,p,q,r,k,iref,nref,nrow,nb,is
  double precision              :: mu,gap,delta,Omin,Omax,Rpool,D,h,dZ,wall0,wall1
  character(len=32)             :: envval
  integer                       :: envstat
  double precision,allocatable  :: e(:),wk(:),Omtest(:),P_lsq(:,:),vMAT(:,:),coef(:,:,:,:)
  integer,allocatable           :: nblk(:,:)
  complex(kind=8),allocatable   :: omref(:),Wref(:,:,:,:),Fref(:),Wtmp(:,:,:),zk(:)
  complex(kind=8)               :: zz(2),Sg(2),a
  complex(kind=8),external      :: GWiw_AC_sym_even_eval

  if(nC /= 0 .or. nR /= 0) error stop 'RGW_iw_Z: frozen core / removed virtuals not implemented'
  if(nfreqs < 4) error stop 'RGW_iw_Z: nfreqs too small'

  call wall_time(wall0)

  n2 = nOrb*nOrb

  allocate(e(nOrb))
  mu  = 0.5d0*(eHF(nO)+eHF(nO+1))
  e   = eHF - mu
  gap = eHF(nO+1) - eHF(nO)

  delta = 1d0/sqrt(flow)
  Omin  = 0.5d0*gap
  Omax  = 1.3d0*(maxval(e)-minval(e))
  Rpool = max(30d0,1.5d0*(maxval(e)-minval(e)))
  D     = delta
  nrow  = int(Rpool/D)
  nref  = nfreqs + nrow
  h     = 0.002d0
  dZ    = 0.01d0
  call get_environment_variable('GWIW_ZH',envval,status=envstat)
  if(envstat == 0) read(envval,*) h
  call get_environment_variable('GWIW_ZDELTA',envval,status=envstat)
  if(envstat == 0) read(envval,*) dZ
  write(*,'(A,2F10.5)') '   Finite-difference step h and contour height for Z = ',h,dZ

  write(*,*)
  write(*,*) '*** Z factors from imaginary-axis data (finite difference of Sigma_c(z)) ***'
  write(*,*)

  allocate(vMAT(n2,n2))
  do p=1,nOrb
    do q=1,nOrb
      do r=1,nOrb
        do k=1,nOrb
          vMAT(1+(r-1)+(p-1)*nOrb,1+(k-1)+(q-1)*nOrb) = ERI(p,q,r,k)
        end do
      end do
    end do
  end do

  allocate(wk(nfreqs),Omtest(nTest),P_lsq(nfreqs,nTest))
  call GWiw_read_wcoord(nfreqs,wk)
  call GWiw_LSQ_setup(nfreqs,wk,Omin,Omax,nTest,Omtest,P_lsq)

  allocate(omref(nref),Wref(nOrb,nOrb,nOrb,nref))
  do k=1,nfreqs
    omref(k) = cmplx(0d0,wk(k),kind=8)
  end do
  do k=1,nrow
    omref(nfreqs+k) = cmplx(D*k,delta,kind=8)
  end do

  allocate(Wtmp(nOrb,nOrb,nOrb))
  !$OMP PARALLEL DO PRIVATE(iref,Wtmp) SHARED(Wref,omref,nref,nOrb,nO,e,vMAT) DEFAULT(NONE) SCHEDULE(DYNAMIC)
  do iref=1,nref
    call GWiw_Wc_MO_RHF(nOrb,nO,e,omref(iref),vMAT,Wtmp)
    Wref(:,:,:,iref) = Wtmp(:,:,:)
  end do
  !$OMP END PARALLEL DO
  deallocate(Wtmp)

! Continued-fraction fits of the diagonal elements W_pp;r(w)

  allocate(coef(5,nblkmax,nOrb,nOrb),nblk(nOrb,nOrb))
  coef = 0d0
  nblk = 0
  !$OMP PARALLEL DO PRIVATE(r,p,Fref,nb) SHARED(Wref,omref,coef,nblk,nref,nOrb) DEFAULT(NONE) SCHEDULE(DYNAMIC) COLLAPSE(2)
  do r=1,nOrb
    do p=1,nOrb
      allocate(Fref(nref))
      Fref(:) = Wref(p,p,r,:)
      call GWiw_AC_sym_even_coeff(nref,omref,Fref,nblkmax,nb,coef(:,:,p,r))
      nblk(p,r) = nb
      deallocate(Fref)
    end do
  end do
  !$OMP END PARALLEL DO

! Z_p from dSigma_c,pp/dw at w = e_p

  Z = 0d0

  !$OMP PARALLEL DO PRIVATE(p,is,r,k,a,zz,Sg,zk) SHARED(nOrb,nO,nfreqs,e,dZ,h,P_lsq,Omtest,Wref,coef,nblk,Z) &
  !$OMP DEFAULT(NONE) SCHEDULE(DYNAMIC)
  do p=1,nOrb
    allocate(zk(nfreqs))
    zz(1) = cmplx(e(p)-h,dZ,kind=8)
    zz(2) = cmplx(e(p)+h,dZ,kind=8)
    Sg    = (0d0,0d0)
    do is=1,2
      do r=1,nOrb
        a = zz(is) - e(r)
        call GWiw_LSQ_weights(nfreqs,nTest,P_lsq,Omtest,a,zk)
        do k=1,nfreqs
          Sg(is) = Sg(is) + zk(k)*Wref(p,p,r,k)
        end do
        if(r <= nO .and. real(a) < 0d0) then
          Sg(is) = Sg(is) - GWiw_AC_sym_even_eval(nblk(p,r),coef(1,1,p,r),e(r)-zz(is))
        else if(r > nO .and. real(a) >= 0d0) then
          Sg(is) = Sg(is) + GWiw_AC_sym_even_eval(nblk(p,r),coef(1,1,p,r),zz(is)-e(r))
        end if
      end do
    end do
    Z(p) = 1d0/(1d0 - real(Sg(2)-Sg(1))/(2d0*h))
    deallocate(zk)
  end do
  !$OMP END PARALLEL DO

  call wall_time(wall1)
  write(*,'(A,F10.3,A)') '   Wall time for the Z factors = ',wall1-wall0,' seconds'
  write(*,*)

  deallocate(e,vMAT,wk,Omtest,P_lsq,omref,Wref,coef,nblk)

end subroutine
