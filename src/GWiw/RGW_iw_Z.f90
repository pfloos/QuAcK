subroutine RGW_iw_Z(flow,nOrb,nC,nO,nV,nR,eHF,ERI,nfreqs,Z)

! SRG quasiparticle renormalization factor of Marie and Loos (restricted case), from imaginary-axis data only:
!
!   1/Z_p = 1 + sum_{r,nu} M_pr,nu^2 g_s(e_p - w_r,nu) ,   g_s(D) = (1 - exp(-2 s D^2))/D^2   (entire function)
!         = 1 - (1/pi) Im  int dx  g_s(e_p - z) Sigma_c,pp(z) ,   z = x + i delta
!
! (same contour representation as the static SRG self-energy, with the kernel K_s(D,D) replaced by g_s).
! Sigma_c,pp(z) is evaluated by contour deformation exactly as in RGW_SRG_iw_self_energy (energy-specific
! least-squares weights on the greenX grid + residues of the continued W). No finite differences.
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
  integer                       :: n2,p,q,r,k,iref,nref,nrow,nb

  integer,parameter             :: nnodes_max = 3000000
  integer                       :: j,nnodes
  double precision              :: mu,gap,delta,Omin,Omax,Rpool,D,wall0,wall1
  double precision,allocatable  :: xnode(:),wnode(:),S(:)
  complex(kind=8)               :: zc,Dl,u,gz
  complex(kind=8),allocatable   :: Sg(:)
  double precision,allocatable  :: e(:),wk(:),Omtest(:),P_lsq(:,:),vMAT(:,:),coef(:,:,:,:)
  integer,allocatable           :: nblk(:,:)
  complex(kind=8),allocatable   :: omref(:),Wref(:,:,:,:),Fref(:),Wtmp(:,:,:),zk(:)
  complex(kind=8)               :: a
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
  write(*,*)
  write(*,*) '*** SRG Z factors from imaginary-axis data (contour integral) ***'
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

! Contour integral

  allocate(xnode(nnodes_max),wnode(nnodes_max))
  call GWiw_SRG_contour_nodes(nOrb,e,flow,delta,Omax,nnodes_max,nnodes,xnode,wnode)
  write(*,'(A,I8)') '   Contour nodes                  = ',nnodes

  allocate(S(nOrb))
  S = 0d0

  !$OMP PARALLEL PRIVATE(j,zc,a,Sg,zk,r,p,k,Dl,u,gz) SHARED(nnodes,xnode,wnode,delta,e,nOrb,nO,nfreqs,P_lsq,Omtest,Wref,coef,nblk,flow) &
  !$OMP REDUCTION(+:S) DEFAULT(NONE)
  allocate(Sg(nOrb),zk(nfreqs))
  !$OMP DO SCHEDULE(STATIC)
  do j=1,nnodes
    zc = cmplx(xnode(j),delta,kind=8)
    Sg = (0d0,0d0)
    do r=1,nOrb
      a = zc - e(r)
      call GWiw_LSQ_weights(nfreqs,nTest,P_lsq,Omtest,a,zk)
      do p=1,nOrb
        do k=1,nfreqs
          Sg(p) = Sg(p) + zk(k)*Wref(p,p,r,k)
        end do
        if(r <= nO .and. real(a) < 0d0) then
          Sg(p) = Sg(p) - GWiw_AC_sym_even_eval(nblk(p,r),coef(1,1,p,r),e(r)-zc)
        else if(r > nO .and. real(a) >= 0d0) then
          Sg(p) = Sg(p) + GWiw_AC_sym_even_eval(nblk(p,r),coef(1,1,p,r),zc-e(r))
        end if
      end do
    end do
    do p=1,nOrb
      Dl = e(p) - zc
      u  = 2d0*flow*Dl*Dl
      if(abs(u) < 1d-6) then
        gz = 2d0*flow*(1d0 - u/2d0 + u*u/6d0)
      else
        gz = (1d0 - exp(-u))/(Dl*Dl)
      end if
      S(p) = S(p) - aimag(gz*Sg(p))*wnode(j)/pi
    end do
  end do
  !$OMP END DO
  deallocate(Sg,zk)
  !$OMP END PARALLEL

  Z(:) = 1d0/(1d0 + S(:))

  call wall_time(wall1)
  write(*,'(A,F10.3,A)') '   Wall time for the Z factors = ',wall1-wall0,' seconds'
  write(*,*)

  deallocate(e,vMAT,wk,Omtest,P_lsq,omref,Wref,coef,nblk,xnode,wnode,S)

end subroutine
