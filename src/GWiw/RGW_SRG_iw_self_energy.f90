subroutine RGW_SRG_iw_self_energy(flow,nOrb,nC,nO,nV,nR,eHF,ERI,nfreqs,SigC)

! SRG static self-energy of Marie and Loos (JCTC 19, 3943 (2023)), restricted case, evaluated WITHOUT the poles
! of W, from imaginary-frequency data only:
!
!   St_pq(s) = sum_{r,nu} M_pr,nu M_qr,nu K_s(e_p - w_r,nu, e_q - w_r,nu)
!            = -(1/pi) Im  int dx  K_s(e_p - z, e_q - z)  Sigma_c,pq(z) ,   z = x + i delta ,  delta = 1/sqrt(s)
!
! (K_s is entire, so the integral over the real axis is deformed to the line Im z = delta: Cauchy's theorem).
! Sigma_c(z) is obtained by contour deformation (Duchemin and Blase, JCTC 16, 1742 (2020)):
!   - imaginary-axis integral G(z+iw) W_c(iw) with energy-specific least-squares quadrature weights
!     on the greenX imaginary frequencies (./grids/wcoord.txt),
!   - residue terms W_c(e_r - z) obtained by analytic continuation of W_c (symmetric-even continued fraction)
!     from reference frequencies on the imaginary axis and on a row parallel to the real axis at height delta.
! The chemical potential is placed in the middle of the HOMO-LUMO gap.
!
! Frozen core / removed virtual orbitals are not implemented (nC = nR = 0 is required).

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: flow
  integer,intent(in)            :: nOrb,nC,nO,nV,nR,nfreqs
  double precision,intent(in)   :: eHF(nOrb)
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)

! Output variables

  double precision,intent(out)  :: SigC(nOrb,nOrb)

! Local variables

  integer,parameter             :: nTest = 600
  integer,parameter             :: nblkmax = 400
  integer,parameter             :: nnodes_max = 3000000

  integer                       :: n2,p,q,r,k,iref,nref,nrow,nnodes,j,nb
  double precision              :: mu,gap,delta,Omin,Omax,Rpool,D,wall0,wall1,dfac
  character(len=32)             :: envval
  integer                       :: envstat
  double precision,allocatable  :: e(:),wk(:),Omtest(:),P_lsq(:,:),vMAT(:,:),xnode(:),wnode(:),coef(:,:,:,:,:),St(:,:)
  integer,allocatable           :: nblk(:,:,:)
  complex(kind=8),allocatable   :: omref(:),Wref(:,:,:,:),Fref(:),Wtmp(:,:,:),Sig(:,:),zk(:)
  complex(kind=8)               :: z,a,Wac
  complex(kind=8),external      :: GWiw_AC_sym_even_eval,GWiw_Ks

  if(nC /= 0 .or. nR /= 0) error stop 'RGW_SRG_iw_self_energy: frozen core / removed virtuals not implemented'
  if(nfreqs < 4) error stop 'RGW_SRG_iw_self_energy: nfreqs too small'

  call wall_time(wall0)

  write(*,*)
  write(*,*) '*** SRG self-energy from imaginary-axis data (contour deformation) ***'
  write(*,*)

  n2 = nOrb*nOrb

! Energies relative to the middle of the HOMO-LUMO gap

  allocate(e(nOrb))
  mu  = 0.5d0*(eHF(nO)+eHF(nO+1))
  e   = eHF - mu
  gap = eHF(nO+1) - eHF(nO)

! Parameters

  delta = 1d0/sqrt(flow)
  Omin  = 0.5d0*gap
  Omax  = 1.3d0*(maxval(e)-minval(e))
  Rpool = max(30d0,1.5d0*(maxval(e)-minval(e)))
  dfac  = 3d0
  call get_environment_variable('GWIW_DFAC',envval,status=envstat)
  if(envstat == 0) read(envval,*) dfac
  D     = dfac*delta
  nrow  = int(Rpool/D)
  nref  = nfreqs + nrow

  write(*,'(A,F12.5)')  '   SRG flow parameter s            = ',flow
  write(*,'(A,F12.5)')  '   Contour height delta = 1/sqrt(s) = ',delta
  write(*,'(A,2F12.5)') '   Quadrature range [Omin,Omax] for the poles of W = ',Omin,Omax
  write(*,'(A,I8)')     '   Imaginary frequencies           = ',nfreqs
  write(*,'(A,I8)')     '   Additional W reference points   = ',nrow
  write(*,'(A,F8.2)')   '   Reference spacing / delta       = ',dfac

! Two-electron integrals as a supermatrix  v(1+(r-1)+(p-1)n, 1+(s-1)+(q-1)n) = <pq|rs>

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

! Quadrature weights setup

  allocate(wk(nfreqs),Omtest(nTest),P_lsq(nfreqs,nTest))
  call GWiw_read_wcoord(nfreqs,wk)
  call GWiw_LSQ_setup(nfreqs,wk,Omin,Omax,nTest,Omtest,P_lsq)

! Reference frequencies for the continuation of W: imaginary axis, then a row at height delta

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

! Analytic continuation fits of W_pq;r(w) for p <= q

  allocate(coef(5,nblkmax,nOrb,nOrb,nOrb),nblk(nOrb,nOrb,nOrb))
  coef = 0d0
  nblk = 0
  !$OMP PARALLEL DO PRIVATE(r,p,q,Fref,nb) SHARED(Wref,omref,coef,nblk,nref,nOrb) DEFAULT(NONE) SCHEDULE(DYNAMIC) COLLAPSE(2)
  do r=1,nOrb
    do p=1,nOrb
      allocate(Fref(nref))
      do q=p,nOrb
        Fref(:) = Wref(p,q,r,:)
        call GWiw_AC_sym_even_coeff(nref,omref,Fref,nblkmax,nb,coef(:,:,p,q,r))
        nblk(p,q,r) = nb
      end do
      deallocate(Fref)
    end do
  end do
  !$OMP END PARALLEL DO

  write(*,'(A,I8)') '   Maximum number of continued-fraction blocks = ',maxval(nblk)

! Contour nodes

  allocate(xnode(nnodes_max),wnode(nnodes_max))
  call GWiw_SRG_contour_nodes(nOrb,e,flow,delta,Omax,nnodes_max,nnodes,xnode,wnode)
  write(*,'(A,I8)') '   Contour nodes                  = ',nnodes

! Contour integral

  allocate(St(nOrb,nOrb))
  St = 0d0

  !$OMP PARALLEL PRIVATE(j,z,a,Sig,zk,r,p,q,k,Wac) SHARED(nnodes,xnode,wnode,delta,e,nOrb,nO,nfreqs,P_lsq,Omtest,Wref,coef,nblk,flow) &
  !$OMP REDUCTION(+:St) DEFAULT(NONE)
  allocate(Sig(nOrb,nOrb),zk(nfreqs))
  !$OMP DO SCHEDULE(STATIC)
  do j=1,nnodes
    z = cmplx(xnode(j),delta,kind=8)
    Sig = (0d0,0d0)
    do r=1,nOrb
      a = z - e(r)
      call GWiw_LSQ_weights(nfreqs,nTest,P_lsq,Omtest,a,zk)
      do q=1,nOrb
        do p=1,q
          do k=1,nfreqs
            Sig(p,q) = Sig(p,q) + zk(k)*Wref(p,q,r,k)
          end do
        end do
      end do
      if(r <= nO .and. real(a) < 0d0) then
        do q=1,nOrb
          do p=1,q
            Sig(p,q) = Sig(p,q) - GWiw_AC_sym_even_eval(nblk(p,q,r),coef(1,1,p,q,r),e(r)-z)
          end do
        end do
      else if(r > nO .and. real(a) >= 0d0) then
        do q=1,nOrb
          do p=1,q
            Sig(p,q) = Sig(p,q) + GWiw_AC_sym_even_eval(nblk(p,q,r),coef(1,1,p,q,r),z-e(r))
          end do
        end do
      end if
    end do
    do q=1,nOrb
      do p=1,q
        St(p,q) = St(p,q) - aimag(GWiw_Ks(e(p)-z,e(q)-z,flow)*Sig(p,q))*wnode(j)/pi
      end do
    end do
  end do
  !$OMP END DO
  deallocate(Sig,zk)
  !$OMP END PARALLEL

  do q=1,nOrb
    do p=1,q
      SigC(p,q) = St(p,q)
      SigC(q,p) = St(p,q)
    end do
  end do

  call wall_time(wall1)
  write(*,'(A,F10.3,A)') '   Wall time for the imaginary-axis SRG self-energy = ',wall1-wall0,' seconds'
  write(*,*)

  deallocate(e,vMAT,wk,Omtest,P_lsq,omref,Wref,coef,nblk,xnode,wnode,St)

end subroutine
