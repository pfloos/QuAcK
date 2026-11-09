subroutine RGW_COT(nOrb,nC,nO,nR,nS,eHF,Om,rho,hCOT,hHerm,eCOT,info)

! One-shot GW-based correlated-orbital theory (COT), restricted closed shell.
!
! This is the dense, reference implementation: diagonalize the Hermitian
! upfolded G0W0 supermatrix, select one quasiparticle root per active orbital,
! and construct the Bloch effective one-particle Hamiltonian.
!
! Sigma_c(p,q;w) = 2 sum_{im} rho(p,i,m)rho(q,i,m)/(w-e_i+Om_m)
!                 + 2 sum_{am} rho(p,a,m)rho(q,a,m)/(w-e_a-Om_m)
!
! H = [ A   B ]        A = diag(eHF)
!     [ B^T D ]        D = diag(e_i-Om_m,e_a+Om_m)
!
! Given selected eigenvectors R=[RP;RQ], S=RQ RP^{-1}, and
! hCOT=A+B*S. Its Hermitian representative is
! hHerm=N^{-1/2} (Omega^T H Omega) N^{-1/2}, N=I+S^T S.
!
! This is a one-shot G0W0 calculation: no qsGW orbital updates, no extra
! density/static correction, and no broadening/SRG regularization.
!
! The active orbitals are nC+1 : nOrb-nR. Frozen orbitals remain diagonal.
! Dense diagonalization scales as O((nAct*(nS+1))^3) and should be used
! only as a numerical reference for modest systems.
!
! info = 0 success
!        1 invalid dimensions
!        2 dense problem exceeds safety limit
!        3 eigensolver/SVD/linear solve failed
!        4 selected one-particle projection ill-conditioned
!        5 invariance/eigenvalue residual failed

  implicit none
  include 'parameters.h'

  integer,intent(in)            :: nOrb,nC,nO,nR,nS
  double precision,intent(in)   :: eHF(nOrb),Om(nS)
  double precision,intent(in)   :: rho(nOrb,nOrb,nS)
  double precision,intent(out)  :: hCOT(nOrb,nOrb),hHerm(nOrb,nOrb)
  double precision,intent(out)  :: eCOT(nOrb)
  integer,intent(out)           :: info

  integer,parameter             :: maxDense = 4000
  double precision,parameter    :: condLimit = 1d10
  double precision,parameter    :: residualTol = 1d-7

  integer                       :: nAct,nOcc,nVir,nAux,nH
  integer                       :: p,q,k,mu,i,a,m,j,ibest,pbest,istat,lwork
  integer,allocatable           :: picked(:),ipiv(:)
  logical,allocatable           :: usedState(:),usedOrb(:)
  double precision              :: best,val,condition,resid,ediff,swap
  double precision,allocatable  :: H(:,:),E(:),B(:,:),D(:)
  double precision,allocatable  :: RP(:,:),RQ(:,:),RPwork(:,:),rhs(:,:),S(:,:)
  double precision,allocatable  :: sig(:),work(:),N(:,:),eN(:),Nhalf(:,:)
  double precision,allocatable  :: NormVec(:,:),hProj(:,:),hBlock(:,:),bs(:,:),ds(:,:)
  double precision,allocatable  :: eCheck(:),eSelected(:),eigWork(:,:),rTest(:,:)
  double precision              :: dummyU(1,1),dummyVT(1,1)

  info=0
  hCOT=0d0
  hHerm=0d0
  eCOT=eHF
  do p=1,nOrb
    hCOT(p,p)=eHF(p)
    hHerm(p,p)=eHF(p)
  end do

  nAct=nOrb-nC-nR
  nOcc=nO-nC
  nVir=nOrb-nR-nO
  if(nAct < 1 .or. nOcc < 0 .or. nVir < 0 .or. nS < 0) then
    info=1
    return
  end if
  if(nS == 0) return

  nAux=(nOcc+nVir)*nS
  nH=nAct+nAux
  if(nH > maxDense) then
    write(*,'(A,I0,A,I0)') 'GW-COT: dense dimension ',nH,' exceeds limit ',maxDense
    info=2
    return
  end if

  allocate(H(nH,nH),E(nH),B(nAct,nAux),D(nAux),picked(nAct))
  allocate(usedState(nH),usedOrb(nAct))
  H=0d0
  B=0d0
  do p=1,nAct
    H(p,p)=eHF(nC+p)
  end do

! The sqrt(2) is required by the QuAcK restricted-GW Lehmann residues.
  mu=0
  do i=nC+1,nO
    do m=1,nS
      mu=mu+1
      D(mu)=eHF(i)-Om(m)
      do p=1,nAct
        B(p,mu)=sqrt(2d0)*rho(nC+p,i,m)
      end do
    end do
  end do
  do a=nO+1,nOrb-nR
    do m=1,nS
      mu=mu+1
      D(mu)=eHF(a)+Om(m)
      do p=1,nAct
        B(p,mu)=sqrt(2d0)*rho(nC+p,a,m)
      end do
    end do
  end do

  H(1:nAct,nAct+1:nH)=B
  H(nAct+1:nH,1:nAct)=transpose(B)
  do mu=1,nAux
    H(nAct+mu,nAct+mu)=D(mu)
  end do

! QuAcK's symmetric diagonalizer overwrites H with normalized eigenvectors.
  call diagonalize_matrix(nH,H,E)

! Greedy one-to-one orbital/root assignment. This selects the connected
! quasiparticle branch when the 1p spectral weights are dominant.
! Near strong mixing / intruder states the resulting selection is not unique;
! the condition-number and residual checks below detect part of this problem.
  usedState=.false.
  usedOrb=.false.
  picked=0
  do k=1,nAct
    best=-1d0
    ibest=0
    pbest=0
    do p=1,nAct
      if(usedOrb(p)) cycle
      do j=1,nH
        if(usedState(j)) cycle
        val=H(p,j)**2
        if(val > best) then
          best=val
          ibest=j
          pbest=p
        end if
      end do
    end do
    if(ibest == 0) then
      info=3
      return
    end if
    picked(pbest)=ibest
    usedState(ibest)=.true.
    usedOrb(pbest)=.true.
  end do

  allocate(RP(nAct,nAct),RQ(nAux,nAct),RPwork(nAct,nAct))
  allocate(sig(nAct),work(max(1,5*nAct)))
  do p=1,nAct
    RP(:,p)=H(1:nAct,picked(p))
    RQ(:,p)=H(nAct+1:nH,picked(p))
    eCOT(nC+p)=E(picked(p))
  end do

! Check rank / conditioning of the projection before constructing S.
  RPwork=RP
  lwork=size(work)
  call dgesvd('N','N',nAct,nAct,RPwork,nAct,sig,dummyU,1,dummyVT,1,work,lwork,istat)
  if(istat /= 0) then
    info=3
    return
  end if
  if(sig(nAct) < 1d-10) then
    info=4
    return
  end if
  condition=sig(1)/sig(nAct)
  if(condition > condLimit) then
    info=4
    return
  end if

! Solve RP^T S^T = RQ^T, without explicitly computing RP^{-1}.
  allocate(rhs(nAct,nAux),ipiv(nAct),S(nAux,nAct))
  RPwork=transpose(RP)
  rhs=transpose(RQ)
  call dgesv(nAct,nAux,RPwork,nAct,ipiv,rhs,nAct,istat)
  if(istat /= 0) then
    info=3
    return
  end if
  S=transpose(rhs)

  allocate(hBlock(nAct,nAct),bs(nAct,nAct),N(nAct,nAct),eN(nAct))
  allocate(ds(nAux,nAct),hProj(nAct,nAct),Nhalf(nAct,nAct))
  bs=matmul(B,S)
  hBlock=bs
  do p=1,nAct
    hBlock(p,p)=hBlock(p,p)+eHF(nC+p)
  end do
  hCOT(nC+1:nOrb-nR,nC+1:nOrb-nR)=hBlock

! Invariant-subspace residual: B^T + D S - S hCOT = 0.
  allocate(rTest(nAux,nAct))
  rTest=transpose(B)-matmul(S,hBlock)
  do mu=1,nAux
    rTest(mu,:)=rTest(mu,:)+D(mu)*S(mu,:)
  end do
  resid=maxval(abs(rTest))

! Hermitian model-space Hamiltonian in an orthonormal metric:
! Omega^T H Omega = A+BS+(BS)^T+S^T D S.
  N=matmul(transpose(S),S)
  do p=1,nAct
    N(p,p)=N(p,p)+1d0
  end do
  ds=S
  do mu=1,nAux
    ds(mu,:)=D(mu)*S(mu,:)
  end do
  hProj=bs+transpose(bs)+matmul(transpose(S),ds)
  do p=1,nAct
    hProj(p,p)=hProj(p,p)+eHF(nC+p)
  end do

  call diagonalize_matrix(nAct,N,eN)
  if(minval(eN) <= 0d0) then
    info=3
    return
  end if
  allocate(NormVec(nAct,nAct))
  NormVec=N
  do k=1,nAct
    NormVec(:,k)=NormVec(:,k)/sqrt(eN(k))
  end do
  Nhalf=matmul(NormVec,transpose(N))
  hBlock=matmul(Nhalf,matmul(hProj,Nhalf))
  hBlock=0.5d0*(hBlock+transpose(hBlock))
  hHerm(nC+1:nOrb-nR,nC+1:nOrb-nR)=hBlock

! Confirm that hHerm has exactly the selected GW eigenvalues.
  allocate(eigWork(nAct,nAct),eCheck(nAct),eSelected(nAct))
  eigWork=hBlock
  call diagonalize_matrix(nAct,eigWork,eCheck)
  eSelected=eCOT(nC+1:nOrb-nR)
  do p=1,nAct-1
    do q=p+1,nAct
      if(eSelected(q) < eSelected(p)) then
        swap=eSelected(p)
        eSelected(p)=eSelected(q)
        eSelected(q)=swap
      end if
    end do
  end do
  ediff=maxval(abs(eCheck-eSelected))

  write(*,*)
  write(*,'(A)') '-----------------------------------------------'
  write(*,'(A)') '     G0W0-based correlated orbital theory     '
  write(*,'(A)') '-----------------------------------------------'
  write(*,'(A,I0)') 'Auxiliary dimension      : ',nAux
  write(*,'(A,ES12.3)') 'One-particle condition # : ',condition
  write(*,'(A,ES12.3)') 'Bloch residual (Eh)      : ',resid
  write(*,'(A,ES12.3)') 'Hermitian eigenvalue err : ',ediff
  write(*,'(A)') 'Orb    HF (eV)      COT (eV)      1p weight'
  do p=1,nAct
    write(*,'(I3,2F14.6,F14.6)') nC+p,eHF(nC+p)*HaToeV, &
             eCOT(nC+p)*HaToeV,sum(H(1:nAct,picked(p))**2)
  end do
  write(*,*)

  if(resid > residualTol .or. ediff > residualTol) then
    write(*,*) 'GW-COT: warning, block-decoupling validation failed.'
    info=5
  end if

end subroutine RGW_COT
