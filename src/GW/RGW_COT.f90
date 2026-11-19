subroutine RGW_COT(TDA_COT,useDavidson,davidsonThreshold,nOrb,nC,nO,nV,nR,nS, &
                   ENuc,T,V,eHF,ERI,hCOT,hHerm,eCOT,totalEnergy,info)

! One-shot GW-based correlated-orbital theory (COT), restricted closed shell.
!
! Build the reference screening and solve the Hermitian G0W0 eigenproblem.
! Select one quasiparticle root per active orbital and construct the Bloch
! effective one-particle Hamiltonian.
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
! Evaluate the restricted mean-field energy of the determinant formed from
! the occupied eigenvectors of hHerm. T, V and ERI are in the RHF MO basis.
! Evaluate full direct-RPA correlation energies for both orbital references.
!
! The active orbitals are nC+1 : nOrb-nR. Frozen orbitals remain diagonal.
! Dense diagonalization scales as O((nAct*(nS+1))^3) and should be used
! only as a numerical reference for modest systems.
!
! info = 0 success
!        1 invalid dimensions
!        2 supermatrix dimension exceeds the applicable limit
!        3 eigensolver/SVD/linear solve failed
!        4 selected one-particle projection ill-conditioned
!        5 invariance/eigenvalue residual failed
!        6 screening modes are nonfinite or have nonpositive energies
!        7 direct-RPA correlation energy could not be evaluated
!        8 Davidson root following did not converge

  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  include 'parameters.h'

  logical,intent(in)            :: TDA_COT,useDavidson
  integer,intent(in)            :: nOrb,nC,nO,nV,nR,nS
  double precision,intent(in)   :: davidsonThreshold,ENuc,T(nOrb,nOrb),V(nOrb,nOrb)
  double precision,intent(in)   :: eHF(nOrb),ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(out)  :: hCOT(nOrb,nOrb),hHerm(nOrb,nOrb)
  double precision,intent(out)  :: eCOT(nOrb),totalEnergy
  integer,intent(out)           :: info

  integer,parameter             :: maxDense = 4000
  double precision,parameter    :: condLimit = 1d10
  double precision,parameter    :: residualTol = 1d-7

  integer                       :: nAct,nOcc,nVir,nAux,nH
  integer(int64)                :: supermatrixDimension
  integer                       :: p,q,k,mu,i,a,m,istat,lwork
  integer,allocatable           :: ipiv(:)
  double precision              :: condition,resid,ediff,swap,blochTol,eigenvalueTol
  double precision              :: EcRHF,EcCOT
  double precision,allocatable  :: Aph(:,:),Bph(:,:),Om(:),XpY(:,:),XmY(:,:),rho(:,:,:)
  double precision,allocatable  :: B(:,:),D(:),roots(:,:)
  double precision,allocatable  :: RP(:,:),RQ(:,:),RPwork(:,:),rhs(:,:),S(:,:)
  double precision,allocatable  :: sig(:),work(:),N(:,:),eN(:),Nhalf(:,:)
  double precision,allocatable  :: NormVec(:,:),hProj(:,:),hBlock(:,:),bs(:,:),ds(:,:)
  double precision,allocatable  :: eCheck(:),eSelected(:),eigWork(:,:),rTest(:,:)
  double precision,allocatable  :: density(:,:)
  double precision,allocatable  :: ERI_COT(:,:,:,:)
  double precision              :: dummyU(1,1),dummyVT(1,1)

  info=0
  hCOT=0d0
  hHerm=0d0
  eCOT=eHF
  totalEnergy=0d0
  do p=1,nOrb
    hCOT(p,p)=eHF(p)
    hHerm(p,p)=eHF(p)
  end do

  nAct=nOrb-nC-nR
  nOcc=nO-nC
  nVir=nV-nR
  if(nC < 0 .or. nR < 0 .or. nAct < 1 .or. nOcc < 0 .or. nVir < 0 .or. &
     nV /= nOrb-nO .or. nS /= int(nOcc,int64)*nVir) then
    info=1
    return
  end if

! Initialize the spin-summed density with the reference occupations.
  allocate(density(nOrb,nOrb))
  density=0d0
  do i=1,nO
    density(i,i)=2d0
  end do

  if(nS == 0) then
    write(*,'(A)') 'GW-COT: empty screening space. The COT Hamiltonians reduce to the RHF reference.'
    call compute_total_energy()
    call print_correlation_energies(0d0,0d0)
    return
  end if

  supermatrixDimension=int(nAct,int64)*(int(nS,int64)+1)
  if(supermatrixDimension > huge(nH)) then
    write(*,'(A,I0)') 'GW-COT: supermatrix dimension exceeds the integer range: ',supermatrixDimension
    info=2
    return
  end if
  if(.not.useDavidson .and. supermatrixDimension > maxDense) then
    write(*,'(A,I0,A,I0)') 'GW-COT: dense dimension ',supermatrixDimension,' exceeds limit ',maxDense
    info=2
    return
  end if
  nH=int(supermatrixDimension)
  nAux=nH-nAct

! Singlet direct RPA screening evaluated at the RHF reference.
  allocate(Aph(nS,nS),Bph(nS,nS),Om(nS),XpY(nS,nS),XmY(nS,nS))
  call phRLR_A(1,.true.,nOrb,nC,nO,nV,nR,nS,1d0,eHF,ERI,Aph)
  call phRLR_B(1,.true.,nOrb,nC,nO,nV,nR,nS,1d0,ERI,Bph)
  call phRLR(TDA_COT,nS,Aph,Bph,EcRHF,Om,XpY,XmY)
  if(.not.all(ieee_is_finite(Om)) .or. any(Om <= 0d0) .or. .not.all(ieee_is_finite(XpY))) then
    write(*,'(A)') 'GW-COT: invalid screening modes at the RHF reference.'
    info=6
    return
  end if
  if(TDA_COT) then
    call compute_drpa_energy(eHF(nC+1:nOrb-nR),Bph,EcRHF)
    if(info /= 0) return
  end if
  allocate(rho(nOrb,nOrb,nS))
  call RGW_excitation_density(nOrb,nC,nO,nR,nS,ERI,XpY,rho)
  deallocate(Aph,XpY,XmY)

  allocate(B(nAct,nAux),D(nAux))

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
  deallocate(Om,rho)

  allocate(roots(nH,nAct))
  call RGW_COT_roots(useDavidson,davidsonThreshold,nAct,nAux,eHF(nC+1:nOrb-nR),B,D, &
                     eCOT(nC+1:nOrb-nR),roots,info)
  if(info /= 0) return

  allocate(RP(nAct,nAct),RQ(nAux,nAct),RPwork(nAct,nAct))
  allocate(sig(nAct),work(max(1,5*nAct)))
  RP=roots(1:nAct,:)
  RQ=roots(nAct+1:nH,:)
  deallocate(roots)

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
             eCOT(nC+p)*HaToeV,sum(RP(:,p)**2)
  end do
  write(*,*)

  blochTol=residualTol
  eigenvalueTol=residualTol
  if(useDavidson) then
    ! The residual matrix has nAct columns bounded by the Davidson threshold.
    ! The Bloch residual transforms it from both sides through the projection.
    eigenvalueTol=eigenvalueTol+sqrt(dble(nAct))*davidsonThreshold
    blochTol=blochTol+sqrt(dble(nAct))*davidsonThreshold/sig(nAct)**2
  end if
  if(resid > blochTol .or. ediff > eigenvalueTol) then
    write(*,*) 'GW-COT: warning, block-decoupling validation failed.'
    info=5
    return
  end if

! Fill the lowest nOcc eigenvectors of the Hermitian active block. The frozen
! core stays doubly occupied and the excluded virtual orbitals stay empty.
  density(nC+1:nOrb-nR,nC+1:nOrb-nR)=2d0*matmul(eigWork(:,1:nOcc),transpose(eigWork(:,1:nOcc)))
  call compute_total_energy()

! Rotate the active integrals with the Hermitian COT eigenvectors, whose
! columns follow the same energy ordering as eCheck and the density above.
  allocate(ERI_COT(nAct,nAct,nAct,nAct))
  call AOtoMO_ERI_RHF(nAct,nAct,eigWork, &
    ERI(nC+1:nOrb-nR,nC+1:nOrb-nR,nC+1:nOrb-nR,nC+1:nOrb-nR),ERI_COT)
  call phRLR_B(1,.true.,nAct,0,nOcc,nVir,0,nS,1d0,ERI_COT,Bph)
  deallocate(ERI_COT)
  call compute_drpa_energy(eCheck,Bph,EcCOT)
  if(info /= 0) return
  call print_correlation_energies(EcRHF,EcCOT)

contains

  subroutine compute_total_energy()

    double precision :: ET,EV,EJ,EK,electronicEnergy
    double precision,allocatable :: J(:,:),K(:,:)

    allocate(J(nOrb,nOrb),K(nOrb,nOrb))
    call Hartree_matrix_AO_basis(nOrb,density,ERI,J)
    call exchange_matrix_AO_basis(nOrb,density,ERI,K)

    ! For the symmetric density, sum(P*M) equals Tr(P M). K includes its sign.
    ET=sum(density*T)
    EV=sum(density*V)
    EJ=0.5d0*sum(density*J)
    EK=0.25d0*sum(density*K)
    electronicEnergy=ET+EV+EJ+EK
    totalEnergy=ENuc+electronicEnergy

    write(*,'(A)') 'COT mean-field energy (Eh)'
    write(*,'(A30,F21.12)') 'Kinetic energy = ',ET
    write(*,'(A30,F21.12)') 'Electron-nuclear energy = ',EV
    write(*,'(A30,F21.12)') 'Hartree energy = ',EJ
    write(*,'(A30,F21.12)') 'Exchange energy = ',EK
    write(*,'(A30,F21.12)') 'COT electronic energy = ',electronicEnergy
    write(*,'(A30,F21.12)') 'Nuclear repulsion = ',ENuc
    write(*,'(A30,F21.12)') 'COT total energy = ',totalEnergy
    write(*,*)

  end subroutine compute_total_energy

  subroutine compute_drpa_energy(orbitalEnergies,responseB,correlationEnergy)

    double precision,intent(in)  :: orbitalEnergies(nAct),responseB(nS,nS)
    double precision,intent(out) :: correlationEnergy

    integer :: occ,vir,ia
    double precision,external :: trace_matrix
    double precision :: gaps(nS),sqrtGaps(nS),omegaSquared(nS),traceA,tolerance
    double precision,allocatable :: casida(:,:)

    do occ=1,nOcc
      do vir=1,nVir
        ia=(occ-1)*nVir+vir
        gaps(ia)=orbitalEnergies(nOcc+vir)-orbitalEnergies(occ)
      end do
    end do
    if(any(gaps < 0d0)) then
      write(*,'(A)') 'GW-COT: negative occupied-virtual gap in the direct-RPA energy calculation.'
      info=7
      return
    end if

    ! For direct RPA, A-B is the diagonal gap matrix. Use only the squared
    ! frequencies, so a zero gap needs no inverse square root.
    traceA=sum(gaps)+trace_matrix(nS,responseB)
    sqrtGaps=sqrt(gaps)
    allocate(casida(nS,nS))
    casida=2d0*responseB
    do ia=1,nS
      casida(ia,ia)=casida(ia,ia)+gaps(ia)
      casida(:,ia)=sqrtGaps*casida(:,ia)*sqrtGaps(ia)
    end do
    call diagonalize_matrix(nS,casida,omegaSquared)
    tolerance=100d0*epsilon(1d0)*maxval(abs(omegaSquared))
    if(.not.all(ieee_is_finite(omegaSquared)) .or. any(omegaSquared < -tolerance)) then
      write(*,'(A)') 'GW-COT: invalid squared frequencies in the direct-RPA energy calculation.'
      info=7
      return
    end if
    correlationEnergy=0.5d0*(sum(sqrt(max(omegaSquared,0d0)))-traceA)

  end subroutine compute_drpa_energy

  subroutine print_correlation_energies(referenceEnergy,cotEnergy)

    double precision,intent(in) :: referenceEnergy,cotEnergy

    write(*,'(A32,F21.12,A)') 'dRPA@RHF correlation energy = ',referenceEnergy,' Eh'
    write(*,'(A32,F21.12,A)') 'dRPA@COT correlation energy = ',cotEnergy,' Eh'
    write(*,*)

  end subroutine print_correlation_energies

end subroutine RGW_COT
