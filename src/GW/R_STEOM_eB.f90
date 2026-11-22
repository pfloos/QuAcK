subroutine R_STEOM_eB(TDA_W,TDA,singlet,triplet,Davidson,nRoots,tolerance,nOrb,nC,nO,nV,nR,nS,eHF,ERI)

! Static second-order electron-boson STEOM in the canonical RHF basis.
! Singlet dRPA modes supply the screening for both neutral spin manifolds.
! The minimal eB Hamiltonian contributes no bare electron-hole kernel.

  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

  logical,intent(in) :: TDA_W,TDA,singlet,triplet,Davidson
  integer,intent(in) :: nRoots,nOrb,nC,nO,nV,nR,nS
  double precision,intent(in) :: tolerance,eHF(nOrb),ERI(nOrb,nOrb,nOrb,nOrb)

  integer :: nAct,nOcc,i,j,a,b,ia,jb,p,q,mode
  double precision :: ec,omega,delta,plus,minus,poleTolerance,denominator
  double precision,allocatable :: aph(:,:),bph(:,:),singletShift(:,:),om(:),xpY(:,:),xmY(:,:),rho(:,:,:)
  double precision,allocatable :: energies(:),coupling(:,:),weighted(:,:),vertex(:,:),selfEnergy(:,:)

  write(*,'(/,A)') 'STEOM-eB (second order, canonical RHF)'
  write(*,'(A,L1,A,L1)') '  TDA response: ',TDA,'   TDA screening: ',TDA_W
  if(nS == 0) then
    write(*,'(A)') '  No active electron-hole pairs.'
    return
  end if

  nAct=nOrb-nC-nR
  nOcc=nO-nC
  energies=eHF(nC+1:nOrb-nR)

  allocate(aph(nS,nS),bph(nS,nS),om(nS),xpY(nS,nS),xmY(nS,nS))
  call phRLR_A(1,.true.,nOrb,nC,nO,nV,nR,nS,1d0,eHF,ERI,aph)
  call phRLR_B(1,.true.,nOrb,nC,nO,nV,nR,nS,1d0,ERI,bph)
  call phRLR(TDA_W,nS,aph,bph,ec,om,xpY,xmY)
  if(.not.all(ieee_is_finite(om)) .or. any(om <= 0d0) .or. .not.all(ieee_is_finite(xpY))) &
    error stop 'STEOM-eB requires finite, positive screening modes'

  allocate(rho(nOrb,nOrb,nS))
  call RGW_excitation_density(nOrb,nC,nO,nR,nS,ERI,xpY,rho)
  deallocate(xpY,xmY)
  aph=0d0
  if(TDA) then
    deallocate(bph)
  else
    bph=0d0
  end if
  allocate(singletShift(merge(nS,0,singlet),merge(nS,0,singlet)),source=0d0)

  allocate(coupling(nAct,nAct),weighted(nAct,nAct),vertex(nAct,nAct),selfEnergy(nAct,nAct))
  selfEnergy=0d0
  do mode=1,nS
    omega=om(mode)
    coupling=sqrt(2d0)*rho(nC+1:nOrb-nR,nC+1:nOrb-nR,mode)
    if(.not.all(ieee_is_finite(coupling))) error stop 'Nonfinite STEOM-eB mode couplings'

    weighted=0d0
    vertex=0d0
    do q=1,nAct
      do p=1,nAct
        if(coupling(p,q) == 0d0) cycle
        delta=energies(p)-energies(q)
        plus=delta+omega
        minus=delta-omega
        poleTolerance=sqrt(epsilon(1d0))*max(1d0,abs(delta),omega)

        ! The one-body term is the symmetrized on-shell GW self-energy.
        denominator=merge(plus,minus,q <= nOcc)
        if(abs(denominator) <= poleTolerance) call resonance('GW self-energy',p,q,mode,denominator)
        vertex(p,q)=coupling(p,q)/denominator

        ! The screened exchange in A uses same-space pairs only.
        ! Cross-space pairs enter T in the singlet and Z in the full response.
        if(TDA .and. .not.singlet .and. ((p <= nOcc) .neqv. (q <= nOcc))) cycle
        if(min(abs(plus),abs(minus)) <= poleTolerance) &
          call resonance('neutral kernel',p,q,mode,merge(plus,minus,abs(plus) < abs(minus)))
        weighted(p,q)=coupling(p,q)*omega/(plus*minus)
      end do
    end do
    call dgemm('N','T',nAct,nAct,nAct,1d0,vertex,nAct,coupling,nAct,1d0,selfEnergy,nAct)

    jb=0
    do j=1,nOcc
      do b=nOcc+1,nAct
        jb=jb+1
        ia=0
        do i=1,nOcc
          do a=nOcc+1,nAct
            ia=ia+1

            ! Spin-independent contributions -U and -Z.
            aph(ia,jb)=aph(ia,jb)-coupling(a,b)*weighted(i,j)-weighted(a,b)*coupling(i,j)
            if(.not.TDA) bph(ia,jb)=bph(ia,jb)-coupling(a,j)*weighted(i,b)-weighted(a,j)*coupling(i,b)

            ! The singlet adds 2T to both response blocks.
            if(singlet) singletShift(ia,jb)=singletShift(ia,jb) &
              +2d0*(coupling(a,i)*weighted(b,j)+weighted(a,i)*coupling(b,j))
          end do
        end do
      end do
    end do
  end do
  deallocate(rho)

  selfEnergy=0.5d0*(selfEnergy+transpose(selfEnergy))
  jb=0
  do j=1,nOcc
    do b=nOcc+1,nAct
      jb=jb+1
      ia=0
      do i=1,nOcc
        do a=nOcc+1,nAct
          ia=ia+1
          if(i == j) aph(ia,jb)=aph(ia,jb)+selfEnergy(a,b)
          if(a == b) aph(ia,jb)=aph(ia,jb)-selfEnergy(j,i)
          if(ia == jb) aph(ia,jb)=aph(ia,jb)+energies(a)-energies(i)
        end do
      end do
    end do
  end do

  if(triplet) call print_spectrum('triplet')
  if(singlet) then
    aph=aph+singletShift
    if(.not.TDA) bph=bph+singletShift
    call print_spectrum('singlet')
  end if

contains

  subroutine resonance(term,p,q,mode,denominator)
    character(*),intent(in) :: term
    integer,intent(in) :: p,q,mode
    double precision,intent(in) :: denominator

    write(*,'(A,A,A,2(1X,I0),A,I0,A,ES13.5)') 'STEOM-eB ',term, &
      ' resonance for orbitals',nC+p,nC+q,', mode ',mode,', denominator (Eh) = ',denominator
    error stop 'Resonant electron-boson channel in the static STEOM-eB approximation'
  end subroutine resonance

  subroutine print_spectrum(manifold)
    character(*),intent(in) :: manifold
    double precision,allocatable :: matrix(:,:),factor(:,:),roots(:),work(:),diagonal(:)
    double precision :: residual
    integer :: lwork,info,nReport,iterations,ia

    nReport=nS
    if(Davidson) nReport=min(nRoots,nS)
    allocate(matrix(nS,nS),roots(nReport))
    allocate(factor(nS,merge(0,nS,TDA)))
    if(TDA) then
      matrix=aph
    else
      ! Stable real response: A-B = L L^T and omega^2 is an eigenvalue
      ! of the symmetric matrix L^T (A+B) L.
      factor=aph-bph
      call dpotrf('L',nS,factor,nS,info)
      if(info /= 0) then
        write(*,'(A)') 'STEOM-eB '//manifold//' response: A-B is not positive definite'
        error stop 'Unstable full STEOM-eB response'
      end if
      matrix=aph+bph
    end if

    if(Davidson) then
      allocate(diagonal(nS))
      do ia=1,nS
        diagonal(ia)=matrix(ia,ia)
        if(.not.TDA) diagonal(ia)=diagonal(ia)*factor(ia,ia)**2
      end do
      call STEOM_eB_davidson(TDA,nS,nReport,tolerance,matrix,factor,diagonal,roots,iterations,residual,info)
      write(*,'(A,A,A,I0,A,I0,A,ES12.4)') 'STEOM-eB Davidson (',manifold,'): ',nReport, &
        ' roots, iterations = ',iterations,', max residual = ',residual
      if(info /= 0) error stop 'STEOM-eB Davidson failed to converge'
    else
      if(.not.TDA) then
        call dtrmm('R','L','N','N',nS,nS,1d0,factor,nS,matrix,nS)
        call dtrmm('L','L','T','N',nS,nS,1d0,factor,nS,matrix,nS)
      end if
      lwork=max(1,3*nS-1)
      allocate(work(lwork))
      call dsyev('N','U',nS,matrix,nS,roots,work,lwork,info)
      if(info /= 0 .or. .not.all(ieee_is_finite(roots))) error stop 'STEOM-eB diagonalization failed'
    end if
    if(.not.TDA) then
      if(any(roots <= 0d0)) then
        write(*,'(A)') 'STEOM-eB '//manifold//' response: nonpositive squared excitation energy'
        error stop 'Unstable full STEOM-eB response'
      end if
      roots=sqrt(roots)
    end if
    call print_excitation_energies('STEOM-eB',manifold,nReport,roots)
  end subroutine print_spectrum

end subroutine R_STEOM_eB
