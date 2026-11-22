subroutine R_STEOM_eB(TDA_W,TDA,singlet,triplet,Davidson,nRoots,tolerance,nOrb,nC,nO,nV,nR,nS,eHF,ERI,COT2)

! Static second-order electron-boson STEOM in the canonical RHF basis.
! Singlet dRPA modes supply the screening for both neutral spin manifolds.
! Build the bare electron-hole kernel from the residual two-electron interaction.

  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  include 'parameters.h'

  logical,intent(in) :: TDA_W,TDA,singlet,triplet,Davidson,COT2
  integer,intent(in) :: nRoots,nOrb,nC,nO,nV,nR,nS
  double precision,intent(in) :: tolerance,eHF(nOrb),ERI(nOrb,nOrb,nOrb,nOrb)

  integer :: nAct,nOcc,i,j,a,b,ia,jb,p,q,mode
  double precision :: ec,omega,delta,plus,minus,poleTolerance,denominator
  double precision,allocatable :: aph(:,:),bph(:,:),singletShift(:,:),om(:),xpY(:,:),xmY(:,:),rho(:,:,:)
  double precision,allocatable :: energies(:),coupling(:,:),weighted(:,:),vertex(:,:),selfEnergy(:,:)

  if(COT2) then
    write(*,'(/,A)') 'STEOM-eB (second-order osCOT, canonical RHF, forward singles)'
    write(*,'(A,L1)') '  TDA screening: ',TDA_W
  else
    write(*,'(/,A)') 'STEOM-eB (legacy response, canonical RHF)'
    write(*,'(A,L1,A,L1)') '  TDA response: ',TDA,'   TDA screening: ',TDA_W
  end if
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
  ! Start from the triplet CIS matrix. The bare singlet shift is 2(ia|jb).
  call phRLR_A(2,.false.,nOrb,nC,nO,nV,nR,nS,1d0,eHF,ERI,aph)
  allocate(singletShift(merge(nS,0,singlet),merge(nS,0,singlet)))
  if(singlet) call phRLR_B(1,.true.,nOrb,nC,nO,nV,nR,nS,1d0,ERI,singletShift)
  if(COT2) then
    deallocate(bph)
    call assemble_COT_singles()
    deallocate(rho)
    if(triplet) call print_COT_spectrum('triplet')
    if(singlet) then
      aph=aph+singletShift
      call print_COT_spectrum('singlet')
    end if
    return
  end if

  if(TDA) then
    deallocate(bph)
  else
    call phRLR_B(2,.false.,nOrb,nC,nO,nV,nR,nS,1d0,ERI,bph)
  end if

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

  subroutine assemble_COT_singles()
    double precision,allocatable :: Vmode(:,:),amplitude(:,:),dressing(:,:)
    integer :: nVir

    nVir=nAct-nOcc
    allocate(Vmode(nAct,nAct),amplitude(nAct,nAct),dressing(nAct,nAct))
    dressing=0d0
    do mode=1,nS
      Vmode=sqrt(2d0)*rho(nC+1:nOrb-nR,nC+1:nOrb-nR,mode)
      if(.not.all(ieee_is_finite(Vmode))) error stop 'Nonfinite STEOM-eB mode couplings'

      ! Boson-creation amplitude L. The occupied/virtual block is fixed to zero.
      amplitude=0d0
      do q=1,nAct
        do p=1,nAct
          if(p <= nOcc .and. q > nOcc) cycle
          if(Vmode(p,q) == 0d0) cycle
          delta=energies(p)-energies(q)
          denominator=om(mode)+delta
          poleTolerance=sqrt(epsilon(1d0))*max(1d0,abs(delta),om(mode))
          if(abs(denominator) <= poleTolerance) call resonance('osCOT amplitude',p,q,mode,denominator)
          amplitude(p,q)=-Vmode(p,q)/denominator
        end do
      end do

      ! D = sum_modes (V P_virtual L - L P_occupied V).
      call dgemm('N','N',nAct,nAct,nVir,1d0,Vmode(1,nOcc+1),nAct,amplitude(nOcc+1,1),nAct,1d0,dressing,nAct)
      call dgemm('N','N',nAct,nAct,nOcc,-1d0,amplitude,nAct,Vmode,nAct,1d0,dressing,nAct)

      !$omp parallel do default(none) shared(nS,nOcc,nVir,aph,singlet,singletShift,Vmode,amplitude) &
      !$omp private(jb,j,b,ia,i,a)
      do jb=1,nS
        j=(jb-1)/nVir+1
        b=mod(jb-1,nVir)+nOcc+1
        do ia=1,nS
          i=(ia-1)/nVir+1
          a=mod(ia-1,nVir)+nOcc+1
          aph(ia,jb)=aph(ia,jb)-Vmode(a,b)*amplitude(j,i)-amplitude(a,b)*Vmode(j,i)
          if(singlet) singletShift(ia,jb)=singletShift(ia,jb)+2d0*amplitude(a,i)*Vmode(j,b)
        end do
      end do
      !$omp end parallel do
    end do

    do jb=1,nS
      j=(jb-1)/nVir+1
      b=mod(jb-1,nVir)+nOcc+1
      do ia=1,nS
        i=(ia-1)/nVir+1
        a=mod(ia-1,nVir)+nOcc+1
        if(i == j) aph(ia,jb)=aph(ia,jb)+dressing(a,b)
        if(a == b) aph(ia,jb)=aph(ia,jb)-dressing(j,i)
      end do
    end do
  end subroutine assemble_COT_singles

  subroutine print_COT_spectrum(manifold)
    character(*),intent(in) :: manifold
    double precision,allocatable :: matrix(:,:),roots(:),imaginary(:),work(:),diagonal(:)
    integer,allocatable :: order(:)
    double precision :: residual,leftDummy(1,1),rightDummy(1,1),workQuery(1),noFactor(nS,0)
    integer :: nReport,info,iterations,k

    if(.not.all(ieee_is_finite(aph))) error stop 'Nonfinite STEOM-eB singles matrix'
    nReport=nS
    if(Davidson) then
      nReport=min(nRoots,nS)
      allocate(roots(nReport),diagonal(nS))
      do k=1,nS
        diagonal(k)=aph(k,k)
      end do
      call STEOM_eB_davidson(.true.,.false.,nS,nReport,tolerance,aph,noFactor,diagonal,roots,iterations,residual,info)
      if(info == 0) then
        write(*,'(A,A,A,I0,A,I0,A,ES12.4)') 'STEOM-eB Davidson (',manifold,'): ',nReport, &
          ' roots, iterations = ',iterations,', max residual = ',residual
        call print_excitation_energies('STEOM-eB',manifold,nReport,roots)
        return
      end if
      write(*,'(A,A,A)') 'STEOM-eB Davidson (',manifold,'): switching to dense diagonalization'
      deallocate(roots)
    end if

    matrix=aph
    allocate(roots(nS),imaginary(nS),order(nS))
    call dgeev('N','N',nS,matrix,nS,roots,imaginary,leftDummy,1,rightDummy,1,workQuery,-1,info)
    if(info /= 0) error stop 'STEOM-eB nonsymmetric eigensolver workspace query failed'
    allocate(work(int(workQuery(1))))
    call dgeev('N','N',nS,matrix,nS,roots,imaginary,leftDummy,1,rightDummy,1,work,size(work),info)
    if(info /= 0 .or. .not.all(ieee_is_finite(roots)) .or. .not.all(ieee_is_finite(imaginary))) &
      error stop 'STEOM-eB nonsymmetric diagonalization failed'
    order=[(k,k=1,nS)]
    call quick_sort(roots,order,nS)
    imaginary=imaginary(order)

    if(any(abs(imaginary(1:nReport)) > 1d-8)) then
      write(*,'(/,A,A,A)') 'STEOM-eB calculation: ',manifold,' manifold (complex roots present)'
      write(*,'(A)') '   State           Re (Eh)           Im (Eh)           Re (eV)           Im (eV)'
      do k=1,min(50,nReport)
        write(*,'(I8,4(1X,F17.8))') k,roots(k),imaginary(k),roots(k)*HaToeV,imaginary(k)*HaToeV
      end do
    else
      call print_excitation_energies('STEOM-eB',manifold,nReport,roots)
    end if
  end subroutine print_COT_spectrum

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
      call STEOM_eB_davidson(TDA,.true.,nS,nReport,tolerance,matrix,factor,diagonal,roots,iterations,residual,info)
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
