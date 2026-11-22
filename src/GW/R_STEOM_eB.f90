subroutine R_STEOM_eB(TDA_W,TDA,singlet,triplet,Davidson,nRoots,tolerance,nOrb,nC,nO,nV,nR,nS,eHF,ERI, &
                      model,chargedDavidson,chargedTolerance)

! Static second-order electron-boson STEOM in the canonical RHF basis.
! Singlet dRPA modes supply the screening for both neutral spin manifolds.
! Build the bare electron-hole kernel from the residual two-electron interaction.

  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  include 'parameters.h'

  logical,intent(in) :: TDA_W,TDA,singlet,triplet,Davidson,chargedDavidson
  character(*),intent(in) :: model
  integer,intent(in) :: nRoots,nOrb,nC,nO,nV,nR,nS
  double precision,intent(in) :: tolerance,chargedTolerance,eHF(nOrb),ERI(nOrb,nOrb,nOrb,nOrb)

  integer :: nAct,nOcc,i,j,a,b,ia,jb,p,q,mode
  double precision :: ec,omega,delta,plus,minus,poleTolerance,denominator
  double precision,allocatable :: aph(:,:),bph(:,:),singletShift(:,:),om(:),xpY(:,:),xmY(:,:),rho(:,:,:)
  double precision,allocatable :: energies(:),coupling(:,:),weighted(:,:),vertex(:,:),selfEnergy(:,:)

  select case(model)
  case('osCOT')
    write(*,'(/,A)') 'STEOM-eB (joint osCOT amplitudes, quadratic transformation, forward singles)'
    write(*,'(A,L1)') '  TDA screening: ',TDA_W
  case('COT2')
    write(*,'(/,A)') 'STEOM-eB (second-order osCOT, canonical RHF, forward singles)'
    write(*,'(A,L1)') '  TDA screening: ',TDA_W
  case default
    write(*,'(/,A)') 'STEOM-eB (legacy response, canonical RHF)'
    write(*,'(A,L1,A,L1)') '  TDA response: ',TDA,'   TDA screening: ',TDA_W
  end select
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
  if(model /= 'legacy') then
    if(model == 'osCOT' .and. singlet) then
      call phRLR_B(1,.false.,nOrb,nC,nO,nV,nR,nS,1d0,ERI,bph)
    else
      deallocate(bph)
    end if
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

  subroutine joint_COT_amplitudes(graph)
    double precision,allocatable,intent(out) :: graph(:,:)
    double precision,allocatable :: fock(:,:),B(:,:),D(:),poles(:),weights(:)
    double precision :: condition,residual,crossFromHoles,crossFromParticles
    integer :: nAux,info,q,first,last
    integer(int64) :: nH

    nH=int(nAct,int64)*(int(nS,int64)+1)
    if(nH > huge(nAux)) error stop 'STEOM-eB joint GW supermatrix exceeds the integer range'
    if(.not.chargedDavidson .and. nH > 4000) error stop 'Enable COT Davidson for joint GW dimensions above 4000'
    nAux=int(nH)-nAct
    allocate(fock(nAct,nAct),B(nAct,nAux),D(nAux),graph(nAux,nAct),poles(nAct),weights(nAct))
    fock=0d0
    do q=1,nAct
      fock(q,q)=energies(q)
    end do
    call RGW_COT_auxiliaries(nOrb,nAct,nC,nOcc,nS,eHF,om,rho,B,D)
    call RGW_COT_wave_operator(chargedDavidson,chargedTolerance,nAct,nAux,fock,B,D,energies,0d0, &
                              poles,graph,weights,condition,residual,info)
    if(info /= 0) then
      write(*,'(A,I0)') 'STEOM-eB joint osCOT wave operator failed with info = ',info
      error stop 'Cannot construct joint osCOT amplitudes for STEOM-eB'
    end if
    first=nOcc*nS+1
    last=nOcc*nS
    crossFromHoles=norm2(graph(first:nAux,1:nOcc))
    crossFromParticles=norm2(graph(1:last,nOcc+1:nAct))
    write(*,'(A,ES12.4)') '  Joint osCOT Bloch residual (Eh) = ',residual
    write(*,'(A,2ES12.4)') '  Joint osCOT cross-block norms (hole, particle) = ',crossFromHoles,crossFromParticles
  end subroutine joint_COT_amplitudes

  subroutine assemble_COT_singles()
    double precision,allocatable :: Vmode(:,:),creation(:,:),annihilation(:,:),commutator(:,:),dressing(:,:),graph(:,:)
    double precision,allocatable :: rightMixed(:),leftMixed(:),bareCoupling(:)
    double precision :: forwardNorm,reverseNorm
    integer :: nVir,mu,nMixed

    nVir=nAct-nOcc
    allocate(Vmode(nAct,nAct),creation(nAct,nAct),dressing(nAct,nAct))
    allocate(annihilation(nAct,nAct),commutator(nAct,nAct))
    if(model == 'osCOT') then
      call joint_COT_amplitudes(graph)
    else
      allocate(graph(0,0))
    end if
    nMixed=merge(nS,0,model == 'osCOT' .and. singlet)
    allocate(rightMixed(nMixed),leftMixed(nMixed),bareCoupling(nMixed))
    dressing=0d0
    forwardNorm=0d0
    reverseNorm=0d0
    do mode=1,nS
      Vmode=sqrt(2d0)*rho(nC+1:nOrb-nR,nC+1:nOrb-nR,mode)
      if(.not.all(ieee_is_finite(Vmode))) error stop 'Nonfinite STEOM-eB mode couplings'

      creation=0d0
      if(model == 'osCOT') then
        ! Embed all four graph blocks as the coefficients of a one-boson
        ! generator. Leave the unprescribed neutral conversion blocks zero.
        annihilation=0d0
        do q=1,nOcc
          mu=(q-1)*nS+mode
          creation(1:nOcc,q)=-graph(mu,1:nOcc)
          annihilation(q,:)=graph(mu,:)
        end do
        do p=nOcc+1,nAct
          mu=(p-1)*nS+mode
          creation(p,:)=graph(mu,:)
          annihilation(nOcc+1:nAct,p)=-graph(mu,nOcc+1:nAct)
        end do
        do q=1,nAct
          commutator(:,q)=Vmode(:,q)+(om(mode)+energies-energies(q))*creation(:,q)
        end do
        ! The zero-boson correction is E(V) E(A) - E(B) E(C),
        ! with A and B the creation and annihilation amplitudes and
        ! C = V + [F,A] + Omega*A. Keep terms through degree two in V and K.
        call add_one_body_product(annihilation,commutator,-1d0,nVir,dressing)
        if(singlet) then
          do i=1,nOcc
            do a=nOcc+1,nAct
              ia=(i-1)*nVir+a-nOcc
              rightMixed(ia)=creation(a,i)
              leftMixed(ia)=annihilation(i,a)
              bareCoupling(ia)=Vmode(i,a)
            end do
          end do
          ! [W_N,K] contracts the mixed amplitudes with the bare singlet
          ! TDHF B block. The pure bosons are outside the final singles space.
          forwardNorm=hypot(forwardNorm,norm2(bareCoupling+matmul(bph,rightMixed)))
          reverseNorm=hypot(reverseNorm,norm2(bareCoupling-matmul(bph,leftMixed)))
        end if
      else
        do q=1,nAct
          do p=1,nAct
            if(p <= nOcc .and. q > nOcc) cycle
            if(Vmode(p,q) == 0d0) cycle
            delta=energies(p)-energies(q)
            denominator=om(mode)+delta
            poleTolerance=sqrt(epsilon(1d0))*max(1d0,abs(delta),om(mode))
            if(abs(denominator) <= poleTolerance) call resonance('osCOT amplitude',p,q,mode,denominator)
            creation(p,q)=-Vmode(p,q)/denominator
          end do
        end do
      end if
      call add_one_body_product(Vmode,creation,1d0,nVir,dressing)
    end do

    if(nMixed > 0) write(*,'(A,2ES12.4)') &
      '  Singlet singles/pure-boson coupling norms (forward, reverse, Eh) = ', &
      sqrt(2d0)*forwardNorm,sqrt(2d0)*reverseNorm

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

  subroutine add_one_body_product(X,Y,scale,nVir,dressing)
    integer,intent(in) :: nVir
    double precision,intent(in) :: X(nAct,nAct),Y(nAct,nAct),scale
    double precision,intent(inout) :: dressing(nAct,nAct)
    integer :: jb,j,b,ia,i,a

    ! Reference-subtracted singles matrix of the product of two normal-
    ! ordered one-body operators. The contraction is X Pv Y - Y Po X.
    call dgemm('N','N',nAct,nAct,nVir,scale,X(1,nOcc+1),nAct,Y(nOcc+1,1),nAct,1d0,dressing,nAct)
    call dgemm('N','N',nAct,nAct,nOcc,-scale,Y,nAct,X,nAct,1d0,dressing,nAct)
    !$omp parallel do default(none) shared(nS,nOcc,nVir,aph,singlet,singletShift,X,Y,scale) &
    !$omp private(jb,j,b,ia,i,a)
    do jb=1,nS
      j=(jb-1)/nVir+1
      b=mod(jb-1,nVir)+nOcc+1
      do ia=1,nS
        i=(ia-1)/nVir+1
        a=mod(ia-1,nVir)+nOcc+1
        aph(ia,jb)=aph(ia,jb)-scale*(X(a,b)*Y(j,i)+Y(a,b)*X(j,i))
        if(singlet) singletShift(ia,jb)=singletShift(ia,jb)+2d0*scale*(X(a,i)*Y(j,b)+Y(a,i)*X(j,b))
      end do
    end do
    !$omp end parallel do
  end subroutine add_one_body_product

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
