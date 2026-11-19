subroutine RGW_COT_roots(useDavidson,tolerance,nAct,nAux,f,B,D,energies,roots,info)

! Select one distinct GW pole for each active HF orbital, for the symmetric
! supermatrix [diag(f), B; transpose(B), diag(D)].

  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

  logical,intent(in) :: useDavidson
  integer,intent(in) :: nAct,nAux
  double precision,intent(in) :: tolerance,f(nAct),B(nAct,nAux),D(nAux)
  double precision,intent(out) :: energies(nAct),roots(nAct+nAux,nAct)
  integer,intent(out) :: info

  integer :: n

  n=nAct+nAux
  info=0
  if(useDavidson) then
    call davidson_roots()
  else
    call dense_roots()
  end if

contains

  subroutine dense_roots()

    integer :: i,picked(nAct)
    double precision,allocatable :: H(:,:),allEnergies(:),weights(:,:)

    write(*,'(A,I0)') 'GW-COT dense supermatrix dimension: ',n
    allocate(H(n,n),allEnergies(n),weights(n,nAct))
    H=0d0
    do i=1,nAct
      H(i,i)=f(i)
    end do
    do i=1,nAux
      H(nAct+i,nAct+i)=D(i)
    end do
    H(1:nAct,nAct+1:n)=B
    H(nAct+1:n,1:nAct)=transpose(B)
    call diagonalize_matrix(n,H,allEnergies)
    weights=transpose(H(1:nAct,:))**2
    if(.not.all(ieee_is_finite(weights)) .or. .not.all(ieee_is_finite(allEnergies))) then
      info=3
      return
    end if
    call assign_roots(weights,picked)
    roots=H(:,picked)
    energies=allEnergies(picked)

  end subroutine dense_roots

  subroutine davidson_roots()

    integer,parameter :: maxIterations=256
    double precision,parameter :: denominatorFloor=0.001d0
    integer :: limit,m,iteration,i,j,k,oldSize,nMatvec,nRestart,picked(nAct)
    double precision :: rnorm(nAct),maxResidual
    double precision,allocatable :: basis(:,:),images(:,:),projection(:,:),theta(:),weights(:,:)
    double precision,allocatable :: rootImages(:,:),corrections(:,:),trial(:),diagonal(:),denominator(:)

    ! Allow room for nearby satellites before restarting. In particular,
    ! high virtual roots can share their HF weight with several Ritz states.
    limit=min(n,max(128,10*nAct))
    allocate(basis(n,limit),images(n,limit),projection(limit,limit),theta(limit),weights(limit,nAct))
    allocate(rootImages(n,nAct),corrections(n,nAct),trial(n),diagonal(n),denominator(n))
    diagonal(1:nAct)=f
    diagonal(nAct+1:n)=D

    ! One unit vector per active occupied or virtual HF orbital.
    basis(:,1:nAct)=0d0
    images(1:nAct,1:nAct)=0d0
    do i=1,nAct
      basis(i,i)=1d0
      images(i,i)=f(i)
    end do
    images(nAct+1:n,1:nAct)=transpose(B)
    m=nAct
    nMatvec=nAct
    nRestart=0
    info=8
    maxResidual=huge(1d0)
    write(*,'(A,I0,A,I0)') 'GW-COT Davidson supermatrix dimension: ',n,'  max subspace: ',limit

    write(*,'(A,ES12.4)') 'GW-COT Davidson residual threshold (Eh): ',tolerance

    do iteration=1,maxIterations
      call dgemm('T','N',m,m,n,1d0,basis,n,images,n,0d0,projection,limit)
      projection(1:m,1:m)=0.5d0*(projection(1:m,1:m)+transpose(projection(1:m,1:m)))
      call diagonalize_matrix(m,projection(1:m,1:m),theta(1:m))
      if(.not.all(ieee_is_finite(theta(1:m)))) exit

      ! Use the physical components of the Ritz vectors for root following.
      ! Form full Ritz vectors only for the roots assigned to HF orbitals.
      weights(1:m,:)=matmul(transpose(projection(1:m,1:m)),transpose(basis(1:nAct,1:m)))**2
      if(.not.all(ieee_is_finite(weights(1:m,:)))) exit
      call assign_roots(weights(1:m,:),picked)
      do j=1,nAct
        k=picked(j)
        energies(j)=theta(k)
        roots(:,j)=matmul(basis(:,1:m),projection(1:m,k))
        rootImages(:,j)=matmul(images(:,1:m),projection(1:m,k))
        corrections(:,j)=rootImages(:,j)-energies(j)*roots(:,j)
        rnorm(j)=norm2(corrections(:,j))
      end do
      if(.not.all(ieee_is_finite(rnorm))) exit
      maxResidual=maxval(rnorm)
      if(maxResidual <= tolerance) then
        do j=1,nAct
          call apply(roots(:,j),trial)
          rootImages(:,j)=trial
          corrections(:,j)=trial-energies(j)*roots(:,j)
          rnorm(j)=norm2(corrections(:,j))
        end do
        nMatvec=nMatvec+nAct
        maxResidual=maxval(rnorm)
        if(.not.all(ieee_is_finite(rnorm))) exit
        if(maxResidual <= tolerance) then
          info=0
          exit
        end if
      end if

      ! Keep the followed Ritz vectors and restore the original HF target
      ! space at each restart so the orbital-character test stays anchored.
      if(m+nAct > limit .and. limit < n) then
        basis(:,1:nAct)=roots
        images(:,1:nAct)=rootImages
        m=nAct
        nRestart=nRestart+1
        do j=1,nAct
          trial=0d0
          trial(j)=1d0
          call append_vector(trial,basis,images,m,limit,nMatvec)
        end do
      end if

      oldSize=m
      do j=1,nAct
        if(rnorm(j) <= tolerance .or. m == limit) cycle
        denominator=energies(j)-diagonal
        denominator=sign(max(abs(denominator),denominatorFloor),denominator)
        trial=corrections(:,j)/denominator
        i=m
        call append_vector(trial,basis,images,m,limit,nMatvec)
        if(i /= m) cycle
        trial=corrections(:,j)
        call append_vector(trial,basis,images,m,limit,nMatvec)
      end do
      if(m == oldSize) exit
    end do

    write(*,'(A,I0,A,I0,A,I0)') 'GW-COT Davidson iterations: ',min(iteration,maxIterations), &
      '  matrix-vector products: ',nMatvec,'  restarts: ',nRestart
    write(*,'(A,ES12.4)') 'GW-COT largest eigenpair residual (Eh): ',maxResidual
    if(info /= 0) write(*,'(A)') 'GW-COT: Davidson root following did not converge.'

  end subroutine davidson_roots

  subroutine apply(x,y)

    double precision,intent(in) :: x(n)
    double precision,intent(out) :: y(n)

    y(1:nAct)=f*x(1:nAct)
    y(nAct+1:n)=D*x(nAct+1:n)
    call dgemv('N',nAct,nAux,1d0,B,nAct,x(nAct+1),1,1d0,y(1),1)
    call dgemv('T',nAct,nAux,1d0,B,nAct,x(1),1,1d0,y(nAct+1),1)

  end subroutine apply

  subroutine append_vector(trial,basis,images,m,limit,nMatvec)

    integer,intent(in) :: limit
    integer,intent(inout) :: m,nMatvec
    double precision,intent(inout) :: trial(n),basis(n,limit),images(n,limit)
    integer :: pass,j
    double precision :: norm

    do pass=1,2
      do j=1,m
        trial=trial-basis(:,j)*dot_product(basis(:,j),trial)
      end do
    end do
    norm=norm2(trial)
    if(norm <= 1d-13) return
    m=m+1
    basis(:,m)=trial/norm
    call apply(basis(:,m),images(:,m))
    nMatvec=nMatvec+1

  end subroutine append_vector

  subroutine assign_roots(weights,picked)

    double precision,intent(inout) :: weights(:,:)
    integer,intent(out) :: picked(nAct)
    integer :: j,location(2)

    ! Match the original dense selector's greedy one-to-one assignment,
    ! including the order of orbital and root indices used to break ties.
    do j=1,nAct
      location=maxloc(weights)
      picked(location(2))=location(1)
      weights(location(1),:)=-1d0
      weights(:,location(2))=-1d0
    end do

  end subroutine assign_roots

end subroutine RGW_COT_roots
