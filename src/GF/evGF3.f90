 subroutine evGF3(maxSCF,thresh,max_diis,renormalization,nBas,nC,nO,nV,nR,V,e0)

! Perform third-order Green function calculation in diagonal approximation

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: thresh
  integer,intent(in)            :: maxSCF,max_diis,renormalization
  integer,intent(in)            :: nBas,nC,nO,nV,nR
  double precision,intent(in)   :: e0(nBas),V(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: nSCF
  integer                       :: n_diis
  double precision              :: eps,eps1,eps2,Conv
  double precision              :: rcond
  double precision,allocatable  :: Sig2(:),SigInf(:),Sig3(:),eGF3(:),eOld(:)
  double precision,allocatable  :: App(:,:),Bpp(:,:),Cpp(:,:),Dpp(:,:)
  double precision,allocatable  :: Z(:),X2h1p(:),X1h2p(:),Sig2h1p(:),Sig1h2p(:)
  double precision,allocatable  :: error_diis(:,:),e_diis(:,:)

  integer                       :: i,j,k,l,a,b,c,d,p

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'|    Third-order Green function calculation    |'
  write(*,*)'************************************************'
  write(*,*)

! Memory allocation

  allocate(eGF3(nBas),eOld(nBas),                                       & 
           Sig2(nBas),SigInf(nBas),Sig3(nBas),                          &
           App(nBas,6),Bpp(nBas,2),Cpp(nBas,6),Dpp(nBas,6),             &
           Z(nBas),X2h1p(nBas),X1h2p(nBas),Sig2h1p(nBas),Sig1h2p(nBas), &
           error_diis(nBas,max_diis),e_diis(nBas,max_diis))

!------------------------------------------------------------------------
! Compute third-order frequency-independent contribution
!------------------------------------------------------------------------

  App(:,:) = 0d0

  do p=nC+1,nBas-nR
    do i=nC+1,nO
    do j=nC+1,nO
    do k=nC+1,nO
      do a=nO+1,nBas-nR
      do b=nO+1,nBas-nR

          eps1 = e0(j) + e0(i) - e0(a) - e0(b)
          eps2 = e0(k) + e0(i) - e0(a) - e0(b)

          App(p,1) = App(p,1) & 
                   - (2d0*V(p,k,p,j) - V(p,k,j,p))*(2d0*V(j,i,a,b) - V(j,i,b,a))*V(a,b,k,i)/(eps1*eps2)

      enddo
      enddo
    enddo
    enddo
    enddo
  enddo

  do p=nC+1,nBas-nR
    do i=nC+1,nO
    do j=nC+1,nO
      do a=nO+1,nBas-nR
      do b=nO+1,nBas-nR
      do c=nO+1,nBas-nR

          eps1 = e0(j) + e0(i) - e0(a) - e0(b)
          eps2 = e0(j) + e0(i) - e0(a) - e0(c)

          App(p,2) = App(p,2) & 
                   + (2d0*V(p,c,p,b) - V(p,c,b,p))*(2d0*V(j,i,a,b) - V(j,i,b,a))*V(j,i,c,a)/(eps1*eps2)

      enddo
      enddo
      enddo
    enddo
    enddo
  enddo
  
  do p=nC+1,nBas-nR
    do i=nC+1,nO
    do j=nC+1,nO
      do a=nO+1,nBas-nR
      do b=nO+1,nBas-nR
      do c=nO+1,nBas-nR

          eps1 = e0(j) + e0(i) - e0(a) - e0(b)
          eps2 = e0(j)         - e0(c)

          App(p,3) = App(p,3) & 
                   + (2d0*V(p,c,p,j) - V(p,c,j,p))*(2d0*V(j,i,a,b) - V(j,i,b,a))*V(a,b,c,i)/(eps1*eps2)

      enddo
      enddo
      enddo
    enddo
    enddo
  enddo

  App(:,4) = App(:,3)

  do p=nC+1,nBas-nR
    do i=nC+1,nO
    do j=nC+1,nO
    do k=nC+1,nO
      do a=nO+1,nBas-nR
      do b=nO+1,nBas-nR

          eps1 = e0(j) + e0(i) - e0(a) - e0(b)
          eps2 = e0(k)         - e0(b)

          App(p,5) = App(p,5) & 
                   - (2d0*V(p,b,p,k) - V(p,b,k,p))*(2d0*V(j,i,a,b) - V(j,i,b,a))*V(i,j,k,a)/(eps1*eps2)

      enddo
      enddo
    enddo
    enddo
    enddo
  enddo
  
  App(:,6) = App(:,5)

! Frequency-independent part of the third-order self-energy

  SigInf(:) = App(:,1) + App(:,2) + App(:,3) + App(:,4) + App(:,5) + App(:,6)

!------------------------------------------------------------------------
! Main SCF loop
!------------------------------------------------------------------------

  nSCF            = 0
  n_diis          = 0
  Conv            = 1d0
  Sig2(:)         = 0d0
  Sig3(:)         = 0d0
  e_diis(:,:)     = 0d0
  error_diis(:,:) = 0d0
  eGF3(:)         = e0(:)
  eOld(:)         = e0(:)

  do while(Conv > thresh .and. nSCF < maxSCF)

    ! Frequency-dependent second-order contribution

    Bpp(:,:) = 0d0

    do p=nC+1,nBas-nR
      do i=nC+1,nO
      do j=nC+1,nO
        do a=nO+1,nBas-nR

          eps = eGF3(p) + e0(a) - e0(i) - e0(j)

          Bpp(p,1) = Bpp(p,1) & 
                   + (2d0*V(p,a,i,j) - V(p,a,j,i))*V(p,a,i,j)/eps

        enddo
      enddo
      enddo
    enddo
 
    do p=nC+1,nBas-nR
      do i=nC+1,nO
        do a=nO+1,nBas-nR
        do b=nO+1,nBas-nR

          eps = eGF3(p) + e0(i) - e0(a) - e0(b)

          Bpp(p,2) = Bpp(p,2) &
                   + (2d0*V(p,i,a,b) - V(p,i,b,a))*V(p,i,a,b)/eps

        enddo
        enddo
      enddo
    enddo
 
    ! Total second-order Green function

    Sig2(:) = Bpp(:,1) + Bpp(:,2)

    ! Frequency-dependent third-order contribution: "C" terms

    Cpp(:,:) = 0d0

    do p=nC+1,nBas-nR
      do i=nC+1,nO
        do a=nO+1,nBas-nR
        do b=nO+1,nBas-nR
        do c=nO+1,nBas-nR
        do d=nO+1,nBas-nR

            eps1 = eGF3(p) + e0(i) - e0(a) - e0(b)
            eps2 = eGF3(p) + e0(i) - e0(c) - e0(d)

            Cpp(p,1) = Cpp(p,1) & 
                     + (2d0*V(p,i,a,b) - V(p,i,b,a))*V(a,b,c,d)*V(p,i,c,d)/(eps1*eps2)

        enddo
        enddo
        enddo
        enddo
      enddo
    enddo

    do p=nC+1,nBas-nR
      do i=nC+1,nO
      do j=nC+1,nO
      do k=nC+1,nO
        do a=nO+1,nBas-nR
        do b=nO+1,nBas-nR

            eps1 = eGF3(p) + e0(i) - e0(a) - e0(b)
            eps2 = e0(j)   + e0(k) - e0(a) - e0(b)

            Cpp(p,2) = Cpp(p,2) & 
                     + (2d0*V(p,i,a,b) - V(p,i,b,a))*V(a,b,j,k)*V(p,i,j,k)/(eps1*eps2)

        enddo
        enddo
      enddo
      enddo
      enddo
    enddo

    Cpp(:,3) = Cpp(:,2)

    do p=nC+1,nBas-nR
      do i=nC+1,nO
      do j=nC+1,nO
        do a=nO+1,nBas-nR
        do b=nO+1,nBas-nR
        do c=nO+1,nBas-nR

            eps1 = eGF3(p) + e0(a) - e0(i) - e0(j)
            eps2 = e0(i)   + e0(j) - e0(b) - e0(c)

            Cpp(p,4) = Cpp(p,4) & 
                     + (2d0*V(p,a,i,j) - V(p,a,j,i))*V(i,j,b,c)*V(p,a,b,c)/(eps1*eps2)
        enddo
        enddo
        enddo
      enddo
      enddo
    enddo

    Cpp(:,5) = Cpp(:,4)

    do p=nC+1,nBas-nR
      do i=nC+1,nO
      do j=nC+1,nO
      do k=nC+1,nO
      do l=nC+1,nO
        do a=nO+1,nBas-nR

            eps1 = eGF3(p) + e0(a) - e0(i) - e0(j)
            eps2 = eGF3(p) + e0(a) - e0(k) - e0(l)

            Cpp(p,6) = Cpp(p,6) & 
                     - (2d0*V(p,a,k,l) - V(p,a,l,k))*V(k,l,i,j)*V(p,a,i,j)/(eps1*eps2)
        enddo
      enddo
      enddo
      enddo
      enddo
    enddo

    ! Frequency-dependent third-order contribution: "D" terms

    Dpp(:,:) = 0d0

    do p=nC+1,nBas-nR
      do i=nC+1,nO
      do j=nC+1,nO
        do a=nO+1,nBas-nR
        do b=nO+1,nBas-nR
        do c=nO+1,nBas-nR

            eps1 = eGF3(p) + e0(i) - e0(a) - e0(b)
            eps2 = eGF3(p) + e0(j) - e0(b) - e0(c)

            Dpp(p,1) = Dpp(p,1) &
                     + V(p,i,a,b)*(V(a,j,i,c)*(    V(p,j,c,b) - 2d0*V(p,j,b,c)) &
                                 + V(a,j,c,i)*(    V(p,j,b,c) - 2d0*V(p,j,c,b)))/(eps1*eps2)

            Dpp(p,1) = Dpp(p,1) &
                     + V(p,i,b,a)*(V(a,j,i,c)*(4d0*V(p,j,b,c) - 2d0*V(p,j,c,b)) &
                                 + V(a,j,c,i)*(    V(p,j,c,b) - 2d0*V(p,j,b,c)))/(eps1*eps2)

        enddo
        enddo
        enddo
      enddo
      enddo
    enddo

    do p=nC+1,nBas-nR
      do i=nC+1,nO
      do j=nC+1,nO
        do a=nO+1,nBas-nR
        do b=nO+1,nBas-nR
        do c=nO+1,nBas-nR

            eps1 = eGF3(p) + e0(i) - e0(a) - e0(c)
            eps2 = e0(i)   + e0(j) - e0(a) - e0(b)

            Dpp(p,2) = Dpp(p,2) &
                     + V(p,i,c,a)*(V(a,b,i,j)*(4d0*V(p,b,c,j) - 2d0*V(p,b,j,c)) &
                                 + V(a,b,j,i)*(    V(p,b,j,c) - 2d0*V(p,b,c,j)))/(eps1*eps2)

            Dpp(p,2) = Dpp(p,2) &
                     + V(p,i,a,c)*(V(a,b,i,j)*(    V(p,b,j,c) - 2d0*V(p,b,c,j)) &
                                 + V(a,b,j,i)*(    V(p,b,c,j) - 2d0*V(p,b,j,c)))/(eps1*eps2)

        enddo
        enddo
        enddo
      enddo
      enddo
    enddo

    Dpp(:,3) = Dpp(:,2)

    do p=nC+1,nBas-nR
      do i=nC+1,nO
      do j=nC+1,nO
      do k=nC+1,nO
        do a=nO+1,nBas-nR
        do b=nO+1,nBas-nR

            eps1 = eGF3(p) + e0(a) - e0(j) - e0(k)
            eps2 = e0(i)   + e0(j) - e0(a) - e0(b)

            Dpp(p,4) = Dpp(p,4) &
                     + V(p,a,k,j)*(V(j,i,a,b)*(4d0*V(p,i,k,b) - 2d0*V(p,i,b,k)) &
                                 + V(j,i,b,a)*(    V(p,i,b,k) - 2d0*V(p,i,k,b)))/(eps1*eps2)

            Dpp(p,4) = Dpp(p,4) &
                     + V(p,a,j,k)*(V(j,i,a,b)*(    V(p,i,b,k) - 2d0*V(p,i,k,b)) &
                                 + V(j,i,b,a)*(    V(p,i,k,b) - 2d0*V(p,i,b,k)))/(eps1*eps2)

        enddo
        enddo
      enddo
      enddo
      enddo
    enddo

    Dpp(:,5) = Dpp(:,4)

    do p=nC+1,nBas-nR
      do i=nC+1,nO
      do j=nC+1,nO
      do k=nC+1,nO
        do a=nO+1,nBas-nR
        do b=nO+1,nBas-nR

            eps1 = eGF3(p) + e0(a) - e0(i) - e0(k)
            eps2 = eGF3(p) + e0(b) - e0(j) - e0(k)

            Dpp(p,6) = Dpp(p,6) &
                     - V(p,a,k,i)*(V(i,b,a,j)*(4d0*V(p,b,k,j) - 2d0*V(p,b,j,k)) &
                                 + V(i,b,j,a)*(    V(p,b,j,k) - 2d0*V(p,b,k,j)))/(eps1*eps2)

            Dpp(p,6) = Dpp(p,6) &
                     - V(p,a,i,k)*(V(i,b,a,j)*(    V(p,b,j,k) - 2d0*V(p,b,k,j)) &
                                 + V(i,b,j,a)*(    V(p,b,k,j) - 2d0*V(p,b,j,k)))/(eps1*eps2)

        enddo
        enddo
      enddo
      enddo
      enddo
    enddo

!   Compute renormalization factor (if required)

    Z(:) = 1d0

    if(renormalization == 0) then

      Sig3(:) = SigInf(:) &
              + Cpp(:,1) + Cpp(:,2) + Cpp(:,3) + Cpp(:,4) + Cpp(:,5) + Cpp(:,6) &
              + Dpp(:,1) + Dpp(:,2) + Dpp(:,3) + Dpp(:,4) + Dpp(:,5) + Dpp(:,6)

    elseif(renormalization == 1) then
  
      Sig3(:) = SigInf(:) &
              + Cpp(:,1) + Cpp(:,2) + Cpp(:,3) + Cpp(:,4) + Cpp(:,5) + Cpp(:,6) &
              + Dpp(:,1) + Dpp(:,2) + Dpp(:,3) + Dpp(:,4) + Dpp(:,5) + Dpp(:,6)

      Z(:) = Cpp(:,2) + Cpp(:,3) + Cpp(:,4) + Cpp(:,5) & 
           + Dpp(:,2) + Dpp(:,3) + Dpp(:,4) + Dpp(:,5)

      Z(nC+1:nBas-nR) = Z(nC+1:nBas-nR)/Sig2(nC+1:nBas-nR)
      Z(:) = 1d0/(1d0 - Z(:))
    
      Sig3(:) = Z(:)*Sig3(:)

    elseif(renormalization == 2) then

      Sig2h1p(:) = Cpp(:,4) + Cpp(:,5) + Cpp(:,6) + Dpp(:,4) + Dpp(:,5) + Dpp(:,6)
      Sig1h2p(:) = Cpp(:,1) + Cpp(:,2) + Cpp(:,3) + Dpp(:,1) + Dpp(:,2) + Dpp(:,3)

      X2h1p(:) = Cpp(:,4) + Cpp(:,5) + Dpp(:,4) + Dpp(:,5)
      X1h2p(:) = Cpp(:,2) + Cpp(:,3) + Dpp(:,2) + Dpp(:,3)
 
      X2h1p(nC+1:nBas-nR) = X2h1p(nC+1:nBas-nR)/Bpp(nC+1:nBas-nR,1)
      X1h2p(nC+1:nBas-nR) = X1h2p(nC+1:nBas-nR)/Bpp(nC+1:nBas-nR,2)

      Sig3(:) = SigInf(:) +                     &
              + 1d0/(1d0 - X2h1p(:))*Sig2h1p(:) &
              + 1d0/(1d0 - X1h2p(:))*Sig1h2p(:)
 
    elseif(renormalization == 3) then

      Sig3(:) = SigInf(:) &
              + Cpp(:,1) + Cpp(:,2) + Cpp(:,3) + Cpp(:,4) + Cpp(:,5) + Cpp(:,6) &
              + Dpp(:,1) + Dpp(:,2) + Dpp(:,3) + Dpp(:,4) + Dpp(:,5) + Dpp(:,6)

      Sig2h1p(:) = Cpp(:,4) + Cpp(:,5) + Cpp(:,6) + Dpp(:,4) + Dpp(:,5) + Dpp(:,6)
      Sig1h2p(:) = Cpp(:,1) + Cpp(:,2) + Cpp(:,3) + Dpp(:,1) + Dpp(:,2) + Dpp(:,3)

      X2h1p(:) = Cpp(:,4) + Cpp(:,5) + Dpp(:,4) + Dpp(:,5)
      X1h2p(:) = Cpp(:,2) + Cpp(:,3) + Dpp(:,2) + Dpp(:,3)

      X2h1p(nC+1:nBas-nR) = X2h1p(nC+1:nBas-nR)/Bpp(nC+1:nBas-nR,1)
      X1h2p(nC+1:nBas-nR) = X1h2p(nC+1:nBas-nR)/Bpp(nC+1:nBas-nR,2)

      Z(:) = X2h1p(:)*Sig2h1p(:) + X1h2p(:)*Sig1h2p(:)
      Z(nC+1:nBas-nR) = Z(nC+1:nBas-nR)/(Sig3(nC+1:nBas-nR) - SigInf(nC+1:nBas-nR))
      Z(:) = 1d0/(1d0 - Z(:))

      Sig3(:) = Z(:)*Sig3(:)

    endif

    ! Total third-order Green function

     eGF3(:) = e0(:) + Sig2(:) + Sig3(:)

    ! Convergence criteria

    Conv = maxval(abs(eGF3 - eOld))

    ! Print results

    call print_evGF3(nBas,nO,nSCF,Conv,e0,Z,eGF3)

    ! DIIS extrapolation

    n_diis = min(n_diis+1,max_diis)
    call DIIS_extrapolation(rcond,nBas,nBas,n_diis,error_diis,e_diis,eGF3-eOld,eGF3)

    if(abs(rcond) < 1d-15) n_diis = 0

    ! Store result for next iteration

    eOld(:) = eGF3(:)

    ! Increment

    nSCF = nSCF + 1

  enddo
!------------------------------------------------------------------------
! End main SCF loop
!------------------------------------------------------------------------

! Did it actually converge?

  if(nSCF == maxSCF+1) then

    write(*,*)
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)'                 Convergence failed                 '
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)

    stop

  endif

end subroutine 
