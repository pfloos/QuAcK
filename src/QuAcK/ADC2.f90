subroutine ADC2(ispin,maxSCF,thresh,max_diis,nBas,nC,nO,nV,nR,e,ERI)

! Compute ADC(2) excitation energies: see Schirmer, Cederbaum & Walter, PRA, 28 (1983) 1237

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: ispin
  integer,intent(in)            :: maxSCF
  double precision,intent(in)   :: thresh
  integer,intent(in)            :: max_diis
  integer,intent(in)            :: nBas,nC,nO,nV,nR
  double precision,intent(in)   :: e(nBas),ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: nH,nP,nHH,nPP,nSCF,n_diis
  double precision              :: Conv
  double precision,external     :: Kronecker_delta
  double precision,allocatable  :: B_ADC(:,:),X_ADC(:,:),e_ADC(:),SigInf(:,:),G_ADC(:,:)
  double precision,allocatable  :: db_ERI(:,:,:,:),eOld(:),error_diis(:,:),e_diis(:,:)

  integer                       :: i,j,k,l
  integer                       :: a,b,c,d
  integer                       :: p,q,r,s
  integer                       :: nADC,iADC,jADC


! Hello world

  write(*,*)
  write(*,*)'***********************************'
  write(*,*)'|    2nd-order ADC calculation    |'
  write(*,*)'***********************************'
  write(*,*)

! Number of holes 

  nH  = nO - nC
  nHH = nH*(nH+1)/2

! Number of particles

  nP  = nV - nR
  nPP = nP*(nP+1)/2

  write(*,*) 'Total    states: ',nH+nP
  write(*,*) 'Hole     states: ',nH
  write(*,*) 'Particle states: ',nP

! Size of ADC(2) matrices

  nADC = nH + nP + nH*nPP + nHH*nP
  write(*,'(1X,A25,I3,A6,I6)') 'Size of ADC(2) matrices: ',nADC,' x ',nADC

! Memory allocation

  allocate(db_ERI(nBas,nBas,nBas,nBas),error_diis(nBas,max_diis),e_diis(nBas,max_diis),eOld(nADC), &
           B_ADC(nADC,nADC),X_ADC(nADC,nADC),e_ADC(nADC),G_ADC(nADC,nADC),SigInf(nADC,nADC))

! Create double-bar MO integrals

  call antisymmetrize_ERI(ispin,nBas,ERI,db_ERI)

! Initialization

  Conv            = 1d0
  nSCF            = 0
  n_diis          = 0
  e_diis(:,:)     = 0d0
  error_diis(:,:) = 0d0
  SigInf(:,:)     = 0d0
  eOld(:)         = 0d0

!------------------------------------------------------------------------
! Main SCF loop
!------------------------------------------------------------------------
!
!     | e + SigInf    (U^I)^t    (U^II)^t  |
!     |                                    |
! B = |    U^I      K^I + C^I       0      |
!     |                                    |
!     |    U^II         0      K^II + C^II |
!
!

  do while(Conv > thresh .and. nSCF < maxSCF)

    !
    ! Build ADC(2) B matrix -- Eq. (38b) -- 
    !

    write(*,'(1X,A7,1X,I4)') 'Cycle: ',nSCF

    !
    ! Diagonal part: static self-energy and epsilon
    !

    B_ADC(:,:) = 0d0
    B_ADC(nC+1:nV,nC+1:nV) = SigInf(nC+1:nV,nC+1:nV)

    jADC = 0

    do p=nC+1,nV

      jADC = jADC + 1
      B_ADC(jADC,jADC) = e(p)

    enddo

    !
    ! U matrices -- Eq. (40a) --
    !

    do p=nC+1,nV

      iADC = p - nC
      jADC = nH + nP

      ! U^I: 2p-1h -- Eqs. (40a) and (41a) -- 

      do i=nC+1,nO
        do a=nO+1,nV-nR
          do b=a,nV-nR

            jADC = jADC + 1
            B_ADC(iADC,jADC) = db_ERI(p,i,a,b)

          enddo
        enddo
      enddo

      ! U^II: 2h-1p -- Eqs. (40a) and (41b) -- 

      do i=nC+1,nO
        do j=i,nO
          do a=nO+1,nV-nR

            jADC = jADC + 1
            B_ADC(iADC,jADC) = db_ERI(p,a,i,j)

          enddo
        enddo
      enddo
   
    enddo

    !  
    ! K matrices -- Eq. (40b) -- 
    !

    ! K^I: 2p-1h -- Eqs. (40b) and (41a) -- 

    jADC = nH + nP

    do i=nC+1,nO
      do a=nO+1,nV-nR
        do b=a,nV-nR

          jADC = jADC + 1
          B_ADC(jADC,jADC) = e(a) + e(b) - e(i)

        enddo
      enddo
    enddo

    ! K^II: 2h-1p -- Eqs. (40b) and (41b) -- 

    do i=nC+1,nO
      do j=i,nO
        do a=nO+1,nV

          jADC = jADC + 1
          B_ADC(jADC,jADC) = e(i) + e(j) - e(a)

        enddo
      enddo
    enddo
   
    !
    ! C matrices -- Eq. (42c)
    !

    ! C^I: 2p-1h-TDA -- Eqs. (42a) and (42c) -- 

    iADC = nH + nP

    do i=nC+1,nO
      do a=nO+1,nV-nR
        do b=a,nV-nR

          iADC = iADC + 1
          jADC = nH + nP

          do j=nC+1,nO
            do c=nO+1,nV
              do d=c,nV-nR

                jADC = jADC + 1
                B_ADC(iADC,jADC) = B_ADC(iADC,jADC)                     &
                                 + db_ERI(a,b,c,d)*Kronecker_delta(i,j) &
                                 - db_ERI(j,b,i,d)*Kronecker_delta(a,c) &
                                 - db_ERI(j,a,i,c)*Kronecker_delta(b,d) &
                                 + db_ERI(b,a,c,d)*Kronecker_delta(i,j) &
                                 - db_ERI(j,a,i,d)*Kronecker_delta(b,c) &
                                 - db_ERI(j,b,i,c)*Kronecker_delta(a,d)

              enddo
            enddo
          enddo
    
        enddo
      enddo
    enddo
     
    ! C^II: 2p-1h-TDA -- Eqs. (42b) and (42c) -- 

    iADC = nH + nP + nH * nPP

    do i=nC+1,nO
      do j=i,nO
        do a=nO+1,nV-nR

          iADC = iADC + 1
          jADC = nH + nP + nH*nPP
    
          do k=nC+1,nO
            do l=k,nO
              do b=nO+1,nV-nR

                jADC = jADC + 1
                B_ADC(iADC,jADC) = B_ADC(iADC,jADC)                     &
                                 - db_ERI(i,j,k,l)*Kronecker_delta(a,b) &
                                 + db_ERI(b,j,a,l)*Kronecker_delta(i,k) &
                                 + db_ERI(b,i,a,k)*Kronecker_delta(j,l) &
                                 - db_ERI(j,i,k,l)*Kronecker_delta(a,b) &
                                 + db_ERI(b,i,a,l)*Kronecker_delta(j,k) &
                                 + db_ERI(b,j,a,k)*Kronecker_delta(i,l)

              enddo
            enddo
          enddo
    
        enddo
      enddo
    enddo
   
    ! Fold B onto itself

    do iADC=1,nADC
      do jADC=iADC+1,nADC

        B_ADC(jADC,iADC) = B_ADC(iADC,jADC)

      enddo
    enddo
  
   ! Diagonalize B to obtain X and E -- Eq. (38a) -- 
 
    X_ADC(:,:) = B_ADC(:,:)
    call diagonalize_matrix(nADC,X_ADC,e_ADC)
   
    ! print results

   
    write(*,*) '================================='
    write(*,*) 'ADC(2) excitation energies (eV)'

    do iADC=1,nADC

      if(NORM2(X_ADC(1:nH+nP,iADC)) > 0.1d0 ) & 
        write(*,'(2(2X,F12.6))') e_ADC(iADC)*HaToeV,NORM2(X_ADC(1:nH+nP,iADC))

    enddo

    write(*,*) '================================='
 
    ! Convergence criteria

     Conv = maxval(abs(e_ADC - eOld))

    ! Store result for next iteration

    eOld(:) = e_ADC(:)

    ! Compute W -- Eq (11) --

    SigInf(:,:) = 0d0

    do i=nC+1,nO
      do p=nC+1,nV-nR
        do q=nC+1,nV-nR

          SigInf(p,q) = SigInf(p,q) - db_ERI(p,i,q,i)

        enddo
      enddo
    enddo

    ! Compute the one-particle Greeen function -- Eq. (28) -- 

    G_ADC(:,:) = 0d0

    do iADC=1,nADC

      if(e_ADC(iADC) > 0d0 ) cycle

      do p=nC+1,nV-nR
        do q=nC+1,nV-nR

          G_ADC(p,q) = G_ADC(p,q) + X_ADC(p,iADC)*X_ADC(q,iADC)

        enddo
      enddo

    enddo

    ! Compute static self-energy for next iteration -- Eq. (25) --
    
    do p=nC+1,nV-nR
      do q=nC+1,nV-nR
        do r=nC+1,nV-nR
          do s=nC+1,nV-nR

            SigInf(p,q) = SigInf(p,q) + db_ERI(p,r,q,s)*G_ADC(r,s)

          enddo
        enddo
      enddo
    enddo

    ! Print results

!   call print_ADC2(nBas,nO,nSCF,Conv,e,eADC)

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

  endif

end subroutine ADC2
