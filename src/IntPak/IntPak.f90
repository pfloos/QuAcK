program IntPak

  implicit none
  include 'parameters.h'

  logical                       :: debug
  logical                       :: chemist_notation

  logical                       :: doOv 
  logical                       :: doKin
  logical                       :: doNuc

  logical                       :: doERI
  logical                       :: doF12
  logical                       :: doYuk
  logical                       :: doErf

  logical                       :: do3eInt(n3eInt)
  logical                       :: do4eInt(n4eInt)

  integer                       :: nNuc,nBas,iType
  integer                       :: nEl,nO,nV,nC,nR
  double precision              :: ExpS
  double precision              :: ENuc
  integer                       :: KG
  double precision,allocatable  :: DG(:),ExpG(:)
  double precision,allocatable  :: ZNuc(:),rNuc(:,:)

  integer                       :: nShell
  integer,allocatable           :: TotAngMomShell(:),KShell(:)
  double precision,allocatable  :: CenterShell(:,:),DShell(:,:),ExpShell(:,:)

  double precision              :: start_1eInt(n1eInt),end_1eInt(n1eInt),t_1eInt(n1eInt)
  double precision              :: start_2eInt(n2eInt),end_2eInt(n2eInt),t_2eInt(n2eInt)
  double precision              :: start_3eInt(n3eInt),end_3eInt(n3eInt),t_3eInt(n3eInt)
  double precision              :: start_4eInt(n4eInt),end_4eInt(n4eInt),t_4eInt(n4eInt)

  integer                       :: np1eInt(n1eInt),nSigp1eInt(n1eInt),nc1eInt(n1eInt),nSigc1eInt(n1eInt)
  integer                       :: np2eInt(n2eInt),nSigp2eInt(n2eInt),nc2eInt(n2eInt),nSigc2eInt(n2eInt)
  integer                       :: np3eInt(n3eInt),nSigp3eInt(n3eInt),nc3eInt(n3eInt),nSigc3eInt(n3eInt)
  integer                       :: np4eInt(n4eInt),nSigp4eInt(n4eInt),nc4eInt(n4eInt),nSigc4eInt(n4eInt)

  double precision,allocatable  :: S(:,:)


! Hello World

  write(*,*)
  write(*,*) '********************************'
  write(*,*) '*            IntPak            *'
  write(*,*) '* Integral Package for dummies *'
  write(*,*) '********************************'
  write(*,*)

! Read options for integral calculations

  call read_options(debug,chemist_notation,ExpS,doOv,doKin,doNuc,doERI,doF12,doYuk,doErf,do3eInt,do4eInt)

!------------------------------------------------------------------------
! Read input information
!------------------------------------------------------------------------

! Read number of atoms, number of electrons of the system
! nO   = number of occupied orbitals
! nV   = number of virtual orbitals (see below)
! nBas = number of basis functions (see below)
!      = nO + nV

  call read_molecule(nNuc,nEl,nO,nC,nR)

  allocate(ZNuc(1:nNuc),rNuc(1:nNuc,1:3))

! Read geometry

  call read_geometry(nNuc,ZNuc,rNuc,ENuc)

  allocate(CenterShell(1:maxShell,1:3),TotAngMomShell(1:maxShell),KShell(1:maxShell), &
           DShell(1:maxShell,1:maxK),ExpShell(1:maxShell,1:maxK))

  call read_basis(nNuc,rNuc,nBas,nO,nV,nShell,TotAngMomShell,CenterShell,KShell,DShell,ExpShell)

!------------------------------------------------------------------------
! Memory allocation
!------------------------------------------------------------------------
  allocate(S(1:nBas,1:nBas))

!------------------------------------------------------------------------
! Compute one-electron overlap integrals
!------------------------------------------------------------------------
  if(doOv) then

    iType = 1

    call cpu_time(start_1eInt(iType))
    call ComputeOv(debug,nBas,nShell,                  &
                    CenterShell,TotAngMomShell,KShell,DShell,ExpShell, &
                    np1eInt(iType),nSigp1eInt(iType),nc1eInt(iType),nSigc1eInt(iType),S)
    call cpu_time(end_1eInt(iType))

    write(*,'(A65,1X,I9)') 'Total number of primitive overlap integrals = ',np1eInt(iType)
    write(*,'(A65,1X,I9)') 'Number of significant primitive overlap integrals = ',nSigp1eInt(iType)

    write(*,'(A65,1X,I9)') 'Total number of contracted overlap integrals = ',nc1eInt(iType)
    write(*,'(A65,1X,I9)') 'Number of significant contracted overlap integrals = ',nSigc1eInt(iType)

    write(*,*)

    t_1eInt(iType) = end_1eInt(iType) - start_1eInt(iType)
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time = ',t_1eInt(iType),' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Compute one-electron kinetic integrals
!------------------------------------------------------------------------

  if(doKin) then

    iType = 2

    call cpu_time(start_1eInt(iType))
    call ComputeKin(debug,nShell,                        &
                    CenterShell,TotAngMomShell,KShell,DShell,ExpShell, &
                    np1eInt(iType),nSigp1eInt(iType),nc1eInt(iType),nSigc1eInt(iType))
    call cpu_time(end_1eInt(iType))

    write(*,'(A65,1X,I9)') 'Total number of primitive kinetic integrals = ',np1eInt(iType)
    write(*,'(A65,1X,I9)') 'Number of significant primitive kinetic integrals = ',nSigp1eInt(iType)

    write(*,'(A65,1X,I9)') 'Total number of contracted kinetic integrals = ',nc1eInt(iType)
    write(*,'(A65,1X,I9)') 'Number of significant contracted kinetic integrals = ',nSigc1eInt(iType)

    write(*,*)

    t_1eInt(iType) = end_1eInt(iType) - start_1eInt(iType)
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time = ',t_1eInt(iType),' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Compute one-electron nuclear attraction integrals
!------------------------------------------------------------------------

  if(doNuc) then

    iType = 3

    call cpu_time(start_1eInt(iType))
    call ComputeNuc(debug,nShell,                        &
                    CenterShell,TotAngMomShell,KShell,DShell,ExpShell, &
                    nNuc,ZNuc,rNuc,                              &
                    np1eInt(iType),nSigp1eInt(iType),nc1eInt(iType),nSigc1eInt(iType))
    call cpu_time(end_1eInt(iType))

    write(*,'(A65,1X,I9)') 'Total number of primitive nuclear integrals = ',np1eInt(iType)
    write(*,'(A65,1X,I9)') 'Number of significant primitive nuclear integrals = ',nSigp1eInt(iType)

    write(*,'(A65,1X,I9)') 'Total number of contracted nuclear integrals = ',nc1eInt(iType)
    write(*,'(A65,1X,I9)') 'Number of significant contracted nuclear integrals = ',nSigc1eInt(iType)

    write(*,*)

    t_1eInt(iType) = end_1eInt(iType) - start_1eInt(iType)
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time = ',t_1eInt(iType),' seconds'
    write(*,*)

  end if

!------------------------------------------------------------------------
! Compute ERIs
!------------------------------------------------------------------------

  if(doERI) then
    
    iType = 1
    KG   = 1
    allocate(DG(1:KG),ExpG(1:KG))
    DG   = (/ 1d0 /) 
    ExpG = (/ 0d0 /)

    call cpu_time(start_2eInt(iType))
    call Compute2eInt(debug,chemist_notation,iType,nShell,               &
                      ExpS,KG,DG,ExpG,                                   &
                      CenterShell,TotAngMomShell,KShell,DShell,ExpShell, &
                      np2eInt(iType),nSigp2eInt(iType),nc2eInt(iType),nSigc2eInt(iType))
    call cpu_time(end_2eInt(iType))

    write(*,'(A65,1X,I9)') 'Total number of primitive ERIs = ',np2eInt(iType)
    write(*,'(A65,1X,I9)') 'Number of significant primitive ERIs = ',nSigp2eInt(iType)

    write(*,'(A65,1X,I9)') 'Total number of contracted ERIs = ',nc2eInt(iType)
    write(*,'(A65,1X,I9)') 'Number of significant contracted ERIs = ',nSigc2eInt(iType)

    write(*,*)

    t_2eInt(iType) = end_2eInt(iType) - start_2eInt(iType)
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time = ',t_2eInt(iType),' seconds'
    write(*,*)

    deallocate(DG,ExpG)

  end if

!------------------------------------------------------------------------
! Compute F12 two-electron integrals 
!------------------------------------------------------------------------

  if(doF12) then

    iType = 2
    KG   = 6
    allocate(DG(1:KG),ExpG(1:KG))
    DG   = (/ 0.3144d0, 0.3037d0, 0.1681d0, 0.09811d0, 0.06024d0, 0.03726d0 /)
    ExpG = (/ 0.2209d0, 1.004d0,  3.622d0, 12.16d0,   45.87d0,  254.4d0     /)

!   KG = 10
!   allocate(DG(1:KG),ExpG(1:KG))

!   DG   = (/ 220.983854141, 18.52358977132, 4.81060044582, 1.892812227999, &
!     0.920641976732, 0.505281134191, 0.295757471525, 0.1753021140139, &
!     0.0969611396173, 0.0386163391551 /)
!   ExpG = (/ 5722.54799330, 191.0413784782, 27.4417708701, 6.39987966572, &
!     1.82203908762, 0.548835646170, 0.156252937904, 0.036440796942, &
!     0.0052344680925, 0.00017474733304 /)

!    KG = 20
!    allocate(DG(1:KG),ExpG(1:KG))

!  DG = (/ 841.88478132, 70.590185207, 18.3616020768, 7.2608642093, &
!3.57483416444, 2.01376031082, 1.24216542801, 0.81754348620, &
!0.564546514023, 0.404228610699, 0.297458536575, 0.223321219537, &
!0.169933732064, 0.130190978230, 0.099652303426, 0.075428246546, &
!0.0555635614051, 0.0386791283055, 0.0237550435652, 0.0100062783874 /)

! ExpG = (/84135.654509, 2971.58727634, 474.716025959, 130.676724560, &
!47.3938388887, 20.2078651631, 9.5411021938, 4.8109546955, &
!2.52795733067, 1.35894103210, 0.73586710268, 0.39557629706, &
!0.20785895177, 0.104809693858, 0.049485682527, 0.021099788990, &
!0.007652472186, 0.0021065225215, 0.0003365204879, 0.00001188556749 /)



    call cpu_time(start_2eInt(iType))
    call Compute2eInt(debug,chemist_notation,iType,nShell,               & 
                      ExpS,KG,DG,ExpG,                                   &
                      CenterShell,TotAngMomShell,KShell,DShell,ExpShell, &
                      np2eInt(iType),nSigp2eInt(iType),nc2eInt(iType),nSigc2eInt(iType))
    call cpu_time(end_2eInt(iType))

    write(*,'(A65,1X,I9)') 'Total number of primitive geminal integrals = ',np2eInt(iType)
    write(*,'(A65,1X,I9)') 'Number of significant primitive geminal integrals = ',nSigp2eInt(iType)

    write(*,'(A65,1X,I9)') 'Total number of contracted geminal integrals = ',nc2eInt(iType)
    write(*,'(A65,1X,I9)') 'Number of significant contracted geminal integrals = ',nSigc2eInt(iType)

    write(*,*)

    t_2eInt(iType) = end_2eInt(iType) - start_2eInt(iType)
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time = ',t_2eInt(iType),' seconds'
    write(*,*)

    deallocate(DG,ExpG)

  end if

!------------------------------------------------------------------------
! Compute Yukawa two-electron integrals 
!------------------------------------------------------------------------

  if(doYuk) then

    iType = 3
    KG   = 6
    allocate(DG(1:KG),ExpG(1:KG))
    DG   = (/ 0.3144d0, 0.3037d0, 0.1681d0, 0.09811d0, 0.06024d0, 0.03726d0 /)
    ExpG = (/ 0.2209d0, 1.004d0,  3.622d0, 12.16d0,   45.87d0,  254.4d0     /)

    call cpu_time(start_2eInt(iType))
    call Compute2eInt(debug,chemist_notation,iType,nShell,               &
                      ExpS,KG,DG,ExpG,                                   &
                      CenterShell,TotAngMomShell,KShell,DShell,ExpShell, &
                      np2eInt(iType),nSigp2eInt(iType),nc2eInt(iType),nSigc2eInt(iType))
    call cpu_time(end_2eInt(iType))

    write(*,'(A65,1X,I9)') 'Total number of primitive Yukawa integrals = ',np2eInt(iType)
    write(*,'(A65,1X,I9)') 'Number of significant primitive Yukawa integrals = ',nSigp2eInt(iType)

    write(*,'(A65,1X,I9)') 'Total number of contracted Yukawa integrals = ',nc2eInt(iType)
    write(*,'(A65,1X,I9)') 'Number of significant contracted Yukawa integrals = ',nSigc2eInt(iType)

    write(*,*)

    t_2eInt(iType) = end_2eInt(iType) - start_2eInt(iType)
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time = ',t_2eInt(iType),' seconds'
    write(*,*)

    deallocate(DG,ExpG)

  end if

!------------------------------------------------------------------------
! Compute long-range Coulomb two-electron integrals 
!------------------------------------------------------------------------

  if(doErf) then

    iType = 4
    KG   = 1
    allocate(DG(1:KG),ExpG(1:KG))
    DG   = (/ 1d0 /)
    ExpG = (/ 1d0 /)
    ExpS = ExpS*ExpS

    call cpu_time(start_2eInt(iType))
    call Compute2eInt(debug,chemist_notation,iType,nShell,               &
                      ExpS,KG,DG,ExpG,                                   &
                      CenterShell,TotAngMomShell,KShell,DShell,ExpShell, &
                      np2eInt(iType),nSigp2eInt(iType),nc2eInt(iType),nSigc2eInt(iType))
    call cpu_time(end_2eInt(iType))

    write(*,'(A65,1X,I9)') 'Total number of primitive long-range Coulomb integrals = ',np2eInt(iType)
    write(*,'(A65,1X,I9)') 'Number of significant primitive long-range Coulomb integrals = ',nSigp2eInt(iType)

    write(*,'(A65,1X,I9)') 'Total number of contracted long-range Coulomb integrals = ',nc2eInt(iType)
    write(*,'(A65,1X,I9)') 'Number of significant contracted long-range Coulomb integrals = ',nSigc2eInt(iType)

    write(*,*)

    t_2eInt(iType) = end_2eInt(iType) - start_2eInt(iType)
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time = ',t_2eInt(iType),' seconds'
    write(*,*)

    deallocate(DG,ExpG)

  end if

!------------------------------------------------------------------------
! Compute three-electron integrals: Type 1 => chain C12 S23
!------------------------------------------------------------------------

  if(do3eInt(1)) then

    iType = 1
!   KG   = 1
     KG   = 6
    allocate(DG(1:KG),ExpG(1:KG))
!   DG   = (/ 1d0 /)
!   ExpG = (/ 1d0 /)
    DG   = (/ 0.3144d0, 0.3037d0, 0.1681d0, 0.09811d0, 0.06024d0, 0.03726d0 /)
    ExpG = (/ 0.2209d0, 1.004d0,  3.622d0, 12.16d0,   45.87d0,  254.4d0     /)

    call cpu_time(start_3eInt(iType))
    call Compute3eInt(debug,iType,nShell,                            &
                      ExpS,KG,DG,ExpG,                                   &
                      CenterShell,TotAngMomShell,KShell,DShell,ExpShell, &
                      np3eInt(iType),nSigp3eInt(iType),nc3eInt(iType),nSigc3eInt(iType))
    call cpu_time(end_3eInt(iType))

    write(*,'(A65,1X,I9)') 'Total number of primitive f23/r12 integrals = ',np3eInt(iType)
    write(*,'(A65,1X,I9)') 'Number of significant primitive f23/r12 integrals = ',nSigp3eInt(iType)

    write(*,'(A65,1X,I9)') 'Total number of contracted f23/r12 integrals = ',nc3eInt(iType)
    write(*,'(A65,1X,I9)') 'Number of significant contracted f23/r12 integrals = ',nSigc3eInt(iType)

    write(*,*)

    t_3eInt(iType) = end_3eInt(iType) - start_3eInt(iType)
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time = ',t_3eInt(iType),' seconds'
    write(*,*)

    deallocate(DG,ExpG)

  end if

!------------------------------------------------------------------------
! Compute three-electron integrals: Type 2 => cyclic C12 S13 S23
!------------------------------------------------------------------------

  if(do3eInt(2)) then

    iType = 2
    KG   = 6
    allocate(DG(1:KG),ExpG(1:KG))
    DG   = (/ 0.3144d0, 0.3037d0, 0.1681d0, 0.09811d0, 0.06024d0, 0.03726d0 /)
    ExpG = (/ 0.2209d0, 1.004d0,  3.622d0, 12.16d0,   45.87d0,  254.4d0     /)

    call cpu_time(start_3eInt(iType))
    call Compute3eInt(debug,iType,nShell,              &
                      ExpS,KG,DG,ExpG,                                   &
                      CenterShell,TotAngMomShell,KShell,DShell,ExpShell, &
                      np3eInt(iType),nSigp3eInt(iType),nc3eInt(iType),nSigc3eInt(iType))
    call cpu_time(end_3eInt(iType))

    write(*,'(A65,1X,I9)') 'Total number of primitive f13.f23/r12 integrals = ',np3eInt(iType)
    write(*,'(A65,1X,I9)') 'Number of significant primitive f13.f23/r12 integrals = ',nSigp3eInt(iType)

    write(*,'(A65,1X,I9)') 'Total number of contracted f13.f23/r12 integrals = ',nc3eInt(iType)
    write(*,'(A65,1X,I9)') 'Number of significant contracted f13.f23/r12 integrals = ',nSigc3eInt(iType)

    write(*,*)

    t_3eInt(iType) = end_3eInt(iType) - start_3eInt(iType)
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time = ',t_3eInt(iType),' seconds'
    write(*,*)

    deallocate(DG,ExpG)

  end if

!------------------------------------------------------------------------
! Compute three-electron integrals: Type 3 => chain S13 S23
!------------------------------------------------------------------------

  if(do3eInt(3)) then

    iType = 3
    KG   = 6
    allocate(DG(1:KG),ExpG(1:KG))
    DG   = (/ 0.3144d0, 0.3037d0, 0.1681d0, 0.09811d0, 0.06024d0, 0.03726d0 /)
    ExpG = (/ 0.2209d0, 1.004d0,  3.622d0, 12.16d0,   45.87d0,  254.4d0     /)

    call cpu_time(start_3eInt(iType))
    call Compute3eInt(debug,iType,nShell,              &
                      ExpS,KG,DG,ExpG,                                   &
                      CenterShell,TotAngMomShell,KShell,DShell,ExpShell, &
                      np3eInt(iType),nSigp3eInt(iType),nc3eInt(iType),nSigc3eInt(iType))
    call cpu_time(end_3eInt(iType))

    write(*,'(A65,1X,I9)') 'Total number of primitive f13.f23 integrals = ',np3eInt(iType)
    write(*,'(A65,1X,I9)') 'Number of significant primitive f13.f23 integrals = ',nSigp3eInt(iType)

    write(*,'(A65,1X,I9)') 'Total number of contracted f13.f23 integrals = ',nc3eInt(iType)
    write(*,'(A65,1X,I9)') 'Number of significant contracted f13.f23 integrals = ',nSigc3eInt(iType)

    write(*,*)

    t_3eInt(iType) = end_3eInt(iType) - start_3eInt(iType)
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time = ',t_3eInt(iType),' seconds'
    write(*,*)

    deallocate(DG,ExpG)

  end if

!------------------------------------------------------------------------
! Compute four-electron integrals: Type 1 => chain C12 S14 S23
!------------------------------------------------------------------------

  if(do4eInt(1)) then

    iType = 1
    KG   = 6
    allocate(DG(1:KG),ExpG(1:KG))
    DG   = (/ 0.3144d0, 0.3037d0, 0.1681d0, 0.09811d0, 0.06024d0, 0.03726d0 /)
    ExpG = (/ 0.2209d0, 1.004d0,  3.622d0, 12.16d0,   45.87d0,  254.4d0     /)

    call cpu_time(start_4eInt(iType))
!    call Compute4eInt(debug,iType,nShell,ExpS,         &
!                      ExpS,KG,DG,ExpG,                                   &
!                      CenterShell,TotAngMomShell,KShell,DShell,ExpShell, &
!                      np4eInt(iType),nSigp4eInt(iType),nc4eInt(iType),nSigc4eInt(iType))
    call cpu_time(end_4eInt(iType))

    write(*,'(A65,1X,I9)') 'Total number of primitive f14.f23/r12 integrals = ',np4eInt(iType)
    write(*,'(A65,1X,I9)') 'Number of significant primitive f14.f23/r12 integrals = ',nSigp4eInt(iType)

    write(*,'(A65,1X,I9)') 'Total number of contracted f14.f23/r12 integrals = ',nc4eInt(iType)
    write(*,'(A65,1X,I9)') 'Number of significant contracted f14.f23/r12 integrals = ',nSigc4eInt(iType)

    write(*,*)

    t_4eInt(iType) = end_4eInt(iType) - start_4eInt(iType)
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time = ',t_4eInt(iType),' seconds'
    write(*,*)

    deallocate(DG,ExpG)

  end if

!------------------------------------------------------------------------
! Compute four-electron integrals: Type 2 => trident C12 S13 S14
!------------------------------------------------------------------------

  if(do4eInt(2)) then

    iType = 2
    KG   = 6
    DG   = (/ 0.3144d0, 0.3037d0, 0.1681d0, 0.09811d0, 0.06024d0, 0.03726d0 /)
    ExpG = (/ 0.2209d0, 1.004d0,  3.622d0, 12.16d0,   45.87d0,  254.4d0     /)

    call cpu_time(start_4eInt(iType))
!    call Compute4eInt(debug,iType,nShell,ExpS,         &
!                      ExpS,KG,DG,ExpG,                                   &
!                      CenterShell,TotAngMomShell,KShell,DShell,ExpShell, &
!                      np4eInt(iType),nSigp4eInt(iType),nc4eInt(iType),nSigc4eInt(iType))
    call cpu_time(end_4eInt(iType))

    write(*,'(A65,1X,I9)') 'Total number of primitive f13.f14/r12 integrals = ',np4eInt(iType)
    write(*,'(A65,1X,I9)') 'Number of significant primitive f13.f14/r12 integrals = ',nSigp4eInt(iType)

    write(*,'(A65,1X,I9)') 'Total number of contracted f13.f14/r12 integrals = ',nc4eInt(iType)
    write(*,'(A65,1X,I9)') 'Number of significant contracted f13.f14/r12 integrals = ',nSigc4eInt(iType)

    write(*,*)

    t_4eInt(iType) = end_4eInt(iType) - start_4eInt(iType)
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time = ',t_4eInt(iType),' seconds'
    write(*,*)

    deallocate(DG,ExpG)

  end if

!------------------------------------------------------------------------
! Compute four-electron integrals: Type 3 => chain C12 S13 S34
!------------------------------------------------------------------------

  if(do4eInt(3)) then

    iType = 3
    KG   = 6
    allocate(DG(1:KG),ExpG(1:KG))
    DG   = (/ 0.3144d0, 0.3037d0, 0.1681d0, 0.09811d0, 0.06024d0, 0.03726d0 /)
    ExpG = (/ 0.2209d0, 1.004d0,  3.622d0, 12.16d0,   45.87d0,  254.4d0     /)

    call cpu_time(start_4eInt(iType))
!    call Compute4eInt(debug,iType,nShell,              &
!                      ExpS,KG,DG,ExpG,                                   &
!                      CenterShell,TotAngMomShell,KShell,DShell,ExpShell, &
!                      np4eInt(iType),nSigp4eInt(iType),nc4eInt(iType),nSigc4eInt(iType))
    call cpu_time(end_4eInt(iType))

    write(*,'(A65,1X,I9)') 'Total number of primitive f13.f34/r12 integrals = ',np4eInt(iType)
    write(*,'(A65,1X,I9)') 'Number of significant primitive f13.f34/r12 integrals = ',nSigp4eInt(iType)

    write(*,'(A65,1X,I9)') 'Total number of contracted f13.f34/r12 integrals = ',nc4eInt(iType)
    write(*,'(A65,1X,I9)') 'Number of significant contracted f13.f34/r12 integrals = ',nSigc4eInt(iType)

    write(*,*)

    t_4eInt(iType) = end_4eInt(iType) - start_4eInt(iType)
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time = ',t_4eInt(iType),' seconds'
    write(*,*)

    deallocate(DG,ExpG)

  end if
!------------------------------------------------------------------------
! End of IntPak
!------------------------------------------------------------------------
end program IntPak
