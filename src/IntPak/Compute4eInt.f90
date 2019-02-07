subroutine Compute4eInt(debug,nEl,iType,nShell,ExpS,         &
                        CenterShell,TotAngMomShell,KShell,DShell,ExpShell, &
                        npErf,nSigpErf,ncErf,nSigcErf)


! Compute long-range Coulomb integrals

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: debug
  integer,intent(in)            :: nEl,iType,nShell
  double precision              :: ExpS
  double precision,intent(in)   :: CenterShell(maxShell,3)
  integer,intent(in)            :: TotAngMomShell(maxShell),KShell(maxShell)
  double precision,intent(in)   :: DShell(maxShell,maxK),ExpShell(maxShell,maxK)

! Local variables

  integer                       :: KA,KB,KC,KD
  double precision              :: CenterA(3),CenterB(3),CenterC(3),CenterD(3)
  integer                       :: TotAngMomA,TotAngMomB,TotAngMomC,TotAngMomD
  integer                       :: AngMomA(3),AngMomB(3),AngMomC(3),AngMomD(3)
  integer                       :: nShellFunctionA,nShellFunctionB, &
                                   nShellFunctionC,nShellFunctionD
  integer,allocatable           :: ShellFunctionA(:,:),ShellFunctionB(:,:), &
                                   ShellFunctionC(:,:),ShellFunctionD(:,:)
  double precision              :: ExpA,ExpB,ExpC,ExpD
  double precision,allocatable  :: DA,DB,DC,DD
  double precision              :: NormCoeff

  integer                       :: iBasA,iBasB,iBasC,iBasD
  integer                       :: iShA,iShB,iShC,iShD
  integer                       :: iShFA,iShFB,iShFC,iShFD
  integer                       :: iKA,iKB,iKC,iKD

  double precision              :: pErf,cErf
  double precision              :: start_cErf,end_cErf,t_cErf

! Output variables

  integer,intent(out)           :: npErf,nSigpErf,ncErf,nSigcErf

! Compute two-electron integrals over long-range Coulomb operator

  write(*,*) '**********************************'
  write(*,*) ' Compute three-electron integrals '
  write(*,*) '**********************************'
  write(*,*)

  npErf = 0
  nSigpErf = 0

  ncErf = 0
  nSigcErf = 0

  iBasA = 0
  iBasB = 0
  iBasC = 0
  iBasD = 0

! Open file to write down integrals

  open(unit=41,file='int/4eInt_Type1.dat')

!------------------------------------------------------------------------
! Loops over shell A
!------------------------------------------------------------------------
  do iShA=1,nShell

    CenterA(1) = CenterShell(iShA,1)
    CenterA(2) = CenterShell(iShA,2)
    CenterA(3) = CenterShell(iShA,3)

    TotAngMomA = TotAngMomShell(iShA)
    nShellFunctionA = (TotAngMomA*TotAngMomA + 3*TotAngMomA + 2)/2
    allocate(ShellFunctionA(1:nShellFunctionA,1:3))
    call GenerateShell(TotAngMomA,nShellFunctionA,ShellFunctionA)

    KA = KShell(iShA)

    do iShFA=1,nShellFunctionA

      iBasA = iBasA + 1
      AngMomA(1) = ShellFunctionA(iShFA,1)
      AngMomA(2) = ShellFunctionA(iShFA,2)
      AngMomA(3) = ShellFunctionA(iShFA,3)

!------------------------------------------------------------------------
! Loops over shell B
!------------------------------------------------------------------------
      do iShB=1,iShA

        CenterB(1) = CenterShell(iShB,1)
        CenterB(2) = CenterShell(iShB,2)
        CenterB(3) = CenterShell(iShB,3)

        TotAngMomB = TotAngMomShell(iShB)
        nShellFunctionB = (TotAngMomB*TotAngMomB + 3*TotAngMomB + 2)/2
        allocate(ShellFunctionB(1:nShellFunctionB,1:3))
        call GenerateShell(TotAngMomB,nShellFunctionB,ShellFunctionB)

        KB = KShell(iShB)

        do iShFB=1,nShellFunctionB

          iBasB = iBasB + 1
          AngMomB(1) = ShellFunctionB(iShFB,1)
          AngMomB(2) = ShellFunctionB(iShFB,2)
          AngMomB(3) = ShellFunctionB(iShFB,3)

!------------------------------------------------------------------------
! Loops over shell C
!------------------------------------------------------------------------
          do iShC=1,iShA

            CenterC(1) = CenterShell(iShC,1)
            CenterC(2) = CenterShell(iShC,2)
            CenterC(3) = CenterShell(iShC,3)

            TotAngMomC = TotAngMomShell(iShC)
            nShellFunctionC = (TotAngMomC*TotAngMomC + 3*TotAngMomC + 2)/2
            allocate(ShellFunctionC(1:nShellFunctionC,1:3))
            call GenerateShell(TotAngMomC,nShellFunctionC,ShellFunctionC)

            KC = KShell(iShC)

            do iShFC=1,nShellFunctionC

              iBasC = iBasC + 1
              AngMomC(1) = ShellFunctionC(iShFC,1)
              AngMomC(2) = ShellFunctionC(iShFC,2)
              AngMomC(3) = ShellFunctionC(iShFC,3)

!------------------------------------------------------------------------
! Loops over shell D
!------------------------------------------------------------------------
              do iShD=1,iShC

                CenterD(1) = CenterShell(iShD,1)
                CenterD(2) = CenterShell(iShD,2)
                CenterD(3) = CenterShell(iShD,3)

                TotAngMomD = TotAngMomShell(iShD)
                nShellFunctionD = (TotAngMomD*TotAngMomD + 3*TotAngMomD + 2)/2
                allocate(ShellFunctionD(1:nShellFunctionD,1:3))
                call GenerateShell(TotAngMomD,nShellFunctionD,ShellFunctionD)

                KD = KShell(iShD)

                do iShFD=1,nShellFunctionD

                  iBasD = iBasD + 1
                  AngMomD(1) = ShellFunctionD(iShFD,1)
                  AngMomD(2) = ShellFunctionD(iShFD,2)
                  AngMomD(3) = ShellFunctionD(iShFD,3)

!------------------------------------------------------------------------
!                 Loops over contraction degrees
!-------------------------------------------------------------------------
                  call cpu_time(start_cErf)

                  cErf = 0d0

                  do iKA=1,KA
                    ExpA = ExpShell(iShA,iKA)
                    DA = DShell(iShA,iKA)*NormCoeff(ExpA,AngMomA)
                    do iKB=1,KB
                      ExpB = ExpShell(iShB,iKB)
                      DB = DShell(iShB,iKB)*NormCoeff(ExpB,AngMomB)
                      do iKC=1,KC
                        ExpC = ExpShell(iShC,iKC)
                        DC = DShell(iShC,iKC)*NormCoeff(ExpC,AngMomC)
                        do iKD=1,KD
                          ExpD = ExpShell(iShD,iKD)
                          DD = DShell(iShD,iKD)*NormCoeff(ExpD,AngMomD)

! Erf module
!                          call ErfInt(debug,npErf,nSigpErf, &
!                                      ExpS,                 &
!                                      ExpA,CenterA,AngMomA, &
!                                      ExpB,CenterB,AngMomB, &
!                                      ExpC,CenterC,AngMomC, &
!                                      ExpD,CenterD,AngMomD, &
!                                      pErf)

!                          cErf = cErf + DA*DB*DC*DD*pErf

                        enddo
                      enddo
                    enddo
                  enddo
                  call cpu_time(end_cErf)

                  ncErf = ncErf + 1
                  if(abs(cErf) > 1d-15) then
                    nSigcErf = nSigcErf + 1
                    t_cErf = end_cErf - start_cErf
                    write(41,'(F20.15,I6,I6,I6,I6)') &
                      cErf,iBasA,iBasB,iBasC,iBasD
                    if(debug) then
                      write(*,'(A10,1X,F16.10,1X,I6,1X,I6,1X,I6,1X,I6)') &
                        '(ab|erf(r)/r|cd) = ',cErf,iBasA,iBasB,iBasC,iBasD
                    endif
                  endif

!------------------------------------------------------------------------
!                 End loops over contraction degrees
!------------------------------------------------------------------------
                enddo
                deallocate(ShellFunctionD)
              enddo
              iBasD = 0
!------------------------------------------------------------------------
! End loops over shell D
!------------------------------------------------------------------------
            enddo
            deallocate(ShellFunctionC)
          enddo
          iBasC = 0
!------------------------------------------------------------------------
! End loops over shell C
!------------------------------------------------------------------------
        enddo
        deallocate(ShellFunctionB)
      enddo
      iBasB = 0
!------------------------------------------------------------------------
! End loops over shell B
!------------------------------------------------------------------------
    enddo
    deallocate(ShellFunctionA)
  enddo
  iBasA = 0
!------------------------------------------------------------------------
! End loops over shell A
!------------------------------------------------------------------------
  write(*,*)

! Close files to write down integrals

  close(unit=41)

end subroutine Compute4eInt
