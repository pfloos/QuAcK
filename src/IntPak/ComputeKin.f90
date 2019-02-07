subroutine ComputeKin(debug,nShell,                        &
                     CenterShell,TotAngMomShell,KShell,DShell,ExpShell, &
                     npKin,nSigpKin,ncKin,nSigcKin)


! Compute one-electron kinetic integrals

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: debug
  integer,intent(in)            :: nShell
  double precision,intent(in)   :: CenterShell(maxShell,3)
  integer,intent(in)            :: TotAngMomShell(maxShell),KShell(maxShell)
  double precision,intent(in)   :: DShell(maxShell,maxK),ExpShell(maxShell,maxK)

! Local variables

  integer                       :: KA,KB
  double precision              :: CenterA(3),CenterB(3)
  integer                       :: TotAngMomA,TotAngMomB
  integer                       :: AngMomA(3),AngMomB(3)
  integer                       :: nShellFunctionA,nShellFunctionB
  integer,allocatable           :: ShellFunctionA(:,:),ShellFunctionB(:,:)
  double precision              :: ExpA,ExpB
  double precision,allocatable  :: DA,DB
  double precision              :: NormCoeff

  integer                       :: iBasA,iBasB
  integer                       :: iShA,iShB
  integer                       :: iShFA,iShFB
  integer                       :: iKA,iKB

  double precision              :: pKin,cKin
  double precision              :: start_cKin,end_cKin,t_cKin

! Output variables

  integer,intent(out)           :: npKin,nSigpKin,ncKin,nSigcKin

! Compute one-electron integrals

  write(*,*) '****************************************'
  write(*,*) ' Compute one-electron kinetic integrals '
  write(*,*) '****************************************'
  write(*,*)

  npKin = 0
  nSigpKin = 0

  ncKin = 0
  nSigcKin = 0

  iBasA = 0
  iBasB = 0

! Open file to write down integrals

  open(unit=9,file='int/Kin.dat')

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
      do iShB=1,nShell

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
!         Loops over contraction degrees
!-------------------------------------------------------------------------
          call cpu_time(start_cKin)

          cKin = 0d0

          do iKA=1,KA
            ExpA = ExpShell(iShA,iKA)
            DA = DShell(iShA,iKA)*NormCoeff(ExpA,AngMomA)
            do iKB=1,KB
              ExpB = ExpShell(iShB,iKB)
              DB = DShell(iShB,iKB)*NormCoeff(ExpB,AngMomB)

              call KinInt(npKin,nSigpKin,      &
                         ExpA,CenterA,AngMomA, &
                         ExpB,CenterB,AngMomB, &
                         pKin)

              cKin = cKin + DA*DB*pKin

            enddo
          enddo
          call cpu_time(end_cKin)

          ncKin = ncKin + 1
          if(abs(cKin) > 1d-15) then
            nSigcKin = nSigcKin + 1
            t_cKin = end_cKin - start_cKin
            write(9,'(I6,I6,F20.15)') iBasA,iBasB,cKin
            if(debug) then
              write(*,'(A10,1X,F16.10,1X,I6,1X,I6)') '(a|T|b) = ',cKin,iBasA,iBasB
            endif
          endif
!------------------------------------------------------------------------
!                 End loops over contraction degrees
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

  close(unit=9)

end subroutine ComputeKin
