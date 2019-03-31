subroutine ComputeOv(debug,nBas,nShell,                 &
                     CenterShell,TotAngMomShell,KShell,DShell,ExpShell, &
                     npOv,nSigpOv,ncOv,nSigcOv,S)


! Compute one-electron overlap integrals

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: debug
  integer,intent(in)            :: nBas,nShell
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
  double precision              :: norm_coeff

  integer                       :: iBasA,iBasB
  integer                       :: iShA,iShB
  integer                       :: iShFA,iShFB
  integer                       :: iKA,iKB

  double precision              :: pOv,cOv
  double precision              :: start_cOv,end_cOv,t_cOv

! Output variables

  integer,intent(out)           :: npOv,nSigpOv,ncOv,nSigcOv
  double precision,intent(out)  :: S(nBas,nBas)


! Compute one-electron integrals

  write(*,*) '****************************************'
  write(*,*) ' Compute one-electron overlap integrals '
  write(*,*) '****************************************'
  write(*,*)

  npOv = 0
  nSigpOv = 0

  ncOv = 0
  nSigcOv = 0

  iBasA = 0
  iBasB = 0

! Open file to write down integrals

  open(unit=8,file='int/Ov.dat')

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
          call cpu_time(start_cOv)

          cOv = 0d0

          do iKA=1,KA
            ExpA = ExpShell(iShA,iKA)
            DA = DShell(iShA,iKA)*norm_coeff(ExpA,AngMomA)
            do iKB=1,KB
              ExpB = ExpShell(iShB,iKB)
              DB = DShell(iShB,iKB)*norm_coeff(ExpB,AngMomB)

              call OvInt(npOv,nSigpOv,         &
                         ExpA,CenterA,AngMomA, &
                         ExpB,CenterB,AngMomB, &
                         pOv)

              cOv = cOv + DA*DB*pOv

            end do
          end do
          call cpu_time(end_cOv)

          ncOv = ncOv + 1
          S(iBasA,iBasB) = cOv
          if(abs(cOv) > 1d-15) then
            nSigcOv = nSigcOv + 1
            t_cOv = end_cOv - start_cOv
            write(8,'(I6,I6,F20.15)') iBasA,iBasB,cOv
!           write(8,'(F20.15,I6,I6)') cOv,iBasA,iBasB
            if(debug) then
              write(*,'(A10,1X,F16.10,1X,I6,1X,I6)') '(a|b) = ',cOv,iBasA,iBasB
            end if
          end if

!------------------------------------------------------------------------
!                 End loops over contraction degrees
!------------------------------------------------------------------------
        end do
        deallocate(ShellFunctionB)
      end do
      iBasB = 0
!------------------------------------------------------------------------
! End loops over shell B
!------------------------------------------------------------------------
    end do
    deallocate(ShellFunctionA)
  end do
  iBasA = 0
!------------------------------------------------------------------------
! End loops over shell A
!------------------------------------------------------------------------
  write(*,*)

! Close files to write down integrals

  close(unit=8)

end subroutine ComputeOv
