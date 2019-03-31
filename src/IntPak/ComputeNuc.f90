subroutine ComputeNuc(debug,nShell,                        &
                      CenterShell,TotAngMomShell,KShell,DShell,ExpShell, &
                      NAtoms,ZNuc,XYZAtoms,                              &
                      npNuc,nSigpNuc,ncNuc,nSigcNuc)


! Compute electron repulsion integrals

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: debug
  integer,intent(in)            :: nShell
  double precision,intent(in)   :: CenterShell(maxShell,3)
  integer,intent(in)            :: TotAngMomShell(maxShell),KShell(maxShell)
  double precision,intent(in)   :: DShell(maxShell,maxK),ExpShell(maxShell,maxK)
  integer                       :: NAtoms
  double precision              :: ZNuc(NAtoms),XYZAtoms(NAtoms,3)

! Local variables

  integer                       :: KA,KB
  double precision              :: CenterA(3),CenterB(3),CenterC(3)
  integer                       :: TotAngMomA,TotAngMomB
  integer                       :: AngMomA(3),AngMomB(3)
  integer                       :: nShellFunctionA,nShellFunctionB
  integer,allocatable           :: ShellFunctionA(:,:),ShellFunctionB(:,:)
  double precision              :: ExpA,ExpB,ZC
  double precision,allocatable  :: DA,DB
  double precision              :: norm_coeff

  integer                       :: iBasA,iBasB
  integer                       :: iShA,iShB,iNucC
  integer                       :: iShFA,iShFB
  integer                       :: iKA,iKB

  double precision              :: pNuc,cNuc
  double precision              :: start_cNuc,end_cNuc,t_cNuc

! Output variables

  integer,intent(out)           :: npNuc,nSigpNuc,ncNuc,nSigcNuc

! Compute one-electron nuclear attraction integrals

  write(*,*) '***************************************************'
  write(*,*) ' Compute one-electron nuclear attraction integrals '
  write(*,*) '***************************************************'
  write(*,*)

  npNuc = 0
  nSigpNuc = 0

  ncNuc = 0
  nSigcNuc = 0

  iBasA = 0
  iBasB = 0
  iNucC = 0

! Open file to write down integrals

  open(unit=10,file='int/Nuc.dat')

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
! Loops over nuclear centers
!------------------------------------------------------------------------
          call cpu_time(start_cNuc)

          cNuc = 0d0

          do iNucC=1,NAtoms

            CenterC(1) = XYZAtoms(iNucC,1)
            CenterC(2) = XYZAtoms(iNucC,2)
            CenterC(3) = XYZAtoms(iNucC,3)

            ZC = ZNuc(iNucC)

!------------------------------------------------------------------------
!             Loops over contraction degrees
!-------------------------------------------------------------------------

            do iKA=1,KA
              ExpA = ExpShell(iShA,iKA)
              DA = DShell(iShA,iKA)*norm_coeff(ExpA,AngMomA)
              do iKB=1,KB
                ExpB = ExpShell(iShB,iKB)
                DB = DShell(iShB,iKB)*norm_coeff(ExpB,AngMomB)

                call NucInt(debug,npNuc,nSigpNuc, &
                            ExpA,CenterA,AngMomA, &
                            ExpB,CenterB,AngMomB, &
                            CenterC,              &
                            pNuc)

                cNuc = cNuc - DA*DB*ZC*pNuc

              end do
            end do
!------------------------------------------------------------------------
!           End loops over contraction degrees
!------------------------------------------------------------------------
          end do
          call cpu_time(end_cNuc)
!------------------------------------------------------------------------
! End loops over nuclear centers C
!------------------------------------------------------------------------

          ncNuc = ncNuc + 1
          if(abs(cNuc) > 1d-15) then
            nSigcNuc = nSigcNuc + 1
            t_cNuc = end_cNuc - start_cNuc
            write(10,'(I6,I6,F20.15)') iBasA,iBasB,cNuc
!           write(10,'(F20.15,I6,I6)') cNuc,iBasA,iBasB
            if(debug) then
              write(*,'(A10,1X,F16.10,1X,I6,1X,I6)') '(a|V|b) = ',cNuc,iBasA,iBasB
              write(*,*)
            end if
          end if

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

  close(unit=10)

end subroutine ComputeNuc
