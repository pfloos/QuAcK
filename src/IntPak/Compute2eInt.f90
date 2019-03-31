subroutine Compute2eInt(debug,chemist_notation,iType,nShell,             &
                      ExpS,KG,DG,ExpG,                                   &
                      CenterShell,TotAngMomShell,KShell,DShell,ExpShell, &
                      np2eInt,nSigp2eInt,nc2eInt,nSigc2eInt)


! Compute various two-electron integrals

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: debug
  logical,intent(in)            :: chemist_notation
  integer,intent(in)            :: iType,nShell
  double precision,intent(in)   :: ExpS
  integer,intent(in)            :: KG
  double precision,intent(in)   :: DG(KG),ExpG(KG)
  double precision,intent(in)   :: CenterShell(maxShell,3)
  integer,intent(in)            :: TotAngMomShell(maxShell),KShell(maxShell)
  double precision,intent(in)   :: DShell(maxShell,maxK),ExpShell(maxShell,maxK)

! Local variables

  integer                       :: KBra(2),KKet(2)
  double precision              :: CenterBra(2,3),CenterKet(2,3)
  integer                       :: TotAngMomBra(2),TotAngMomKet(2)
  integer                       :: AngMomBra(2,3),AngMomKet(2,3)
  integer                       :: nShellFunctionBra(2),nShellFunctionKet(2)
  integer,allocatable           :: ShellFunctionA1(:,:),ShellFunctionA2(:,:)
  integer,allocatable           :: ShellFunctionB1(:,:),ShellFunctionB2(:,:)
  double precision              :: ExpBra(2),ExpKet(2)
  double precision              :: DBra(2),DKet(2)
  double precision              :: norm_coeff

  integer                       :: iBasA1,iBasA2,iBasB1,iBasB2
  integer                       :: iShA1,iShA2,iShB1,iShB2
  integer                       :: iShFA1,iShFA2,iShFB1,iShFB2
  integer                       :: iKA1,iKA2,iKB1,iKB2
  integer                       :: iFile

  double precision              :: p2eInt,c2eInt
  double precision              :: start_c2eInt,end_c2eInt,t_c2eInt

! Output variables

  integer,intent(out)           :: np2eInt,nSigp2eInt,nc2eInt,nSigc2eInt

  np2eInt = 0
  nSigp2eInt = 0

  nc2eInt = 0
  nSigc2eInt = 0

  iBasA1 = 0
  iBasA2 = 0
  iBasB1 = 0
  iBasB2 = 0

! Open file to write down integrals

  iFile = 0

  if(iType == 1) then

! Compute two-electron integrals over the Coulomb operator

    write(*,*) '******************************************'
    write(*,*) ' Compute two-electron repulsion integrals '
    write(*,*) '******************************************'
    write(*,*)

    iFile = 21
    open(unit=iFile,file='int/ERI.dat')

  elseif(iType == 2) then

! Compute two-electron integrals over Slater geminals

    write(*,*) '****************************************'
    write(*,*) ' Compute two-electron geminal integrals '
    write(*,*) '****************************************'
    write(*,*)

    iFile = 22
    open(unit=iFile,file='int/F12.dat')

  elseif(iType == 3) then

! Compute two-electron integrals over the Yukawa operator

    write(*,*) '***************************************'
    write(*,*) ' Compute two-electron Yukawa integrals '
    write(*,*) '***************************************'
    write(*,*)

    iFile = 23
    open(unit=iFile,file='int/Yuk.dat')

  elseif(iType == 4) then

! Compute two-electron integrals over the long-range Coulomb operator

    write(*,*) '**************************************'
    write(*,*) ' Compute long-range Coulomb integrals '
    write(*,*) '**************************************'
    write(*,*)

    iFile = 24
    open(unit=iFile,file='int/Erf.dat')

  end if

!------------------------------------------------------------------------
! Loops over shell A1
!------------------------------------------------------------------------
  do iShA1=1,nShell

    CenterBra(1,1) = CenterShell(iShA1,1)
    CenterBra(1,2) = CenterShell(iShA1,2)
    CenterBra(1,3) = CenterShell(iShA1,3)

    TotAngMomBra(1) = TotAngMomShell(iShA1)
    nShellFunctionBra(1) = (TotAngMomBra(1)*TotAngMomBra(1) + 3*TotAngMomBra(1) + 2)/2
    allocate(ShellFunctionA1(1:nShellFunctionBra(1),1:3))
    call GenerateShell(TotAngMomBra(1),nShellFunctionBra(1),ShellFunctionA1)

    KBra(1) = KShell(iShA1)

    do iShFA1=1,nShellFunctionBra(1)

      iBasA1 = iBasA1 + 1
      AngMomBra(1,1) = ShellFunctionA1(iShFA1,1)
      AngMomBra(1,2) = ShellFunctionA1(iShFA1,2)
      AngMomBra(1,3) = ShellFunctionA1(iShFA1,3)

!------------------------------------------------------------------------
! Loops over shell B1
!------------------------------------------------------------------------
      do iShB1=1,iShA1
   
        CenterKet(1,1) = CenterShell(iShB1,1)
        CenterKet(1,2) = CenterShell(iShB1,2)
        CenterKet(1,3) = CenterShell(iShB1,3)
   
        TotAngMomKet(1) = TotAngMomShell(iShB1)
        nShellFunctionKet(1) = (TotAngMomKet(1)*TotAngMomKet(1) + 3*TotAngMomKet(1) + 2)/2
        allocate(ShellFunctionB1(1:nShellFunctionKet(1),1:3))
        call GenerateShell(TotAngMomKet(1),nShellFunctionKet(1),ShellFunctionB1)
   
        KKet(1) = KShell(iShB1)
   
        do iShFB1=1,nShellFunctionKet(1)
   
          iBasB1 = iBasB1 + 1
          AngMomKet(1,1) = ShellFunctionB1(iShFB1,1)
          AngMomKet(1,2) = ShellFunctionB1(iShFB1,2)
          AngMomKet(1,3) = ShellFunctionB1(iShFB1,3)

!------------------------------------------------------------------------
! Loops over shell A2
!------------------------------------------------------------------------
          do iShA2=1,iShA1
         
            CenterBra(2,1) = CenterShell(iShA2,1)
            CenterBra(2,2) = CenterShell(iShA2,2)
            CenterBra(2,3) = CenterShell(iShA2,3)
         
            TotAngMomBra(2) = TotAngMomShell(iShA2)
            nShellFunctionBra(2) = (TotAngMomBra(2)*TotAngMomBra(2) + 3*TotAngMomBra(2) + 2)/2
            allocate(ShellFunctionA2(1:nShellFunctionBra(2),1:3))
            call GenerateShell(TotAngMomBra(2),nShellFunctionBra(2),ShellFunctionA2)
         
            KBra(2) = KShell(iShA2)
         
            do iShFA2=1,nShellFunctionBra(2)
         
              iBasA2 = iBasA2 + 1
              AngMomBra(2,1) = ShellFunctionA2(iShFA2,1)
              AngMomBra(2,2) = ShellFunctionA2(iShFA2,2)
              AngMomBra(2,3) = ShellFunctionA2(iShFA2,3)

!------------------------------------------------------------------------
! Loops over shell B2
!------------------------------------------------------------------------
              do iShB2=1,iShA2
              
                CenterKet(2,1) = CenterShell(iShB2,1)
                CenterKet(2,2) = CenterShell(iShB2,2)
                CenterKet(2,3) = CenterShell(iShB2,3)
              
                TotAngMomKet(2) = TotAngMomShell(iShB2)
                nShellFunctionKet(2) = (TotAngMomKet(2)*TotAngMomKet(2) + 3*TotAngMomKet(2) + 2)/2
                allocate(ShellFunctionB2(1:nShellFunctionKet(2),1:3))
                call GenerateShell(TotAngMomKet(2),nShellFunctionKet(2),ShellFunctionB2)
              
                KKet(2) = KShell(iShB2)
              
                do iShFB2=1,nShellFunctionKet(2)
              
                  iBasB2 = iBasB2 + 1
                  AngMomKet(2,1) = ShellFunctionB2(iShFB2,1)
                  AngMomKet(2,2) = ShellFunctionB2(iShFB2,2)
                  AngMomKet(2,3) = ShellFunctionB2(iShFB2,3)

!------------------------------------------------------------------------
!                         Loops over contraction degrees
!-------------------------------------------------------------------------
                  call cpu_time(start_c2eInt)
          
                  c2eInt = 0d0
          
                  do iKA1=1,KBra(1)
                    ExpBra(1) = ExpShell(iShA1,iKA1)
                    DBra(1) = DShell(iShA1,iKA1)*norm_coeff(ExpBra(1),AngMomBra(1,1:3))
                    do iKA2=1,KBra(2)
                      ExpBra(2) = ExpShell(iShA2,iKA2)
                      DBra(2) = DShell(iShA2,iKA2)*norm_coeff(ExpBra(2),AngMomBra(2,1:3))
                      do iKB1=1,KKet(1)
                        ExpKet(1) = ExpShell(iShB1,iKB1)
                        DKet(1) = DShell(iShB1,iKB1)*norm_coeff(ExpKet(1),AngMomKet(1,1:3))
                        do iKB2=1,KKet(2)
                          ExpKet(2) = ExpShell(iShB2,iKB2)
                          DKet(2) = DShell(iShB2,iKB2)*norm_coeff(ExpKet(2),AngMomKet(2,1:3))

                          call S2eInt(debug,iType,np2eInt,nSigp2eInt, &
                                      ExpS,KG,DG,ExpG,                &
                                      ExpBra,CenterBra,AngMomBra,     &
                                      ExpKet,CenterKet,AngMomKet,     &
                                      p2eInt)

                          c2eInt = c2eInt + DBra(1)*DBra(2)*DKet(1)*DKet(2)*p2eInt

                        end do
                      end do
                    end do
                  end do
                  call cpu_time(end_c2eInt)

                  nc2eInt = nc2eInt + 1

                  if(abs(c2eInt) > 1d-15) then

                    nSigc2eInt = nSigc2eInt + 1
                    t_c2eInt = end_c2eInt - start_c2eInt

                    if(chemist_notation) then

                      write(iFile,'(I6,I6,I6,I6,F20.15)') iBasA1,iBasB1,iBasA2,iBasB2,c2eInt
!                     write(iFile,'(F20.15,I6,I6,I6,I6)') c2eInt,iBasA1,iBasB1,iBasA2,iBasB2

                      if(debug) then
                        write(*,'(A10,1X,F16.10,1X,I6,1X,I6,1X,I6,1X,I6)') &
                          '(a1b1|a2b2) = ',c2eInt,iBasA1,iBasB1,iBasA2,iBasB2
                      end if

                    else

                      write(iFile,'(I6,I6,I6,I6,F20.15)') iBasA1,iBasA2,iBasB1,iBasB2,c2eInt
!                     write(iFile,'(F20.15,I6,I6,I6,I6)') c2eInt,iBasA1,iBasA2,iBasB1,iBasB2

                      if(debug) then
                        write(*,'(A10,1X,F16.10,1X,I6,1X,I6,1X,I6,1X,I6)') &
                          '<a1a2|b1b2> = ',c2eInt,iBasA1,iBasA2,iBasB1,iBasB2
                      end if
 
                    end if
                  end if

!------------------------------------------------------------------------
!                 End loops over contraction degrees
!------------------------------------------------------------------------
                end do
                deallocate(ShellFunctionB2)
              end do
              iBasB2 = 0
!------------------------------------------------------------------------
! End loops over shell B2
!------------------------------------------------------------------------
            end do
            deallocate(ShellFunctionA2)
          end do
          iBasA2 = 0
!------------------------------------------------------------------------
! End loops over shell A2
!------------------------------------------------------------------------
        end do
        deallocate(ShellFunctionB1)
      end do
      iBasB1 = 0
!------------------------------------------------------------------------
! End loops over shell B1
!------------------------------------------------------------------------
    end do
    deallocate(ShellFunctionA1)
  end do
  iBasA1 = 0
!------------------------------------------------------------------------
! End loops over shell A1
!------------------------------------------------------------------------
  write(*,*)

! Close files to write down integrals

  close(unit=iFile)

end subroutine Compute2eInt
