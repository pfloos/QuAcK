subroutine Compute3eInt(debug,iType,nShell,                                &
                        ExpS,KG,DG,ExpG,                                   &
                        CenterShell,TotAngMomShell,KShell,DShell,ExpShell, &
                        np3eInt,nSigp3eInt,nc3eInt,nSigc3eInt)


! Compute long-range Coulomb integrals

  implicit none
  include 'parameters.h'


! Input variables

  logical,intent(in)            :: debug
  integer,intent(in)            :: iType,nShell
  double precision,intent(in)   :: ExpS
  integer,intent(in)            :: KG
  double precision,intent(in)   :: DG(KG),ExpG(KG)
  double precision,intent(in)   :: CenterShell(maxShell,3)
  integer,intent(in)            :: TotAngMomShell(maxShell),KShell(maxShell)
  double precision,intent(in)   :: DShell(maxShell,maxK),ExpShell(maxShell,maxK)

! Local variables

  integer                       :: KBra(3),KKet(3)
  double precision              :: CenterBra(3,3),CenterKet(3,3)
  integer                       :: TotAngMomBra(3),TotAngMomKet(3)
  integer                       :: AngMomBra(3,3),AngMomKet(3,3)
  integer                       :: nShellFunctionBra(3),nShellFunctionKet(3)
  integer,allocatable           :: ShellFunctionA1(:,:),ShellFunctionA2(:,:),ShellFunctionA3(:,:)
  integer,allocatable           :: ShellFunctionB1(:,:),ShellFunctionB2(:,:),ShellFunctionB3(:,:)
  double precision              :: ExpBra(3),ExpKet(3)
  double precision              :: DBra(3),DKet(3)
  double precision              :: norm_coeff

  integer                       :: iBasA1,iBasA2,iBasA3,iBasB1,iBasB2,iBasB3
  integer                       :: iShA1,iShA2,iShA3,iShB1,iShB2,iShB3
  integer                       :: iShFA1,iShFA2,iShFA3,iShFB1,iShFB2,iShFB3
  integer                       :: iKA1,iKA2,iKA3,iKB1,iKB2,iKB3
  integer                       :: iFile

  double precision              :: p3eInt,c3eInt
  double precision              :: start_c3eInt,end_c3eInt,t_c3eInt

! Output variables

  integer,intent(out)           :: np3eInt,nSigp3eInt,nc3eInt,nSigc3eInt

! Compute three-electron integrals

  write(*,*) '**********************************'
  write(*,*) ' Compute three-electron integrals '
  write(*,*) '**********************************'
  write(*,*)

  np3eInt = 0
  nSigp3eInt = 0

  nc3eInt = 0
  nSigc3eInt = 0

  iBasA1 = 0
  iBasA2 = 0
  iBasA3 = 0
  iBasB1 = 0
  iBasB2 = 0
  iBasB3 = 0

! Open file to write down integrals

  iFile = 0

  if(iType == 1) then
    iFile = 31
    open(unit=iFile,file='int/3eInt_Type1.dat')
  elseif(iType == 2) then
    iFile = 32
    open(unit=iFile,file='int/3eInt_Type2.dat')
  elseif(iType == 3) then
    iFile = 33
    open(unit=iFile,file='int/3eInt_Type3.dat')
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
! Loops over shell A2
!------------------------------------------------------------------------
      do iShA2=1,nShell
  
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
! Loops over shell A3
!------------------------------------------------------------------------
          do iShA3=1,nShell
         
            CenterBra(3,1) = CenterShell(iShA3,1)
            CenterBra(3,2) = CenterShell(iShA3,2)
            CenterBra(3,3) = CenterShell(iShA3,3)
         
            TotAngMomBra(3) = TotAngMomShell(iShA3)
            nShellFunctionBra(3) = (TotAngMomBra(3)*TotAngMomBra(3) + 3*TotAngMomBra(3) + 2)/2
            allocate(ShellFunctionA3(1:nShellFunctionBra(3),1:3))
            call GenerateShell(TotAngMomBra(3),nShellFunctionBra(3),ShellFunctionA3)
         
            KBra(3) = KShell(iShA3)
         
            do iShFA3=1,nShellFunctionBra(3)
         
              iBasA3 = iBasA3 + 1
              AngMomBra(3,1) = ShellFunctionA3(iShFA3,1)
              AngMomBra(3,2) = ShellFunctionA3(iShFA3,2)
              AngMomBra(3,3) = ShellFunctionA3(iShFA3,3)

!------------------------------------------------------------------------
! Loops over shell B1
!------------------------------------------------------------------------
              do iShB1=1,nShell
   
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
! Loops over shell B2
!------------------------------------------------------------------------
                  do iShB2=1,nShell
              
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
! Loops over shell B3
!------------------------------------------------------------------------
                      do iShB3=1,nShell
                     
                        CenterKet(3,1) = CenterShell(iShB3,1)
                        CenterKet(3,2) = CenterShell(iShB3,2)
                        CenterKet(3,3) = CenterShell(iShB3,3)
                     
                        TotAngMomKet(3) = TotAngMomShell(iShB3)
                        nShellFunctionKet(3) = (TotAngMomKet(3)*TotAngMomKet(3) + 3*TotAngMomKet(3) + 2)/2
                        allocate(ShellFunctionB3(1:nShellFunctionKet(3),1:3))
                        call GenerateShell(TotAngMomKet(3),nShellFunctionKet(3),ShellFunctionB3)
                     
                        KKet(3) = KShell(iShB3)
                     
                        do iShFB3=1,nShellFunctionKet(3)
                     
                          iBasB3 = iBasB3 + 1
                          AngMomKet(3,1) = ShellFunctionB3(iShFB3,1)
                          AngMomKet(3,2) = ShellFunctionB3(iShFB3,2)
                          AngMomKet(3,3) = ShellFunctionB3(iShFB3,3)

!------------------------------------------------------------------------
!                         Loops over contraction degrees
!-------------------------------------------------------------------------
                          call cpu_time(start_c3eInt)
          
                          c3eInt = 0d0
          
                          do iKA1=1,KBra(1)
                            ExpBra(1) = ExpShell(iShA1,iKA1)
                            DBra(1) = DShell(iShA1,iKA1)*norm_coeff(ExpBra(1),AngMomBra(1,1:3))
                            do iKA2=1,KBra(2)
                              ExpBra(2) = ExpShell(iShA2,iKA2)
                              DBra(2) = DShell(iShA2,iKA2)*norm_coeff(ExpBra(2),AngMomBra(2,1:3))
                              do iKA3=1,KBra(3)
                                ExpBra(3) = ExpShell(iShA3,iKA3)
                                DBra(3) = DShell(iShA3,iKA3)*norm_coeff(ExpBra(3),AngMomBra(3,1:3))
                                do iKB1=1,KKet(1)
                                  ExpKet(1) = ExpShell(iShB1,iKB1)
                                  DKet(1) = DShell(iShB1,iKB1)*norm_coeff(ExpKet(1),AngMomKet(1,1:3))
                                  do iKB2=1,KKet(2)
                                    ExpKet(2) = ExpShell(iShB2,iKB2)
                                    DKet(2) = DShell(iShB2,iKB2)*norm_coeff(ExpKet(2),AngMomKet(2,1:3))
                                    do iKB3=1,KKet(3)
                                      ExpKet(3) = ExpShell(iShB3,iKB3)
                                      DKet(3) = DShell(iShB3,iKB3)*norm_coeff(ExpKet(3),AngMomKet(3,1:3))

                                      call S3eInt(debug,iType,np3eInt,nSigp3eInt, &
                                                  ExpS,KG,DG,ExpG,                &
                                                  ExpBra,CenterBra,AngMomBra,     &
                                                  ExpKet,CenterKet,AngMomKet,     &
                                                  p3eInt)
      
                                      c3eInt = c3eInt + DBra(1)*DBra(2)*DBra(3)*DKet(1)*DKet(2)*DKet(3)*p3eInt

                                    end do
                                  end do
                                end do
                              end do
                            end do
                          end do
                          call cpu_time(end_c3eInt)

                          nc3eInt = nc3eInt + 1
                          if(abs(c3eInt) > 1d-15) then
                            nSigc3eInt = nSigc3eInt + 1
                            t_c3eInt = end_c3eInt - start_c3eInt
                            write(iFile,'(I9,I9,I9,I9,I9,I9,F25.15)') &
                              iBasA1,iBasA2,iBasA3,iBasB1,iBasB2,iBasB3,c3eInt
                            if(.true.) then
                              write(*,'(A15,1X,I6,1X,I6,1X,I6,1X,I6,1X,I6,1X,I6,1X,F16.10)') &
                                '(a1a2a3|b1b2b3) = ',iBasA1,iBasA2,iBasA3,iBasB1,iBasB2,iBasB3,c3eInt
                            end if
                          end if

!------------------------------------------------------------------------
!                 End loops over contraction degrees
!------------------------------------------------------------------------
                        end do
                        deallocate(ShellFunctionB3)
                      end do
                      iBasB3 = 0
!------------------------------------------------------------------------
! End loops over shell B3
!------------------------------------------------------------------------
                    end do
                    deallocate(ShellFunctionB2)
                  end do
                  iBasB2 = 0
!------------------------------------------------------------------------
! End loops over shell B2
!------------------------------------------------------------------------
                end do
                deallocate(ShellFunctionB1)
              end do
              iBasB1 = 0
!------------------------------------------------------------------------
! End loops over shell B1
!------------------------------------------------------------------------
            end do
            deallocate(ShellFunctionA3)
          end do
          iBasA3 = 0
!------------------------------------------------------------------------
! End loops over shell A3
!------------------------------------------------------------------------
        end do
        deallocate(ShellFunctionA2)
      end do
      iBasA2 = 0
!------------------------------------------------------------------------
! End loops over shell A2
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

end subroutine Compute3eInt
