subroutine read_auxiliary_basis(NAtoms,XYZAtoms,nShell,CenterShell, &
                        TotAngMomShell,KShell,DShell,ExpShell)

! Read auxiliary basis set information

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: NAtoms
  double precision,intent(in)   :: XYZAtoms(NAtoms,3)

! Local variables

  integer                       :: nShAt,iAt
  integer                       :: i,j,k
  character                     :: shelltype

! Output variables

  integer,intent(out)           :: nShell
  double precision,intent(out)  :: CenterShell(maxShell,3)
  integer,intent(out)           :: TotAngMomShell(maxShell),KShell(maxShell)
  double precision,intent(out)  :: DShell(maxShell,maxK),ExpShell(maxShell,maxK)

!------------------------------------------------------------------------
! Primary basis set information
!------------------------------------------------------------------------

! Open file with basis set specification

  open(unit=2,file='input/basis')

! Read basis information

  write(*,'(A28)') 'Gaussian basis set'
  write(*,'(A28)') '------------------'

  nShell = 0
  do i=1,NAtoms
    read(2,*) iAt,nShAt
    write(*,'(A28,1X,I16)') 'Atom n. ',iAt
    write(*,'(A28,1X,I16)') 'number of shells ',nShAt
    write(*,'(A28)') '------------------'

!   Basis function centers

    do j=1,nShAt
      nShell = nShell + 1
      do k=1,3
        CenterShell(nShell,k) = XYZAtoms(iAt,k)
      end do

! Shell type and contraction degree

      read(2,*) shelltype,KShell(nShell)
        if(shelltype == "S") then
          TotAngMomShell(nShell) = 0
          write(*,'(A28,1X,I16)') 's-type shell with K = ',KShell(nShell)
        elseif(shelltype == "P") then
          TotAngMomShell(nShell) = 1
          write(*,'(A28,1X,I16)') 'p-type shell with K = ',KShell(nShell)
        elseif(shelltype == "D") then
          TotAngMomShell(nShell) = 2
          write(*,'(A28,1X,I16)') 'd-type shell with K = ',KShell(nShell)
        elseif(shelltype == "F") then
          TotAngMomShell(nShell) = 3
          write(*,'(A28,1X,I16)') 'f-type shell with K = ',KShell(nShell)
        elseif(shelltype == "G") then
          TotAngMomShell(nShell) = 4
          write(*,'(A28,1X,I16)') 'g-type shell with K = ',KShell(nShell)
        elseif(shelltype == "H") then
          TotAngMomShell(nShell) = 5
          write(*,'(A28,1X,I16)') 'h-type shell with K = ',KShell(nShell)
        elseif(shelltype == "I") then
          TotAngMomShell(nShell) = 6
          write(*,'(A28,1X,I16)') 'i-type shell with K = ',KShell(nShell)
        end if

! Read exponents and contraction coefficients

        write(*,'(A28,1X,A16,A16)') '','Exponents','Contraction'
        do k=1,Kshell(nShell)
          read(2,*) ExpShell(nShell,k),DShell(nShell,k)
          write(*,'(A28,1X,F16.10,F16.10)') '',ExpShell(nShell,k),DShell(nShell,k)
        end do
    end do
    write(*,'(A28)') '------------------'
  end do

! Total number of shells

  write(*,'(A28,1X,I16)') 'Number of shells in OBS',nShell
  write(*,'(A28)') '------------------'
  write(*,*)

! Close file with basis set specification

  close(unit=2)

!------------------------------------------------------------------------
! Auxiliary basis set information
!------------------------------------------------------------------------

! Open file with auxilairy basis specification

  open(unit=3,file='input/auxbasis')

! Read basis information

  write(*,'(A28)') 'Auxiliary basis set'
  write(*,'(A28)') '-------------------'

  do i=1,NAtoms
    read(3,*) iAt,nShAt
    write(*,'(A28,1X,I16)') 'Atom n. ',iAt
    write(*,'(A28,1X,I16)') 'number of shells ',nShAt
    write(*,'(A28)') '------------------'

!   Basis function centers

    do j=1,nShAt
      nShell = nShell + 1
      do k=1,3
        CenterShell(nShell,k) = XYZAtoms(iAt,k)
      end do

! Shell type and contraction degree

      read(3,*) shelltype,KShell(nShell)
        if(shelltype == "S") then
          TotAngMomShell(nShell) = 0
          write(*,'(A28,1X,I16)') 's-type shell with K = ',KShell(nShell)
        elseif(shelltype == "P") then
          TotAngMomShell(nShell) = 1
          write(*,'(A28,1X,I16)') 'p-type shell with K = ',KShell(nShell)
        elseif(shelltype == "D") then
          TotAngMomShell(nShell) = 2
          write(*,'(A28,1X,I16)') 'd-type shell with K = ',KShell(nShell)
        elseif(shelltype == "F") then
          TotAngMomShell(nShell) = 3
          write(*,'(A28,1X,I16)') 'f-type shell with K = ',KShell(nShell)
        elseif(shelltype == "G") then
          TotAngMomShell(nShell) = 4
          write(*,'(A28,1X,I16)') 'g-type shell with K = ',KShell(nShell)
        elseif(shelltype == "H") then
          TotAngMomShell(nShell) = 5
          write(*,'(A28,1X,I16)') 'h-type shell with K = ',KShell(nShell)
        elseif(shelltype == "I") then
          TotAngMomShell(nShell) = 6
          write(*,'(A28,1X,I16)') 'i-type shell with K = ',KShell(nShell)
        end if

! Read exponents and contraction coefficients

        write(*,'(A28,1X,A16,A16)') '','Exponents','Contraction'
        do k=1,Kshell(nShell)
          read(3,*) ExpShell(nShell,k),DShell(nShell,k)
          write(*,'(A28,1X,F16.10,F16.10)') '',ExpShell(nShell,k),DShell(nShell,k)
        end do
    end do
    write(*,'(A28)') '------------------'
  end do

! Total number of shells

  write(*,'(A28,1X,I16)') 'Number of shells in ABS',nShell
  write(*,'(A28)') '------------------'
  write(*,*)

! Close file with basis set specification

  close(unit=3)

end subroutine 
