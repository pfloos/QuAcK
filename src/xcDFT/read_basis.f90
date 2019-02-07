subroutine read_basis(nAt,rAt,nBas,nO,nV,nShell,TotAngMomShell,CenterShell,KShell,DShell,ExpShell)

! Read basis set information

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nAt,nO
  double precision,intent(in)   :: rAt(nAt,3)

! Local variables

  integer                       :: nShAt,iAt,iShell
  integer                       :: i,j,k
  character                     :: shelltype

! Output variables

  integer,intent(out)           :: nShell,nBas,nV
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
  do i=1,nAt
    read(2,*) iAt,nShAt
    write(*,'(A28,1X,I16)') 'Atom n. ',iAt
    write(*,'(A28,1X,I16)') 'number of shells ',nShAt
    write(*,'(A28)') '------------------'

!   Basis function centers

    do j=1,nShAt
      nShell = nShell + 1
      do k=1,3
        CenterShell(nShell,k) = rAt(iAt,k)
      enddo

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
        endif

! Read exponents and contraction coefficients

        write(*,'(A28,1X,A16,A16)') '','Exponents','Contraction'
        do k=1,Kshell(nShell)
          read(2,*) ExpShell(nShell,k),DShell(nShell,k)
          write(*,'(A28,1X,F16.10,F16.10)') '',ExpShell(nShell,k),DShell(nShell,k)
        enddo
    enddo
    write(*,'(A28)') '------------------'
  enddo

! Total number of shells

  write(*,'(A28,1X,I16)') 'Number of shells',nShell
  write(*,'(A28)') '------------------'
  write(*,*)

! Close file with basis set specification

  close(unit=2)

! Calculate number of basis functions

  nBas = 0
  do iShell=1,nShell
    nBas = nBas + (TotAngMomShell(iShell)*TotAngMomShell(iShell) + 3*TotAngMomShell(iShell) + 2)/2
  enddo

  write(*,'(A28)') '------------------'
  write(*,'(A28,1X,I16)') 'Number of basis functions',NBas
  write(*,'(A28)') '------------------'
  write(*,*)

! Number of virtual orbitals

  nV = nBas - nO 

end subroutine read_basis
