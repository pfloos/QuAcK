
! ---

subroutine read_restart_HFB(nBas, nOrb, Occ, c, S, chem_pot)

! Write the binary file used to restart calcs

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)                :: nBas
  integer,intent(in)                :: nOrb
  double precision,intent(in)       :: S(nBas,nBas)

! Local variables

  integer                           :: ibas,iorb,iorb1,iunit=667

  integer                           :: nBas_
  integer                           :: nOrb_
  double precision                  :: chem_pot_
  double precision                  :: max_diff
  double precision,allocatable      :: eigVAL(:)
  double precision,allocatable      :: c_tmp(:,:)
  double precision,allocatable      :: S_mol(:,:)
  double precision,allocatable      :: X_mol(:,:)

! Output variables

  double precision,intent(out)      :: chem_pot
  double precision,intent(out)      :: Occ(nOrb)
  double precision,intent(out)      :: c(nBas,nOrb)

! Dump results

  allocate(eigVAL(nOrb),c_tmp(nBas,nOrb),S_mol(nOrb,nOrb),X_mol(nOrb,nOrb))

  open(unit=iunit,form='unformatted',file='hfb_bin',status='old')
  read(iunit) nBas_,nOrb_ 
  read(iunit) chem_pot_
  do iorb=1,nOrb 
   read(iunit) c_tmp(:,iorb) 
  enddo
  do iorb=1,nOrb 
   read(iunit) eigVAL(iorb) 
  enddo
  close(iunit)

  if(nBas==nBas_ .and. nOrb==nOrb_) then
   write(*,*)
   write(*,*)' Reading restart file'
   write(*,*)

   chem_pot = chem_pot_
   Occ(:) = eigVAL(:)
   c(:,:) = c_tmp(:,:)
   
   ! Check the orthonormality

   max_diff=-1d0
   S_mol = matmul(transpose(c),matmul(S,c))
   do iorb=1,nOrb
    do iorb1=1,iorb
     if(iorb==iorb1) then
      if(abs(S_mol(iorb,iorb)-1d0)>1d-8) max_diff=abs(S_mol(iorb,iorb)-1d0)
     else
      if(abs(S_mol(iorb,iorb1))>1d-8) max_diff=abs(S_mol(iorb,iorb1))
     endif
    enddo
   enddo
   if(max_diff>1d-8) then 
    write(*,*) ' '
    write(*,'(a)') ' Orthonormality violations. Applying Lowdin orthonormalization.'
    write(*,*) ' '
    X_mol = S_mol
    S_mol = 0d0
    call diagonalize_matrix(nOrb,X_mol,eigVAL)
    do iorb=1,nOrb
     S_mol(iorb,iorb) = 1d0/(sqrt(eigVAL(iorb)) + 1d-10)
    enddo
    X_mol = matmul(X_mol,matmul(S_mol,transpose(X_mol)))
    c = matmul(c,X_mol)
    max_diff=-1d0
    S_mol = matmul(transpose(c),matmul(S,c))
    do iorb=1,nOrb
     do iorb1=1,iorb
      if(iorb==iorb1) then
       if(abs(S_mol(iorb,iorb)-1d0)>1d-8) max_diff=abs(S_mol(iorb,iorb)-1d0)
      else
       if(abs(S_mol(iorb,iorb1))>1d-8) max_diff=abs(S_mol(iorb,iorb1))
      endif
     enddo
    enddo
    if(max_diff<=1d-8) then 
     write(*,*) ' '
     write(*,'(a)') ' No orthonormality violations.'
     write(*,*) ' '
    else 
     write(*,*) ' '
     write(*,'(a)') ' Error in Lowdin orthonormalization.'
     write(*,*) ' '
     stop
    endif
   else
    write(*,*) ' '
    write(*,'(a)') ' No orthonormality violations.'
    write(*,*) ' '
   endif
 
  endif  

  deallocate(eigVAL,c_tmp,S_mol,X_mol)

end subroutine 
