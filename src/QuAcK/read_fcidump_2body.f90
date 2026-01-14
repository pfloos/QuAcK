subroutine read_fcidump_2body(nBas,ERI_AO)

  implicit none

! Input variables

  integer,intent(in) :: nBas

! Output variables

  double precision,intent(out)  :: ERI_AO(nBas,nBas,nBas,nBas)

! Local variables

  logical                       :: file_exists
  integer                       :: ibas,jbas,kbas,lbas
  double precision              :: Val

  inquire(file='FCIDUMP', exist=file_exists)
  if(file_exists) then
   write(*,*)
   write(*,*) 'Reading FCIDUMP two-body integrals'
   write(*,*)
   ERI_AO=0d0 
   open(unit=314, form='formatted', file='FCIDUMP', status='old')
   do
    read(314,*) Val,ibas,jbas,kbas,lbas
    if(kbas==lbas .and. lbas==0) then
     if(ibas==jbas .and. ibas==0) then
      exit
     endif
    else
     ERI_AO(ibas,jbas,kbas,lbas)=Val
    endif
   enddo
  endif
  close(314)
 

end subroutine
