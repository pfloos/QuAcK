subroutine guess_read_rhf(nBas,P)

! Perform restricted Hartree-Fock calculation

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas

! Local variables

  logical                       :: file_exists
  integer                       :: ibas,jbas
  double precision              :: Val

! Output variables

  double precision,intent(out)  :: P(nBas,nBas)

  ! If guess_type is read density and files P_ao_bin and/or P_ao_form exists
   inquire(file='P_ao_bin', exist=file_exists)
   if(file_exists) then
    write(*,*) 'Reading P_ao_bin matrix...'
    open(unit=314, form='unformatted', file='P_ao_bin', status='old')
    do
     read(314) ibas,jbas,Val
     if(ibas==0 .and. jbas==0) exit
     P(ibas,jbas)=Val 
    enddo
    close(314)
   endif
   inquire(file='P_ao_form', exist=file_exists)
   if(file_exists) then
    write(*,*) 'Reading P_ao_form matrix...'
    open(unit=314, form='formatted', file='P_ao_form', status='old')
    ibas=1;jbas=0;
    do
     read(314,*) Val
     jbas=jbas+1
     if(jbas>nBas) then
      ibas=ibas+1
      jbas=1
     endif
     P(ibas,jbas)=Val 
     if(ibas==jbas .and. ibas==nBas) exit
    enddo
    close(314)
   endif

end subroutine 

