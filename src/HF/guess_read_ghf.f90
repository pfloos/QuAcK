subroutine guess_read_ghf(nBas,P)

! Perform unrestricted Hartree-Fock calculation

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas

! Local variables


  logical                       :: file_exists
  integer                       :: ibas,jbas
  double precision              :: Val

! Output variables

  double precision,intent(out)  :: P(nBas*2,nBas*2)

   ! Read P_aa
   inquire(file='P_ao_form_aa', exist=file_exists)
   if(file_exists) then
    write(*,*) 'Reading P_ao_form matrices...'
    open(unit=314, form='formatted', file='P_ao_form_aa', status='old')
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
   ! Read P_bb
   inquire(file='P_ao_form_bb', exist=file_exists)
   if(file_exists) then
    open(unit=314, form='formatted', file='P_ao_form_bb', status='old')
    ibas=1;jbas=0;
    do
     read(314,*) Val
     jbas=jbas+1
     if(jbas>nBas) then
      ibas=ibas+1
      jbas=1
     endif
     P(nBas+ibas,nBas+jbas)=Val 
     if(ibas==jbas .and. ibas==nBas) exit
    enddo
    close(314)
   endif
   ! Read P_ab
   inquire(file='P_ao_form_ab', exist=file_exists)
   if(file_exists) then
    open(unit=314, form='formatted', file='P_ao_form_ab', status='old')
    ibas=1;jbas=0;
    do
     read(314,*) Val
     jbas=jbas+1
     if(jbas>nBas) then
      ibas=ibas+1
      jbas=1
     endif
     P(ibas,nBas+jbas)=Val 
     if(ibas==jbas .and. ibas==nBas) exit
    enddo
    close(314)
   endif
   ! Read P_ba
   inquire(file='P_ao_form_ba', exist=file_exists)
   if(file_exists) then
    open(unit=314, form='formatted', file='P_ao_form_ba', status='old')
    ibas=1;jbas=0;
    do
     read(314,*) Val
     jbas=jbas+1
     if(jbas>nBas) then
      ibas=ibas+1
      jbas=1
     endif
     P(ibas+nBas,jbas)=Val 
     if(ibas==jbas .and. ibas==nBas) exit
    enddo
    close(314)
   endif

end subroutine
