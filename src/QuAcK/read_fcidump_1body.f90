subroutine read_fcidump_1body(nBas,nOrb,ncart,S,T,V,Hc,X,dipole_int_AO,ENuc,Znuc)

  implicit none

! Input variables

  integer,intent(in) :: nBas,nOrb,ncart

! Output variables

  double precision,intent(out)  :: ENuc,Znuc
  double precision,intent(out)  :: S(nBas,nBas)
  double precision,intent(out)  :: X(nBas,nOrb)
  double precision,intent(out)  :: T(nBas,nBas)
  double precision,intent(out)  :: V(nBas,nBas)
  double precision,intent(out)  :: Hc(nBas,nBas)
  double precision,intent(out)  :: dipole_int_AO(nBas,nBas,ncart)

! Local variables

  logical                       :: file_exists
  integer                       :: ibas,jbas,kbas,lbas
  double precision              :: Val
 
 inquire(file='FCIDUMP', exist=file_exists)
 if(file_exists) then
  write(*,*)
  write(*,*) 'Reading FCIDUMP one-body integrals'
  write(*,*)
  S=0d0; T=0d0; V=0d0; Hc=0d0; X=0d0;
  dipole_int_AO=0d0; ENuc=0d0; Znuc=0d0;
  do ibas=1,nBas
   S(ibas,ibas) = 1d0
   X(ibas,ibas) = 1d0
  enddo
  open(unit=314, form='formatted', file='FCIDUMP', status='old')
  do
   read(314,*) Val,ibas,jbas,kbas,lbas
   if(kbas==lbas .and. lbas==0) then
    if(ibas==jbas .and. ibas==0) then
     ENuc=Val
     exit
    else
     T(ibas,jbas) =Val
     T(jbas,ibas) =T(ibas,jbas)
     Hc(ibas,jbas)=T(ibas,jbas)
     Hc(jbas,ibas)=T(ibas,jbas)
    endif
   endif
  enddo
  close(314)
 endif

end subroutine
