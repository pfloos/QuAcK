
subroutine print_RHFB_AO_rdms(nBas,ENuc,sigma,S,T,V,P,Panom,ERI)

  implicit none

  integer,intent(in)            :: nBas
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: sigma
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: T(nBas,nBas)
  double precision,intent(in)   :: V(nBas,nBas)
  double precision,intent(in)   :: P(nBas,nBas)
  double precision,intent(in)   :: Panom(nBas,nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

  integer                       :: iunit,iunit2
  integer                       :: ibas,jbas,kbas,lbas
  double precision,allocatable  :: AO_1rdm(:,:)
  double precision,allocatable  :: AO_2rdm(:,:,:,:)
  double precision              :: trace_1rdm
  double precision              :: trace_2rdm
  double precision              :: S2_val
  double precision              :: Ecore
  double precision              :: Eee

  write(*,*)
  write(*,'(a)') ' -------------------------------------------'
  write(*,'(a)') ' Computing and printing RDMs in the AO basis'
  write(*,'(a)') ' -------------------------------------------'
  allocate(AO_1rdm(nBas,nBas),AO_2rdm(nBas,nBas,nBas,nBas))
  AO_1rdm(:,:)=P(:,:)
  AO_2rdm(:,:,:,:)=0d0
  do ibas=1,nBas
   do jbas=1,nBas
    do kbas=1,nBas
     do lbas=1,nBas
      ! Hartree
      AO_2rdm(ibas,jbas,kbas,lbas)=AO_2rdm(ibas,jbas,kbas,lbas)+0.5d0*P(ibas,kbas)*P(jbas,lbas)
      ! Exchange
      AO_2rdm(ibas,jbas,kbas,lbas)=AO_2rdm(ibas,jbas,kbas,lbas)-0.25d0*P(ibas,lbas)*P(jbas,kbas)
      ! Pairing
      AO_2rdm(ibas,jbas,kbas,lbas)=AO_2rdm(ibas,jbas,kbas,lbas)+sigma*Panom(ibas,jbas)*Panom(kbas,lbas)
      !AO_2rdm(ibas,jbas,kbas,lbas)=AO_2rdm(ibas,jbas,kbas,lbas)-Panom(ibas,jbas)*Panom(kbas,lbas) ! This is CA NOFA with correct trace
     enddo
    enddo
   enddo
  enddo
  trace_1rdm=0d0
  trace_2rdm=0d0
  iunit=312
  iunit2=313
  open(unit=iunit,form='unformatted',file='rhfb_ao_1rdm')
  open(unit=iunit2,form='unformatted',file='rhfb_ao_2rdm')
  Ecore=0d0; Eee=0d0;
  do ibas=1,nBas
   do jbas=1,nBas
    trace_1rdm=trace_1rdm+AO_1rdm(ibas,jbas)*S(ibas,jbas)
    write(iunit) ibas,jbas,AO_1rdm(ibas,jbas)
    Ecore=Ecore+AO_1rdm(ibas,jbas)*(T(ibas,jbas)+V(ibas,jbas))
    do kbas=1,nBas
     do lbas=1,nBas
      trace_2rdm=trace_2rdm+AO_2rdm(ibas,jbas,kbas,lbas)*S(ibas,kbas)*S(jbas,lbas)
      write(iunit2) ibas,jbas,kbas,lbas,AO_2rdm(ibas,jbas,kbas,lbas)
      Eee=Eee+AO_2rdm(ibas,jbas,kbas,lbas)*ERI(ibas,jbas,kbas,lbas)
     enddo
    enddo
   enddo
  enddo
  write(iunit) 0,0,0d0
  write(iunit2) 0,0,0,0,0d0
  close(iunit)
  close(iunit2)
  ! Compute <S^2>
  S2_val=0d0
   ! Density contribution
   S2_val=-trace_1rdm*(trace_1rdm-4d0)/4d0
   ! Daa = Dbb
   AO_2rdm=0d0
   do ibas=1,nBas
    do jbas=1,nBas
     do kbas=1,nBas
      do lbas=1,nBas
       ! Hartree
       AO_2rdm(ibas,jbas,kbas,lbas)=AO_2rdm(ibas,jbas,kbas,lbas)+0.125d0*P(ibas,kbas)*P(jbas,lbas)
       ! Exchange
       AO_2rdm(ibas,jbas,kbas,lbas)=AO_2rdm(ibas,jbas,kbas,lbas)-0.125d0*P(ibas,lbas)*P(jbas,kbas)
       ! Contribution to <S^2>
       S2_val=S2_val+2d0*AO_2rdm(ibas,jbas,kbas,lbas)*S(ibas,kbas)*S(jbas,lbas)
      enddo
     enddo
    enddo
   enddo
   ! Dab = Dba
   AO_2rdm=0d0
   do ibas=1,nBas
    do jbas=1,nBas
     do kbas=1,nBas
      do lbas=1,nBas
       ! Hartree
       AO_2rdm(ibas,jbas,kbas,lbas)=AO_2rdm(ibas,jbas,kbas,lbas)+0.125d0*P(ibas,kbas)*P(jbas,lbas)
       ! Pairing
       AO_2rdm(ibas,jbas,kbas,lbas)=AO_2rdm(ibas,jbas,kbas,lbas)+0.5d0*sigma*Panom(ibas,jbas)*Panom(kbas,lbas)
       ! Contribution to <S^2>
       S2_val=S2_val-2d0*AO_2rdm(ibas,jbas,kbas,lbas)*S(ibas,lbas)*S(jbas,kbas)
      enddo
     enddo
    enddo
   enddo
  deallocate(AO_1rdm,AO_2rdm)
  write(*,'(a)') '  Energies computed using the 1-RDM and the 2-RDM in the AO basis'
  write(*,*)
  write(*,'(a,f17.8)') '   Hcore (T+V) ',Ecore
  write(*,'(a,f17.8)') '     Vee (Hxc) ',Eee
  write(*,'(a,f17.8)') '   Eelectronic ',Ecore+Eee
  write(*,*)           ' --------------'
  write(*,'(a,f17.8)') '        Etotal ',Ecore+Eee+ENuc
  write(*,*)           ' --------------'
  write(*,'(a,f17.8)') '   Tr[ 1D^AO ] ',trace_1rdm
  write(*,'(a,f17.8)') '   Tr[ 2D^AO ] ',trace_2rdm
  write(*,'(a,f17.8)') '         <S^2> ',S2_val
  write(*,*)


end subroutine
