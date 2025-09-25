
subroutine EcRPA_Bogoliubov_FCIDUMP(nO,nOrb,nOrb_twice,sigma,ERI_MO,vMAT,Fock,Delta,pMAT,panomMAT,eQP_state,U_QP, &
                                    chem_pot,ntimes,nfreqs,wcoord,wweight)
  implicit none
  include 'parameters.h'
!
  integer,intent(in)             :: nO
  integer,intent(in)             :: nOrb
  integer,intent(in)             :: nOrb_twice
  integer,intent(in)             :: ntimes
  integer,intent(in)             :: nfreqs
  double precision,intent(in)    :: sigma
  double precision,intent(in)    :: wcoord(nfreqs),wweight(nfreqs)

! Local variables

  logical                        :: file_exists
  integer                        :: iorb,jorb,korb,lorb
  double precision               :: Ecore,ENuc,EJ,EK,EL,Eelec,EHFB,thrs_N,trace_1rdm,Val,EcRPA,EcGM  
  double precision,external      :: trace_matrix
  double precision,allocatable   :: J(:,:)
  double precision,allocatable   :: K(:,:)
  double precision,allocatable   :: H_HFB(:,:)
  double precision,allocatable   :: Hc(:,:)
  double precision,allocatable   :: R(:,:)

! Output variables

  double precision,intent(inout) :: chem_pot
  double precision,intent(out)   :: ERI_MO(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(out)   :: vMAT(nOrb*nOrb,nOrb*nOrb)
  double precision,intent(out)   :: Fock(nOrb,nOrb)
  double precision,intent(out)   :: Delta(nOrb,nOrb)
  double precision,intent(out)   :: pMAT(nOrb,nOrb)
  double precision,intent(out)   :: panomMAT(nOrb,nOrb)
  double precision,intent(out)   :: eQP_state(nOrb_twice)
  double precision,intent(out)   :: U_QP(nOrb_twice,nOrb_twice)

  character(len=100)             :: arg                              
!

  write(*,*)
  write(*,*) '*********************************'
  write(*,*) '* Bogoliubov EcRPA from FCIDUMP *'
  write(*,*) '*********************************'
  write(*,*)

  ENuc=0d0
  allocate(J(nOrb,nOrb))
  allocate(K(nOrb,nOrb))
  allocate(Hc(nOrb,nOrb))
  allocate(H_HFB(nOrb_twice,nOrb_twice))
  allocate(R(nOrb_twice,nOrb_twice))
  
  inquire(file='FCIDUMP', exist=file_exists)
  if(file_exists) then
   write(*,*) 'Reading FCIDUMP file before computing EcRPA Bogoliubov'
   ERI_MO=0d0; vMAT=0d0; ENuc=0d0; Hc=0d0;
   open(unit=314, form='formatted', file='FCIDUMP', status='old')
   read(314,'(a)') arg
   read(314,'(a)') arg
   read(314,'(a)') arg
   do
    read(314,*) Val,iorb,jorb,korb,lorb
    if(iorb==jorb .and. jorb==korb .and. korb==lorb .and. iorb==0) then
     ENuc=Val
     exit
    endif
    if(korb==lorb .and. lorb==0) then
     Hc(iorb,jorb)=Val
     Hc(jorb,iorb)=Hc(iorb,jorb)
    else
     ERI_MO(iorb,korb,jorb,lorb)=Val
     ERI_MO(iorb,lorb,jorb,korb)=Val
     ERI_MO(jorb,lorb,iorb,korb)=Val
     ERI_MO(jorb,korb,iorb,lorb)=Val
     ERI_MO(korb,iorb,lorb,jorb)=Val
     ERI_MO(lorb,iorb,korb,jorb)=Val
     ERI_MO(lorb,jorb,korb,iorb)=Val
     ERI_MO(korb,jorb,lorb,iorb)=Val
    endif
   enddo
   close(314)
  endif
  do iorb=1,nOrb
   do jorb=1,nOrb
    do korb=1,nOrb
     do lorb=1,nOrb
      vMAT(1+(korb-1)+(iorb-1)*nOrb,1+(lorb-1)+(jorb-1)*nOrb)=ERI_MO(iorb,jorb,korb,lorb)
     enddo
    enddo
   enddo
  enddo

   inquire(file='occupancies', exist=file_exists)
   if(file_exists) then
    write(*,*) 'Reading occupancies file before computing EcRPA Bogoliubov'
    pMAT=0d0
    panomMAT=0d0
    open(unit=314, form='formatted', file='occupancies', status='old')
    do iorb=1,nOrb
      read(314,*) Val
      pMAT(iorb,iorb)=Val
      Val=0.5d0*Val
      panomMAT(iorb,iorb)=sqrt(abs(Val*(1.0d0-Val)))
    enddo
    close(314)
   endif
   
   call Hartree_matrix_AO_basis(nOrb,pMAT,ERI_MO,J)
   call exchange_matrix_AO_basis(nOrb,pMAT,ERI_MO,K)
   call anomalous_matrix_AO_basis(nOrb,sigma,panomMAT,ERI_MO,Delta)
   thrs_N=1d-6
   Fock(:,:)=Hc(:,:)+J(:,:)+0.5d0*K(:,:)
   do iorb=1,nOrb
   Fock(iorb,iorb)=Fock(iorb,iorb)-chem_pot
   enddo
   ! Hcore energy
   Ecore = trace_matrix(nOrb,matmul(pMAT,Hc))
   ! Coulomb energy
   EJ = 0.5d0*trace_matrix(nOrb,matmul(pMAT,J))
   ! Exchange energy
   EK = 0.25d0*trace_matrix(nOrb,matmul(pMAT,K))
   ! Anomalous energy
   EL = trace_matrix(nOrb,matmul(panomMAT,Delta))
   ! Total electronic energy
   Eelec = Ecore + EJ + EK + EL
   ! Total energy
   EHFB = Eelec + ENuc
   H_HFB(1:nOrb,1:nOrb)=Fock(1:nOrb,1:nOrb)
   H_HFB(nOrb+1:nOrb_twice,nOrb+1:nOrb_twice)=-Fock(1:nOrb,1:nOrb)
   H_HFB(1:nOrb,nOrb+1:nOrb_twice)=Delta(1:nOrb,1:nOrb)
   H_HFB(nOrb+1:nOrb_twice,1:nOrb)=Delta(1:nOrb,1:nOrb)
   call diagonalize_matrix(nOrb_twice,U_QP,eQP_state)
   ! Build R 
   R(:,:)     = 0d0
   do iorb=1,nOrb
    R(:,:) = R(:,:) + matmul(U_QP(:,iorb:iorb),transpose(U_QP(:,iorb:iorb)))
   enddo
   trace_1rdm = 0d0
   do iorb=1,nOrb
    trace_1rdm = trace_1rdm+R(iorb,iorb)
   enddo
   ! Adjust the chemical potential 
   if( abs(trace_1rdm-nO) > thrs_N ) &
    call fix_chem_pot(nO,nOrb,nOrb_twice,0,thrs_N,trace_1rdm,chem_pot,H_HFB,U_QP,R,eQP_state)
 
   write(*,*)
   write(*,'(A50)')           '---------------------------------------'
   write(*,'(A33)')           ' Summary               '
   write(*,'(A50)')           '---------------------------------------'
   write(*,'(A33,1X,F16.10,A3)') ' One-electron energy = ',Ecore,' au'
   write(*,'(A50)')           '---------------------------------------'
   write(*,'(A33,1X,F16.10,A3)') ' Two-electron energy = ',EJ + EK + EL,' au'
   write(*,'(A33,1X,F16.10,A3)') ' Hartree      energy = ',EJ,' au'
   write(*,'(A33,1X,F16.10,A3)') ' Exchange     energy = ',EK,' au'
   write(*,'(A33,1X,F16.10,A3)') ' Anomalous    energy = ',EL,' au'
   write(*,'(A50)')           '---------------------------------------'
   write(*,'(A33,1X,F16.10,A3)') ' Electronic   energy = ',Eelec,' au'
   write(*,'(A33,1X,F16.10,A3)') ' Nuclear   repulsion = ',ENuc,' au'
   write(*,'(A33,1X,F16.10,A3)') ' HFB          energy = ',EHFB,' au'
   write(*,'(A50)')           '---------------------------------------'
   write(*,'(A33,1X,F16.10,A3)') ' Chemical potential  = ',chem_pot,' au'
   write(*,'(A50)')           '---------------------------------------'
   write(*,*)
   write(*,'(A50)') '---------------------------------------'
   write(*,'(A50)') ' HFB QP energies '
   write(*,'(A50)') '---------------------------------------'
   do iorb=1,nOrb_twice
    write(*,'(I7,10F15.8)') iorb,eQP_state(iorb)
   enddo
   write(*,*)
   
   !U_QP(:,1:nOrb)=0d0
   !do iorb=1,nOrb
   ! U_QP(iorb,iorb) = sqrt(abs(0.5d0*pMAT(iorb,iorb)))
   ! U_QP(iorb+nOrb,iorb) = sqrt(abs(1.0d0-0.5d0*pMAT(iorb,iorb)))
   !enddo
   Eelec = 0d0
   call EcRPA_EcGM_w_RHFB(nOrb,nOrb_twice,1,eQP_state,nfreqs,ntimes,wweight,wcoord,vMAT, &
                          U_QP,EHFB,EcRPA,EcGM)
   deallocate(J,K,Hc)
   deallocate(H_HFB,R)

end subroutine

