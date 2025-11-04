subroutine G_Dyson_GW_RHF(nBas,nOrb,nO,c,eHF,nfreqs,wweight,wcoord,ERI,vMAT,&
                         Enuc,EcGM,T,V,S,Pcorr)


! Use the restricted Sigma_c(E) to compute the linnearized approximation to G

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nfreqs
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nO

  double precision,intent(in)   :: Enuc
  double precision,intent(in)   :: EcGM
  double precision,intent(inout):: eHF(nOrb)
  double precision,intent(in)   :: wweight(nfreqs)
  double precision,intent(in)   :: wcoord(nfreqs)
  double precision,intent(in)   :: T(nBas,nBas)
  double precision,intent(in)   :: V(nBas,nBas)
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: c(nBas,nOrb)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: vMAT(nOrb*nOrb,nOrb*nOrb)

! Local variables

  integer                       :: kind_int
  integer                       :: isteps
  integer                       :: ifreq
  integer                       :: iorb,jorb
  integer                       :: nOrb2
  integer                       :: nfreqs2

  double precision              :: thrs_closer,chem_pot_change,grad_electrons,trace_2up,trace_2down
  double precision              :: trace_1rdm,trace_old,trace_up,trace_down,delta_chem_pot,thrs_N
  double precision              :: eta
  double precision              :: lim_inf,lim_sup
  double precision              :: alpha,beta
  double precision              :: ET,EV,EJ,EK,EgG,trace_occ
  double precision              :: chem_pot
  double precision,external     :: trace_matrix
  double precision,allocatable  :: J(:,:)
  double precision,allocatable  :: K(:,:)
  double precision,allocatable  :: Occ(:)
  double precision,allocatable  :: wweight2(:)
  double precision,allocatable  :: wcoord2(:)
  double precision,allocatable  :: Pcorr_mo(:,:)

  complex *16,allocatable       :: wcoord2_cpx(:)
  complex *16,allocatable       :: Sigma_c(:,:,:)
  complex *16,allocatable       :: Tmp_mo(:,:)
  complex *16,allocatable       :: Identity(:,:)

! Ouput variables

  double precision,intent(out)  :: Pcorr(nBas,nBas)

! Allocate and initialize arrays and variables

  write(*,*)
  write(*,*) '********************************************'
  write(*,*) '*          G = Go + Go Sigma G             *'
  write(*,*) '*          full-Dyson equation             *'
  write(*,*) '********************************************'
  write(*,*)

  eta=0d0
  nfreqs2=10*nfreqs
  allocate(Sigma_c(nfreqs2,nOrb,nOrb))
  allocate(Tmp_mo(nOrb,nOrb),Pcorr_mo(nOrb,nOrb),Identity(nOrb,nOrb))
  Identity=czero
  do iorb=1,nOrb
   Identity(iorb,iorb)=1d0
  enddo

! Prepare second quadrature

  kind_int = 1
  lim_inf = 0d0; lim_sup = 1d0;
  alpha = 0d0;   beta  = 0d0;
  allocate(wweight2(nfreqs2),wcoord2(nfreqs2),wcoord2_cpx(nfreqs2))
  call cgqf(nfreqs2,kind_int,alpha,beta,lim_inf,lim_sup,wcoord2,wweight2)
  wweight2(:)=wweight2(:)/((1d0-wcoord2(:))**2d0)
  wcoord2(:)=wcoord2(:)/(1d0-wcoord2(:))
  wcoord2_cpx(:)=wcoord2(:)*im

! Build Sigma_c(iw)

  call build_Sigmac_w_RHF(nOrb,nO,nfreqs2,eta,0,wcoord2_cpx,eHF,nfreqs,0,&
                          wweight,wcoord,vMAT,Sigma_c)

! Fix the chemical potential

  isteps = 0
  delta_chem_pot  = 2d-1
  thrs_closer     = 2d-1
  thrs_N          = 1d-8
  chem_pot_change = 0d0
  grad_electrons  = 1d0
  trace_1rdm      = -1d0
  chem_pot        = 0.1d0

  write(*,*)
  write(*,'(a)') ' Fixing the Tr[1D] at full-Dyson '
  write(*,*)
  write(*,*)'------------------------------------------------------'
  write(*,'(1X,A1,1X,A15,1X,A1,1X,A15,1X,A1A15,2X,A1)') &
          '|','Tr[1D]','|','Chem. Pot.','|','Grad N','|'
  write(*,*)'------------------------------------------------------'

  call trace_1rdm_Gdyson(nOrb,nO,nfreqs2,chem_pot,trace_old,eHF,wcoord2_cpx,wweight2,&
                         Tmp_mo,Sigma_c,Pcorr_mo,Identity)
  write(*,'(1X,A1,F16.10,1X,A1,F16.10,1X,A1F16.10,1X,A1)') &
  '|',trace_old,'|',chem_pot,'|',grad_electrons,'|'
  do while( abs(trace_old-2*nO) > thrs_closer .and. isteps <= 100 )
   isteps = isteps + 1
   call trace_1rdm_Gdyson(nOrb,nO,nfreqs2,chem_pot,trace_old,eHF,wcoord2_cpx,wweight2, &
                          Tmp_mo,Sigma_c,Pcorr_mo,Identity)
   call trace_1rdm_Gdyson(nOrb,nO,nfreqs2,chem_pot-delta_chem_pot,trace_down,eHF, &
                          wcoord2_cpx,wweight2,Tmp_mo,Sigma_c,Pcorr_mo,Identity)
   call trace_1rdm_Gdyson(nOrb,nO,nfreqs2,chem_pot+delta_chem_pot,trace_up,eHF, &
                          wcoord2_cpx,wweight2,Tmp_mo,Sigma_c,Pcorr_mo,Identity)
   if( abs(trace_up-2*nO) > abs(trace_old-2*nO) .and. abs(trace_down-2*nO) > abs(trace_old-2*nO) ) then
     write(*,'(1X,A1,F16.10,1X,A1,F16.10,1X,A1F16.10,1X,A1)') &
     '|',trace_old,'|',chem_pot,'|',grad_electrons,'|'
     delta_chem_pot = 0.75d0*delta_chem_pot
     thrs_closer = 0.5d0*thrs_closer
     write(*,*) "| contracting ...                                     |"
     if(delta_chem_pot<1d-2) exit
   else
     if( abs(trace_up-2*nO) < abs(trace_old-2*nO) ) then
      chem_pot=chem_pot+delta_chem_pot
      write(*,'(1X,A1,F16.10,1X,A1,F16.10,1X,A1F16.10,1X,A1)') &
      '|',trace_up,'|',chem_pot,'|',grad_electrons,'|'
     else
      if( abs(trace_down-2*nO) < abs(trace_old-2*nO) ) then
       chem_pot=chem_pot-delta_chem_pot
       write(*,'(1X,A1,F16.10,1X,A1,F16.10,1X,A1F16.10,1X,A1)') &
       '|',trace_down,'|',chem_pot,'|',grad_electrons,'|'
      endif
     endif
   endif
  enddo

  ! Do  final search

  write(*,*)'------------------------------------------------------'
  isteps = 0
  delta_chem_pot  = 1.0d-3
  do while( abs(trace_1rdm-2*nO) > thrs_N .and. isteps <= 100 )
   isteps = isteps + 1
   chem_pot = chem_pot + chem_pot_change
   call trace_1rdm_Gdyson(nOrb,nO,nfreqs2,chem_pot,trace_1rdm,eHF,wcoord2_cpx,wweight2,&
                          Tmp_mo,Sigma_c,Pcorr_mo,Identity)
   call trace_1rdm_Gdyson(nOrb,nO,nfreqs2,chem_pot+2d0*delta_chem_pot,trace_2up,eHF,   &
                          wcoord2_cpx,wweight2,Tmp_mo,Sigma_c,Pcorr_mo,Identity)
   call trace_1rdm_Gdyson(nOrb,nO,nfreqs2,chem_pot+delta_chem_pot,trace_up,eHF,        &
                          wcoord2_cpx,wweight2,Tmp_mo,Sigma_c,Pcorr_mo,Identity)
   call trace_1rdm_Gdyson(nOrb,nO,nfreqs2,chem_pot-delta_chem_pot,trace_down,eHF,      &
                          wcoord2_cpx,wweight2,Tmp_mo,Sigma_c,Pcorr_mo,Identity)
   call trace_1rdm_Gdyson(nOrb,nO,nfreqs2,chem_pot-2d0*delta_chem_pot,trace_2down,eHF, &
                          wcoord2_cpx,wweight2,Tmp_mo,Sigma_c,Pcorr_mo,Identity)
!   grad_electrons = (trace_up-trace_down)/(2d0*delta_chem_pot)
   grad_electrons = (-trace_2up+8d0*trace_up-8d0*trace_down+trace_2down)/(12d0*delta_chem_pot)
   chem_pot_change = -(trace_1rdm-2*nO)/(grad_electrons+1d-10)
   ! Maximum change is bounded within +/- 0.10
   chem_pot_change = max( min( chem_pot_change , 0.1d0 / real(isteps) ), -0.1d0 / real(isteps) )
   write(*,'(1X,A1,F16.10,1X,A1,F16.10,1X,A1F16.10,1X,A1)') &
   '|',trace_1rdm,'|',chem_pot,'|',grad_electrons,'|'
  enddo
  write(*,*)'------------------------------------------------------'
  write(*,*)
  call trace_1rdm_Gdyson(nOrb,nO,nfreqs2,chem_pot,trace_1rdm,eHF,wcoord2_cpx,wweight2,&
                         Tmp_mo,Sigma_c,Pcorr_mo,Identity)

! Compute  AO densities, new total energy, and Occ numbers

  Pcorr(:,:) = matmul(c,matmul(Pcorr_mo(:,:),transpose(c)))

  allocate(J(nBas,nBas))
  allocate(K(nBas,nBas))
  allocate(Occ(nOrb))

  call diagonalize_matrix(nOrb,Pcorr_mo,Occ)
  call Hartree_matrix_AO_basis(nBas,Pcorr,ERI,J)
  call exchange_matrix_AO_basis(nBas,Pcorr,ERI,K)

  ! Kinetic energy

  ET = trace_matrix(nBas,matmul(Pcorr,T))

  ! Potential energy

  EV = trace_matrix(nBas,matmul(Pcorr,V))

  ! Hartree energy

  EJ = 0.5d0*trace_matrix(nBas,matmul(Pcorr,J))

  ! Exchange energy

  EK = 0.25d0*trace_matrix(nBas,matmul(Pcorr,K))

  ! Total energy (incl. the Galitkii-Migdal contribution)

  EgG = ET + EV + EJ + EK + EcGM

! Print results

  write(*,*)
  write(*,'(A50)')           '---------------------------------------'
  write(*,'(A33)')           ' Summary               '
  write(*,'(A50)')           '---------------------------------------'
  write(*,'(A33,1X,F16.10,A3)') ' One-electron energy = ',ET + EV,' au'
  write(*,'(A33,1X,F16.10,A3)') ' Kinetic      energy = ',ET,' au'
  write(*,'(A33,1X,F16.10,A3)') ' Potential    energy = ',EV,' au'
  write(*,'(A50)')           '---------------------------------------'
  write(*,'(A33,1X,F16.10,A3)') ' Two-electron energy = ',EJ + EK + EcGM,' au'
  write(*,'(A33,1X,F16.10,A3)') ' Hartree      energy = ',EJ,' au'
  write(*,'(A33,1X,F16.10,A3)') ' Exchange     energy = ',EK,' au'
  write(*,'(A33,1X,F16.10,A3)') ' GM           energy = ',EcGM,' au'
  write(*,'(A50)')           '---------------------------------------'
  write(*,'(A33,1X,F16.10,A3)') ' Electronic   energy = ',EgG,' au'
  write(*,'(A33,1X,F16.10,A3)') ' Nuclear   repulsion = ',ENuc,' au'
  write(*,'(A33,1X,F16.10,A3)') ' G-Dyson    energy   = ',EgG + ENuc,' au'
  write(*,'(A50)')           '---------------------------------------'
  write(*,'(A33,1X,F16.10,A3)') ' Chemical potential  = ',chem_pot,' au'
  write(*,'(A50)')           '---------------------------------------'
  write(*,'(A50)') '---------------------------------------'
  write(*,'(A50)') ' G-Dyson occupation numbers '
  write(*,'(A50)') '---------------------------------------'
  trace_occ=0d0
  do iorb=1,nOrb
   if(abs(Occ(iorb))>1d-8) then
    write(*,'(I7,10F15.8)') iorb,Occ(nOrb-(iorb-1))
   endif
   trace_occ=trace_occ+Occ(iorb)
  enddo
  write(*,*)
  write(*,'(A33,1X,F16.10,A3)') ' Trace [ 1D ]        = ',trace_occ,'   '
  write(*,*)

  ! Deallocate arrays

  deallocate(J,K,Occ)
  deallocate(Sigma_c,Tmp_mo,Pcorr_mo,wweight2,wcoord2,wcoord2_cpx)

end subroutine
  
subroutine trace_1rdm_Gdyson(nOrb,nO,nfreqs2,chem_pot,trace_1rdm,eHF,wcoord2_cpx,wweight2, &
                             Tmp_mo,Sigma_c,Pcorr_mo,Identity)
  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nO
  integer,intent(in)            :: nfreqs2

  double precision,intent(in)   :: chem_pot
  double precision,intent(inout):: eHF(nOrb)
  double precision,intent(in)   :: wweight2(nfreqs2)

  complex *16,intent(in)        :: wcoord2_cpx(nfreqs2)
  complex *16,intent(in)        :: Identity(nOrb,nOrb)
  complex *16,intent(in)        :: Sigma_c(nfreqs2,nOrb,nOrb)

! Local variables

  integer                       :: ifreq
  integer                       :: iorb

  complex *16,allocatable       :: Go_mo(:,:)

! Output variables

  double precision,intent(inout):: trace_1rdm
  double precision,intent(inout):: Pcorr_mo(nOrb,nOrb)

  complex *16,intent(inout)     :: Tmp_mo(nOrb,nOrb)

! Integration along imag. freqs
  
  allocate(Go_mo(nOrb,nOrb))

  Pcorr_mo(:,:)=0d0
  do ifreq=1,nfreqs2
  
   call G_MO_RHF(nOrb,nO,0d0,eHF,wcoord2_cpx(ifreq),Go_mo)                 ! Go(iw2)
   Tmp_mo(:,:)=Go_mo(:,:)
   call complex_inverse_matrix(nOrb,Tmp_mo,Tmp_mo)                         ! Go(iw2)^-1
   Tmp_mo(:,:)=Tmp_mo(:,:) - Sigma_c(ifreq,:,:) - chem_pot*Identity(:,:)   ! G(iw2)^-1
   call complex_inverse_matrix(nOrb,Tmp_mo,Tmp_mo)                         ! G(iw2)
   Tmp_mo(:,:)=Tmp_mo(:,:) - Go_mo(:,:)                                    ! Gcorr(iw2) = G(iw2) - Go(iw2)
 
   Pcorr_mo(:,:) = Pcorr_mo(:,:) + 2d0*wweight2(ifreq)*real(Tmp_mo(:,:))   ! Integrate along iw2
  
  enddo
  Pcorr_mo(:,:) = Pcorr_mo(:,:)/(2d0*pi)
  
  ! P = Pcorr + P_HF with P_HF = delta_ij for the occ orbitals
  do iorb=1,nO
   Pcorr_mo(iorb,iorb) = Pcorr_mo(iorb,iorb) + 1d0
  enddo

  Pcorr_mo(:,:) = 2d0*Pcorr_mo(:,:)                                        ! Times 2 because we are adding both spin channels to build the spinless P
  
  trace_1rdm=0d0
  do iorb=1,nOrb
   trace_1rdm=trace_1rdm+Pcorr_mo(iorb,iorb)
  enddo

  deallocate(Go_mo)

end subroutine

