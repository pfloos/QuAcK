subroutine linDyson_GW_RHFB(nBas,nOrb,nOrb_twice,c,eQP_state,nfreqs,wweight,wcoord,ERI,vMAT,U_QP,&
                            Enuc,EcGM,sigma,T,V,S,X,P,Panom,Pcorr,Panomcorr)

! Use the restricted Sigma_c(E) to compute the linearized approx. to the Dyson eq

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nfreqs
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nOrb_twice

  double precision,intent(in)   :: Enuc,EcGM,sigma
  double precision,intent(inout):: eQP_state(nOrb_twice)
  double precision,intent(in)   :: wweight(nfreqs)
  double precision,intent(in)   :: wcoord(nfreqs)
  double precision,intent(in)   :: T(nBas,nBas)
  double precision,intent(in)   :: V(nBas,nBas)
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: X(nBas,nBas)
  double precision,intent(in)   :: c(nBas,nOrb)
  double precision,intent(in)   :: P(nBas,nBas)
  double precision,intent(in)   :: Panom(nBas,nBas)
  double precision,intent(in)   :: U_QP(nOrb_twice,nOrb_twice)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: vMAT(nOrb*nOrb,nOrb*nOrb)

! Local variables

  logical                       :: file_exists

  integer                       :: kind_int
  integer                       :: ifreq
  integer                       :: iunit,iunit2
  integer                       :: iorb,jorb
  integer                       :: ibas,jbas,kbas,lbas
  integer                       :: nfreqs2
  integer                       :: nBas_twice

  double precision              :: trace_1rdm
  double precision              :: trace_2rdm
  double precision              :: Ecore
  double precision              :: Eee
  double precision              :: eta
  double precision              :: lim_inf,lim_sup
  double precision              :: alpha,beta
  double precision              :: ET,EV,EJ,EK,EL,ElinG,trace_occ,N_anom,trace_r_can,dev_Idemp_r_can
  double precision,external     :: trace_matrix
  double precision,allocatable  :: wweight2(:)
  double precision,allocatable  :: wcoord2(:)
  double precision,allocatable  :: Occ(:)
  double precision,allocatable  :: Occ_R(:)
  double precision,allocatable  :: Mat1(:,:)
  double precision,allocatable  :: Mat2(:,:)
  double precision,allocatable  :: Rcorr(:,:)
  double precision,allocatable  :: J(:,:)
  double precision,allocatable  :: K(:,:)
  double precision,allocatable  :: Delta(:,:)
  double precision,allocatable  :: c_inv(:,:)
  double precision,allocatable  :: Pcorr_mo(:,:)
  double precision,allocatable  :: Panomcorr_mo(:,:)
  double precision,allocatable  :: R_ao(:,:)
  double precision,allocatable  :: AO_1rdm(:,:)
  double precision,allocatable  :: AO_2rdm(:,:,:,:)

  complex *16,allocatable       :: Sigma_c_he(:,:,:)
  complex *16,allocatable       :: Sigma_c_hh(:,:,:)
  complex *16,allocatable       :: Sigma_c_eh(:,:,:)
  complex *16,allocatable       :: Sigma_c_ee(:,:,:)
  complex *16,allocatable       :: Tmp_mo(:,:)
  complex *16,allocatable       :: Tmp_QP(:,:)
  complex *16,allocatable       :: Sigma_c_QP(:,:)
  complex *16,allocatable       :: wcoord2_cpx(:)

! Ouput variables

  double precision,intent(out)   :: Pcorr(nBas,nBas)
  double precision,intent(out)   :: Panomcorr(nBas,nBas)

! Allocate and initialize arrays and variables

  write(*,*)
  write(*,*) '******************************************************************'
  write(*,*) '*    G^Gorkov = Go^Gorkov + Go^Gorkov Sigma^Gorkov Go^Gorkov     *'
  write(*,*) '*       Bogoliubov linearized-Dyson equation approximation       *'
  write(*,*) '******************************************************************'
  write(*,*)

  eta=0d0
  nfreqs2=10*nfreqs
  nBas_twice=2*nBas
  allocate(Tmp_mo(nOrb,nOrb),Pcorr_mo(nOrb,nOrb),Panomcorr_mo(nOrb,nOrb),R_ao(nBas_twice,nBas_twice))
  allocate(Tmp_QP(nOrb_twice,nOrb_twice))
  allocate(Sigma_c_QP(nOrb_twice,nOrb_twice))
  allocate(Rcorr(nOrb_twice,nOrb_twice))
  allocate(Sigma_c_he(nfreqs2,nOrb,nOrb))
  allocate(Sigma_c_hh(nfreqs2,nOrb,nOrb))
  allocate(Sigma_c_eh(nfreqs2,nOrb,nOrb))
  allocate(Sigma_c_ee(nfreqs2,nOrb,nOrb))
  allocate(Mat1(nOrb,nOrb),Mat2(nOrb,nOrb))
  Mat1(1:nOrb,1:nOrb)=U_QP(1:nOrb,1:nOrb)
  Mat2(1:nOrb,1:nOrb)=U_QP(nOrb+1:nOrb_twice,1:nOrb)

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

  call build_Sigmac_w_RHFB(nOrb,nOrb_twice,nfreqs2,eta,0,wcoord2_cpx,eQP_state,nfreqs,0,wweight,wcoord, & 
                           vMAT,U_QP,Sigma_c_he,Sigma_c_hh,Sigma_c_eh,Sigma_c_ee,.true.,.true.)

! Integration along imag. freqs contributions

  Rcorr(:,:)=0d0
  do ifreq=1,nfreqs2
  
   ! Sigma_c^Gorkov
   Sigma_c_QP(1:nOrb           ,1:nOrb           ) =  Sigma_c_he(ifreq,1:nOrb,1:nOrb)
   Sigma_c_QP(1:nOrb           ,nOrb+1:nOrb_twice) = -Sigma_c_hh(ifreq,1:nOrb,1:nOrb)
   Sigma_c_QP(nOrb+1:nOrb_twice,1:nOrb           ) = -Sigma_c_ee(ifreq,1:nOrb,1:nOrb)
   Sigma_c_QP(nOrb+1:nOrb_twice,nOrb+1:nOrb_twice) =  Sigma_c_eh(ifreq,1:nOrb,1:nOrb)

   ! G^Gorkov
   call G_MO_RHFB(nOrb,nOrb_twice,eta,eQP_state,wcoord2_cpx(ifreq),Mat1,Mat1,Mat2, Mat2,Tmp_mo) ! G_he(iw2)
   Tmp_QP(1:nOrb           ,1:nOrb           ) = Tmp_mo(1:nOrb,1:nOrb)
   call G_MO_RHFB(nOrb,nOrb_twice,eta,eQP_state,wcoord2_cpx(ifreq),Mat1,Mat2,-Mat2,Mat1,Tmp_mo) ! G_hh(iw2)
   Tmp_QP(1:nOrb           ,nOrb+1:nOrb_twice) = Tmp_mo(1:nOrb,1:nOrb)
   call G_MO_RHFB(nOrb,nOrb_twice,eta,eQP_state,wcoord2_cpx(ifreq),Mat2,Mat1,Mat1,-Mat2,Tmp_mo) ! G_ee(iw2)
   Tmp_QP(nOrb+1:nOrb_twice,1:nOrb           ) = Tmp_mo(1:nOrb,1:nOrb)
   call G_MO_RHFB(nOrb,nOrb_twice,eta,eQP_state,wcoord2_cpx(ifreq),Mat2,Mat2,Mat1, Mat1,Tmp_mo) ! G_eh(iw2)
   Tmp_QP(nOrb+1:nOrb_twice,nOrb+1:nOrb_twice) = Tmp_mo(1:nOrb,1:nOrb)

   Tmp_QP(:,:)=matmul(Tmp_QP(:,:),matmul(Sigma_c_QP(:,:),Tmp_QP(:,:)))  ! This is G^Gorkov(iw2) Sigma_c^Gorkov(iw2) G^Gorkov(iw2)

   Rcorr(:,:) = Rcorr(:,:) + wweight2(ifreq)*real(Tmp_QP(:,:)+conjg(Tmp_QP(:,:))) ! Integrate along iw2 

  enddo

  Rcorr(:,:)     = Rcorr(:,:)/(2d0*pi)
  Pcorr(:,:)     = P(:,:)     + 2d0*matmul(c,matmul(Rcorr(1:nOrb,1:nOrb),transpose(c)))
  Panomcorr(:,:) = Panom(:,:) + matmul(c,matmul(Rcorr(1:nOrb,nOrb+1:nOrb_twice),transpose(c)))

  R_ao(1:nBas           ,1:nBas           )=0.5d0*Pcorr(1:nBas,1:nBas)
  R_ao(1:nBas           ,nBas+1:nBas_twice)=Panomcorr(1:nBas,1:nBas)
  R_ao(nBas+1:nBas_twice,1:nBas           )=Panomcorr(1:nBas,1:nBas)
  R_ao(nBas+1:nBas_twice,nBas+1:nBas_twice)=matmul(X(1:nBas,1:nOrb), transpose(X(1:nBas,1:nOrb)))-0.5d0*Pcorr(1:nBas,1:nBas)
!    write(*,*) 'R_linGWB_ao iter '
!    do ibas=1,nBas_twice
!     write(*,'(*(f10.5))') R_ao(ibas,:)
!    enddo

! Compute new total energy, r^can eigenvalues, and Occ numbers

  allocate(J(nBas,nBas))
  allocate(K(nBas,nBas))
  allocate(Delta(nBas,nBas))
  allocate(Occ(nOrb))
  allocate(Occ_R(nOrb_twice))
  allocate(c_inv(nOrb,nBas))
  c_inv(:,:) = matmul(transpose(c),S)
  Pcorr_mo(1:nOrb,1:nOrb)     = 2d0*Rcorr(1:nOrb,1:nOrb)        &
                              + matmul(matmul(c_inv,P),transpose(c_inv)) 
  Panomcorr_mo(1:nOrb,1:nOrb) = Rcorr(1:nOrb,nOrb+1:nOrb_twice) &
                              + matmul(matmul(c_inv,Panom),transpose(c_inv))
  Rcorr(1:nOrb           ,1:nOrb           ) =  Rcorr(1:nOrb           ,1:nOrb           ) &
                                             + 0.5d0*matmul(matmul(c_inv,P),transpose(c_inv))
  Rcorr(1:nOrb           ,nOrb+1:nOrb_twice) =  Rcorr(1:nOrb           ,nOrb+1:nOrb_twice) &
                                             + matmul(matmul(c_inv,Panom),transpose(c_inv))
  Rcorr(nOrb+1:nOrb_twice,1:nOrb           ) =  Rcorr(nOrb+1:nOrb_twice,1:nOrb           ) &
                                             + matmul(matmul(c_inv,Panom),transpose(c_inv))
  Rcorr(nOrb+1:nOrb_twice,nOrb+1:nOrb_twice) =  Rcorr(nOrb+1:nOrb_twice,nOrb+1:nOrb_twice) &
                                             - 0.5d0*matmul(matmul(c_inv,P),transpose(c_inv))
  do iorb=1,nOrb
   Rcorr(iorb+nOrb,iorb+nOrb) = Rcorr(iorb+nOrb,iorb+nOrb) + 1d0
  enddo
  call diagonalize_matrix(nOrb_twice,Rcorr,Occ_R)
  trace_r_can=0d0
  dev_Idemp_r_can=0d0
  do iorb=1,nOrb
   trace_r_can=trace_r_can+Occ_R(iorb)
   dev_Idemp_r_can=dev_Idemp_r_can+abs(Occ_R(iorb))
  enddo
  do iorb=nOrb+1,nOrb_twice
   trace_r_can=trace_r_can+Occ_R(iorb)
   dev_Idemp_r_can=dev_Idemp_r_can+abs(1d0-Occ_R(iorb))
  enddo

  N_anom = trace_matrix(nBas,matmul(transpose(Panomcorr),Panomcorr))
  call diagonalize_matrix(nOrb,Pcorr_mo,Occ)
  call Hartree_matrix_AO_basis(nBas,Pcorr,ERI,J)
  call exchange_matrix_AO_basis(nBas,Pcorr,ERI,K)
  call anomalous_matrix_AO_basis(nBas,sigma,Panomcorr,ERI,Delta)

  ! Kinetic energy

  ET = trace_matrix(nBas,matmul(Pcorr,T))

  ! Potential energy

  EV = trace_matrix(nBas,matmul(Pcorr,V))

  ! Hartree energy

  EJ = 0.5d0*trace_matrix(nBas,matmul(Pcorr,J))

  ! Exchange energy

  EK = 0.25d0*trace_matrix(nBas,matmul(Pcorr,K))

  ! Anomalous energy

  EL = trace_matrix(nBas,matmul(Panomcorr,Delta))

  ! Total energy (incl. the Galitkii-Migdal contribution)

  ElinG = ET + EV + EJ + EK + EL + EcGM

! Print results

  write(*,*)
  write(*,'(A50)')           '---------------------------------------'
  write(*,'(A33)')           ' Summary               '
  write(*,'(A50)')           '---------------------------------------'
  write(*,'(A33,1X,F16.10,A3)') ' One-electron energy = ',ET + EV,' au'
  write(*,'(A33,1X,F16.10,A3)') ' Kinetic      energy = ',ET,' au'
  write(*,'(A33,1X,F16.10,A3)') ' Potential    energy = ',EV,' au'
  write(*,'(A50)')           '---------------------------------------'
  write(*,'(A33,1X,F16.10,A3)') ' Two-electron energy = ',EJ + EK + EL + EcGM,' au'
  write(*,'(A33,1X,F16.10,A3)') ' Hartree      energy = ',EJ,' au'
  write(*,'(A33,1X,F16.10,A3)') ' Exchange     energy = ',EK,' au'
  write(*,'(A33,1X,F16.10,A3)') ' Anomalous    energy = ',EL,' au'
  write(*,'(A33,1X,F16.10,A3)') ' GM           energy = ',EcGM,' au'
  write(*,'(A50)')           '---------------------------------------'
  write(*,'(A33,1X,F16.10,A3)') ' Electronic   energy = ',ElinG,' au'
  write(*,'(A33,1X,F16.10,A3)') ' Nuclear   repulsion = ',ENuc,' au'
  write(*,'(A33,1X,F16.10,A3)') ' B-lin-Dyson  energy = ',ElinG + ENuc,' au'
  write(*,'(A50)')           '---------------------------------------'
  write(*,'(A33,1X,F16.10,A3)') ' | Anomalous dens |  = ',N_anom,'   '
  write(*,'(A33,1X,F16.10,A3)') ' Dev. Idemp. r^can   = ',dev_Idemp_r_can,'   '
  write(*,'(A33,1X,F16.10,A3)') '        Tr[ r^can ]  = ',trace_r_can,'   '
  write(*,'(A50)')           '---------------------------------------'
  write(*,'(A50)') '-----------------------------------------'
  write(*,'(A50)') ' Bogoliubov lin-Dyson occupation numbers '
  write(*,'(A50)') '-----------------------------------------'
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

  deallocate(Sigma_c_he)
  deallocate(Sigma_c_hh)
  deallocate(Sigma_c_eh)
  deallocate(Sigma_c_ee)
  deallocate(Mat1,Mat2)
  deallocate(J,K,Delta,Occ,Occ_R,c_inv)
  deallocate(Rcorr,R_ao,Sigma_c_QP,Tmp_QP,Tmp_mo,Pcorr_mo,Panomcorr_mo,wweight2,wcoord2,wcoord2_cpx)

end subroutine

