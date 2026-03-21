subroutine BRPAx(nBas,nOrb,TDA,cHFB,Hc,S,ERI,chem_pot,sigma,U_QP,ERHFB,EcRPAx)

! Perform second-order Bogoliubov RPAx calculation

  implicit none

! Input variables

  logical,intent(in)            :: TDA
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  double precision,intent(in)   :: sigma,chem_pot
  double precision,intent(in)   :: ERHFB
  double precision,intent(in)   :: Hc(nBas,nBas) 
  double precision,intent(in)   :: S(nBas,nBas) 
  double precision,intent(in)   :: cHFB(nBas,nOrb)
  double precision,intent(in)   :: U_QP(2*nOrb,2*nOrb)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: iorb
  integer                       :: k,kprime,l,lprime
  integer                       :: iqp_pair,jqp_pair
  integer                       :: nOrb2,nOrb3,nOrb4,nRPA,nRPA2

  double precision              :: Ekkprime
  double precision              :: trace_1rdm
  double precision              :: trace_matrix

  double precision,allocatable  :: Om(:)
  double precision,allocatable  :: eQP_sw(:)
  double precision,allocatable  :: F(:,:),J(:,:),EXX(:,:),Delta(:,:)
  double precision,allocatable  :: P(:,:),Panom(:,:),R(:,:)
  double precision,allocatable  :: Amat(:,:)
  double precision,allocatable  :: Bmat(:,:)
  double precision,allocatable  :: H_RPAx(:,:)
  double precision,allocatable  :: U_QP_sw(:,:)
  double precision,allocatable  :: ERI_MO(:,:,:,:)
  double precision,allocatable  :: ERI_MO_sw(:,:,:,:)
  double precision,allocatable  :: H40(:,:,:,:)
  double precision,allocatable  :: H22(:,:,:,:)
  double precision,allocatable  :: Omega04(:,:,:,:)
  double precision,allocatable  :: Ua(:,:),Va(:,:),U(:,:),V(:,:)

! Output variables

  double precision,intent(out)  :: EcRPAx

! Hello world

  write(*,*)
  write(*,*)'*******************************'
  write(*,*)'* Bogoliubov RPAx Calculation *'
  write(*,*)'*******************************'
  write(*,*)

  nOrb2=2*nOrb
  nOrb3=nOrb2+nOrb
  nOrb4=2*nOrb2
  nRPA=0
  do k=1,nOrb2
   do kprime=k+1,nOrb2
    nRPA=nRPA+1
   enddo
  enddo
  nRPA2=2*nRPA
  allocate(U_QP_sw(nOrb4,nOrb4))
  allocate(eQP_sw(nOrb4))
  allocate(ERI_MO_sw(nOrb2,nOrb2,nOrb2,nOrb2))

! Build P and Panom -> U_QP_sw, eQP_sw
  ! Allocate arrays
  allocate(P(nBas,nBas))
  allocate(Panom(nBas,nBas))
  allocate(R(nOrb2,nOrb2))
  allocate(F(nBas,nBas))
  allocate(J(nBas,nBas))
  allocate(EXX(nBas,nBas))
  allocate(Delta(nBas,nBas))
  trace_1rdm=0d0
  R=0d0
  do iorb=1,nOrb
   R(:,:) = R(:,:) + matmul(U_QP(:,iorb:iorb),transpose(U_QP(:,iorb:iorb))) 
  enddo
  do iorb=1,nOrb
   trace_1rdm = trace_1rdm + R(iorb,iorb) 
  enddo
  P=2d0*matmul(matmul(cHFB,R(1:nOrb,1:nOrb)),transpose(cHFB))
  Panom=matmul(matmul(cHFB,R(1:nOrb,nOrb+1:nOrb2)),transpose(cHFB))
  call Hartree_matrix_AO_basis(nBas,P,ERI,J)
  call exchange_matrix_AO_basis(nBas,P,ERI,EXX)
  call anomalous_matrix_AO_basis(nBas,sigma,Panom,ERI,Delta)
  F(:,:) = Hc(:,:) + J(:,:) + 0.5d0*EXX(:,:) - chem_pot*S(:,:)
  U_QP_sw=0d0
  ! 1 row-block
  U_QP_sw(1:nOrb      ,1:nOrb      )  = matmul(transpose(cHFB),matmul(F,cHFB))
  U_QP_sw(1:nOrb      ,nOrb3+1:nOrb4) = matmul(transpose(cHFB),matmul(Delta,cHFB))
  ! 2 row-block
  U_QP_sw(nOrb+1:nOrb2,nOrb+1:nOrb2)  =  U_QP_sw(1:nOrb      ,1:nOrb)
  U_QP_sw(nOrb+1:nOrb2,nOrb2+1:nOrb3) = -U_QP_sw(1:nOrb      ,nOrb3+1:nOrb4)
  ! 3 row-block
  U_QP_sw(nOrb2+1:nOrb3,nOrb2+1:nOrb3) = -U_QP_sw(1:nOrb      ,1:nOrb)
  U_QP_sw(nOrb2+1:nOrb3,nOrb+1:nOrb2)  =  U_QP_sw(nOrb+1:nOrb2,nOrb2+1:nOrb3)
  ! 4 row-block 
  U_QP_sw(nOrb3+1:nOrb4,nOrb3+1:nOrb4) = -U_QP_sw(1:nOrb      ,1:nOrb)
  U_QP_sw(nOrb3+1:nOrb4,1:nOrb       ) =  U_QP_sw(1:nOrb      ,nOrb3+1:nOrb4)
  call diagonalize_matrix(nOrb4,U_QP_sw,eQP_sw)
  deallocate(P,Panom,R)
  deallocate(F,J,EXX,Delta)

! Allocate NO spin-with ERIs
  allocate(ERI_MO(nOrb,nOrb,nOrb,nOrb))
  call AOtoMO_ERI_RHF(nBas,nOrb,cHFB,ERI,ERI_MO)
  ERI_MO_sw=0d0
  ERI_MO_sw(      1:nOrb,      1:nOrb,      1:nOrb,      1:nOrb)=ERI_MO(1:nOrb,1:nOrb,1:nOrb,1:nOrb) ! a a a a 
  ERI_MO_sw(nOrb+1:nOrb2,nOrb+1:nOrb2,nOrb+1:nOrb2,nOrb+1:nOrb2)=ERI_MO(1:nOrb,1:nOrb,1:nOrb,1:nOrb) ! b b b b 
  ERI_MO_sw(      1:nOrb,nOrb+1:nOrb2,      1:nOrb,nOrb+1:nOrb2)=ERI_MO(1:nOrb,1:nOrb,1:nOrb,1:nOrb) ! a b a b 
  ERI_MO_sw(nOrb+1:nOrb2,      1:nOrb,nOrb+1:nOrb2,      1:nOrb)=ERI_MO(1:nOrb,1:nOrb,1:nOrb,1:nOrb) ! b a b a 
  deallocate(ERI_MO)

! Build Ua, Va, V, U [ Note: we use Va and Ua to define U and V so that k is related to +- Ek ]
  allocate(Ua(nOrb2,nOrb2))
  allocate(Va(nOrb2,nOrb2))
  allocate(V(nOrb2,nOrb2))
  allocate(U(nOrb2,nOrb2))
  Va(1:nOrb2,1:nOrb2) = U_QP_sw(1:nOrb2      ,1:nOrb2)
  Ua(1:nOrb2,1:nOrb2) = U_QP_sw(nOrb2+1:nOrb4,1:nOrb2)
  U(1:nOrb2,1:nOrb2)  = -Va(1:nOrb2,1:nOrb2)
  V(1:nOrb2,1:nOrb2)  =  Ua(1:nOrb2,1:nOrb2)
  deallocate(U_QP_sw)

! Build H40 and H22 (making H40 also be anti-symmetric; using Omega40)
  allocate(H40(nOrb2,nOrb2,nOrb2,nOrb2))
  allocate(H22(nOrb2,nOrb2,nOrb2,nOrb2))
  H40=0d0;  H22=0d0;
  call ERI_MO2QP_H40(nOrb2,ERI_MO_sw,Ua,Va,H40)
  call ERI_MO2QP_H40_2(nOrb2,ERI_MO_sw,Ua,Va,H40)
  call ERI_MO2QP_H40_3(nOrb2,ERI_MO_sw,Ua,Va,H40)
  call ERI_MO2QP_H40_4(nOrb2,ERI_MO_sw,Ua,Va,H40)
  call ERI_MO2QP_H40_5(nOrb2,ERI_MO_sw,Ua,Va,H40)
  call ERI_MO2QP_H40_6(nOrb2,ERI_MO_sw,Ua,Va,H40)
!  H40=0.25d0*H40  ! Remove this factor when using Omega40
  call ERI_MO2QP_H22(nOrb2,ERI_MO_sw,U,Ua,H22)
  call ERI_MO2QP_H22_2(nOrb2,ERI_MO_sw,V,Va,H22)
  call ERI_MO2QP_H22_3(nOrb2,ERI_MO_sw,U,Ua,V,Va,H22)
  call ERI_MO2QP_H22_4(nOrb2,ERI_MO_sw,U,Ua,V,Va,H22)
  call ERI_MO2QP_H22_5(nOrb2,ERI_MO_sw,U,Ua,V,Va,H22)
  call ERI_MO2QP_H22_6(nOrb2,ERI_MO_sw,U,Ua,V,Va,H22)
  deallocate(ERI_MO_sw)
  deallocate(Ua,Va,U,V)

! Prepare A and B
  allocate(Amat(nRPA,nRPA))
  allocate(Bmat(nRPA,nRPA))
  Amat=0d0; Bmat=0d0;
  iqp_pair=1; jqp_pair=1;
  do k=1,nOrb2
   do kprime=k+1,nOrb2
    Ekkprime=-(eQP_sw(k)+eQP_sw(kprime))        ! Using negative energies as positive with a minus
    do l=1,nOrb2
     do lprime=l+1,nOrb2
      Amat(iqp_pair,jqp_pair)=H22(k,kprime,l,lprime)
      Bmat(iqp_pair,jqp_pair)=H40(k,kprime,l,lprime)        ! Version using Omega40
!      Bmat(iqp_pair,jqp_pair)=2.4d1*H40(k,kprime,l,lprime) !         using H40 in Ring and Schuck to be corrected with permutations. 
      if(iqp_pair==jqp_pair) then
       Amat(iqp_pair,jqp_pair)=Amat(iqp_pair,jqp_pair)+Ekkprime
      endif
      jqp_pair=jqp_pair+1
      if(jqp_pair>nRPA) then
       jqp_pair=1
       iqp_pair=iqp_pair+1
      endif
     enddo 
    enddo 
   enddo 
  enddo
  if(TDA) then
   Bmat=0d0
   write(*,*)
   write(*,*) ' Tamm-Dancoff approximation activated!'
   write(*,*)
  endif

! Prepare H_RPAx
  allocate(H_RPAx(nRPA2,nRPA2),Om(nRPA2))
  H_RPAx(1:nRPA      ,1:nRPA      )= Amat(1:nRPA,1:nRPA)
  H_RPAx(1:nRPA      ,nRPA+1:nRPA2)= Bmat(1:nRPA,1:nRPA)
  H_RPAx(nRPA+1:nRPA2,1:nRPA      )=-Bmat(1:nRPA,1:nRPA)
  H_RPAx(nRPA+1:nRPA2,nRPA+1:nRPA2)=-Amat(1:nRPA,1:nRPA)

! Diagonalize QP-RPAx
  call diagonalize_general_matrix(nRPA2,H_RPAx,Om,H_RPAx)
  call sort_ascending(nRPA2,Om)                        

! Compute EcRPAx = 1/2 sum_p  Om_p - A_pp  (for Om_p > 0)
  EcRPAx=0.5d0*( sum(Om(nRPA+1:nRPA2))-trace_matrix(nRPA,Amat) )

! Print TD-HFB excitation energies
  call print_excitation_energies('RPAx@RHFB','Singlet',nRPA,Om(nRPA+1:nRPA2))

  write(*,*)'--------------------------------------------------------------------------'
  write(*,'(A53,F15.6,A3)') 'Tr@RPAx@RHFB correlation energy = ',EcRPAx,' au'
  write(*,'(A53,F15.6,A3)') 'Tr@RPAx@RHFB total energy       = ',ERHFB+EcRPAx,' au'
  write(*,*)'--------------------------------------------------------------------------'
  write(*,*)

  deallocate(H_RPAx)
  deallocate(Amat,Bmat)
  deallocate(H40,H22)
  deallocate(eQP_sw)

end subroutine 
