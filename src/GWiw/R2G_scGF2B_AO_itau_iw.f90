
subroutine Res_2_Gen_scGF2B(nBas,nOrb,nOrb_twice,nfreqs,no_fock,Enuc,sigma,chem_pot,maxSCF_GF,max_diis_GF,thresh_GF,restart_scGF2,verbose_scGF2,chem_pot_scG, &
                            wcoord,wweight,Fock,Delta,Hc,S,X,pMAT,panomMAT,MOCoef,ERI_AO)

  logical,intent(in)             :: no_fock
  logical,intent(in)             :: restart_scGF2
  logical,intent(in)             :: verbose_scGF2
  logical,intent(in)             :: chem_pot_scG
  integer,intent(in)             :: nfreqs
  integer,intent(in)             :: maxSCF_GF,max_diis_GF
  integer,intent(in)             :: nBas,nOrb,nOrb_twice
  double precision,intent(in)    :: ENuc
  double precision,intent(in)    :: sigma
  double precision,intent(in)    :: chem_pot
  double precision,intent(in)    :: thresh_GF
  double precision,intent(in)    :: wcoord(nfreqs)
  double precision,intent(in)    :: wweight(nfreqs)
  double precision,intent(in)    :: S(nBas,nBas)
  double precision,intent(in)    :: X(nBas,nOrb)
  double precision,intent(in)    :: Hc(nBas,nBas)
  double precision,intent(inout) :: Fock(nBas,nBas)
  double precision,intent(inout) :: Delta(nBas,nBas)
  double precision,intent(in)    :: pMAT(nBas,nBas)
  double precision,intent(in)    :: panomMAT(nBas,nBas)
  double precision,intent(in)    :: MOCoef(nBas,nOrb)
  double precision,intent(in)    :: ERI_AO(nBas,nBas,nBas,nBas)

  integer                        :: ibas,jbas,kbas,lbas
  integer                        :: nBas_twice
  integer                        :: nBas_3times
  integer                        :: nBas_4times
  integer                        :: nOrb_3times
  integer                        :: nOrb_4times
  double precision               :: trace_1rdm
  double precision,allocatable   :: Gen_eQP_state(:)
  double precision,allocatable   :: J(:,:)
  double precision,allocatable   :: K(:,:)
  double precision,allocatable   :: Gen_Hc(:,:)
  double precision,allocatable   :: Gen_S(:,:)
  double precision,allocatable   :: Gen_R(:,:)
  double precision,allocatable   :: Gen_R_mo(:,:)
  double precision,allocatable   :: Gen_MOCoef(:,:)
  double precision,allocatable   :: Gen_H_HFB_ao(:,:)
  double precision,allocatable   :: Gen_MOCoef4(:,:)
  double precision,allocatable   :: Gen_U_QP(:,:)
  double precision,allocatable   :: sw_ERI_AO(:,:,:,:)
  double precision,allocatable   :: db_ERI_AO(:,:,:,:)

  nBas_twice=2*nBas
  nBas_3times=nBas+nBas_twice
  nBas_4times=2*nBas_twice
  nOrb_3times=nOrb+nOrb_twice
  nOrb_4times=2*nOrb_twice
  allocate(J(nBas,nBas))
  allocate(K(nBas,nBas))
  allocate(Gen_eQP_state(nOrb_4times))
  allocate(Gen_Hc(nBas_twice,nBas_twice))
  allocate(Gen_S(nBas_twice,nBas_twice))
  allocate(Gen_R(nBas_4times,nBas_4times))
  allocate(Gen_R_mo(nOrb_4times,nOrb_4times))
  allocate(Gen_H_HFB_ao(nBas_4times,nBas_4times))
  allocate(Gen_MOCoef(nBas_twice,nOrb_twice))
  allocate(Gen_MOCoef4(nBas_4times,nOrb_4times))
  allocate(Gen_U_QP(nOrb_4times,nOrb_4times))
  allocate(sw_ERI_AO(nBas_twice,nBas_twice,nBas_twice,nBas_twice))
  allocate(db_ERI_AO(nBas_twice,nBas_twice,nBas_twice,nBas_twice))

  Gen_Hc=0d0
  Gen_Hc(1:nBas           ,1:nBas           ) = Hc(1:nBas,1:nBas)
  Gen_Hc(nBas+1:nBas_twice,nBas+1:nBas_twice) = Hc(1:nBas,1:nBas)
  Gen_S=0d0
  Gen_S(1:nBas           ,1:nBas           ) = S(1:nBas,1:nBas)
  Gen_S(nBas+1:nBas_twice,nBas+1:nBas_twice) = S(1:nBas,1:nBas)
  Gen_R=0d0
  Gen_R(1:nBas                   ,1:nBas                   ) = 0.5d0*pMAT(1:nBas,1:nBas)
  Gen_R(nBas+1:nBas_twice        ,nBas+1:nBas_twice        ) = 0.5d0*pMAT(1:nBas,1:nBas)
  Gen_R(nBas_twice+1:nBas_3times ,nBas_twice+1:nBas_3times ) = matmul(X(1:nBas,1:nOrb), transpose(X(1:nBas,1:nOrb)))-0.5d0*pMAT(1:nBas,1:nBas)
  Gen_R(nBas_3times+1:nBas_4times,nBas_3times+1:nBas_4times) = Gen_R(nBas_twice+1:nBas_3times,nBas_twice+1:nBas_3times)
  Gen_R(1:nBas                   ,nBas_3times+1:nBas_4times) =  panomMAT(1:nBas,1:nBas)
  Gen_R(nBas+1:nBas_twice        ,nBas_twice+1:nBas_3times ) = -panomMAT(1:nBas,1:nBas)
  Gen_R(nBas_twice+1:nBas_3times ,nBas+1:nBas_twice        ) = -panomMAT(1:nBas,1:nBas)
  Gen_R(nBas_3times+1:nBas_4times,1:nBas                   ) =  panomMAT(1:nBas,1:nBas)
  Gen_MOCoef=0d0
  Gen_MOCoef(1:nBas           ,1:nOrb           ) = MOCoef(1:nBas,1:nOrb)
  Gen_MOCoef(nBas+1:nBas_twice,nOrb+1:nOrb_twice) = MOCoef(1:nBas,1:nOrb)
  Gen_MOCoef4=0d0
  Gen_MOCoef4(1:nBas                   ,1:nOrb                   ) = MOCoef(1:nBas,1:nOrb)
  Gen_MOCoef4(nBas+1:nBas_twice        ,nOrb+1:nOrb_twice        ) = MOCoef(1:nBas,1:nOrb)
  Gen_MOCoef4(nBas_twice+1:nBas_3times ,nOrb_twice+1:nOrb_3times ) = MOCoef(1:nBas,1:nOrb)
  Gen_MOCoef4(nBas_3times+1:nBas_4times,nOrb_3times+1:nOrb_4times) = MOCoef(1:nBas,1:nOrb)
  sw_ERI_AO=0d0
  sw_ERI_AO(1:nBas           ,1:nBas           ,1:nBas           ,1:nBas           ) = ERI_AO(1:nBas,1:nBas,1:nBas,1:nBas)  ! a a a a
  sw_ERI_AO(1:nBas           ,nBas+1:nBas_twice,1:nBas           ,nBas+1:nBas_twice) = ERI_AO(1:nBas,1:nBas,1:nBas,1:nBas)  ! a b a b
  sw_ERI_AO(nBas+1:nBas_twice,1:nBas           ,nBas+1:nBas_twice,1:nBas           ) = ERI_AO(1:nBas,1:nBas,1:nBas,1:nBas)  ! b a b a
  sw_ERI_AO(nBas+1:nBas_twice,nBas+1:nBas_twice,nBas+1:nBas_twice,nbas+1:nBas_twice) = ERI_AO(1:nBas,1:nBas,1:nBas,1:nBas)  ! b b b b
  do ibas=1,nBas_twice
   do jbas=1,nBas_twice
    do kbas=1,nBas_twice
     do lbas=1,nBas_twice
      db_ERI_AO(ibas,jbas,kbas,lbas)=sw_ERI_AO(ibas,jbas,kbas,lbas)-sw_ERI_AO(ibas,jbas,lbas,kbas)
     enddo
    enddo
   enddo
  enddo

  call Hartree_matrix_AO_basis(nBas,pMAT,ERI_AO,J)
  call exchange_matrix_AO_basis(nBas,pMAT,ERI_AO,K)
  call anomalous_matrix_AO_basis(nBas,sigma,panomMAT,ERI_AO,Delta)
  Fock(:,:) = Hc(:,:) + J(:,:) + 0.5d0*K(:,:) - chem_pot*S(:,:)
  Gen_H_HFB_ao(:,:)    = 0d0
  Gen_H_HFB_ao(1:nBas                   ,1:nBas                   ) = Fock(1:nBas,1:nBas)
  Gen_H_HFB_ao(nBas+1:nBas_twice        ,nBas+1:nBas_twice        ) = Fock(1:nBas,1:nBas)
  Gen_H_HFB_ao(nBas_twice+1:nBas_3times ,nBas_twice+1:nBas_3times ) =-Fock(1:nBas,1:nBas)
  Gen_H_HFB_ao(nBas_3times+1:nBas_4times,nBas_3times+1:nBas_4times) =-Fock(1:nBas,1:nBas)
  Gen_H_HFB_ao(1:nBas                   ,nBas_3times+1:nBas_4times) = Delta(1:nBas,1:nBas)
  Gen_H_HFB_ao(nBas+1:nBas_twice        ,nBas_twice+1:nBas_3times ) =-Delta(1:nBas,1:nBas)
  Gen_H_HFB_ao(nBas_twice+1:nBas_3times ,nBas+1:nBas_twice        ) =-Delta(1:nBas,1:nBas)
  Gen_H_HFB_ao(nBas_3times+1:nBas_4times,1:nBas                   ) = Delta(1:nBas,1:nBas)
  Gen_U_QP = matmul(transpose(Gen_MOCoef4),matmul(Gen_H_HFB_ao,Gen_MOCoef4))
  call diagonalize_matrix(nOrb_4times,Gen_U_QP,Gen_eQP_state)
  !Gen_R_mo(:,:)     = 0d0
  !do iorb=1,nOrb_twice
  ! Gen_R_mo(:,:) = Gen_R_mo(:,:) + matmul(Gen_U_QP(:,iorb:iorb),transpose(Gen_U_QP(:,iorb:iorb))) 
  !enddo
  !trace_1rdm = 0d0
  !do iorb=1,nOrb_twice
  ! trace_1rdm = trace_1rdm+Gen_R_mo(iorb,iorb) 
  !enddo
  deallocate(Gen_R_mo)
  deallocate(J)
  deallocate(K)
  deallocate(Gen_H_HFB_ao)
  deallocate(Gen_MOCoef4)
  deallocate(sw_ERI_AO)
  call scGGF2B_AO_itau_iw(nBas_twice,nBas_4times,nOrb_twice,nOrb_4times,maxSCF_GF,thresh_GF,max_diis_GF,restart_scGF2,        &
                          verbose_scGF2,chem_pot_scG,no_fock,ENuc,Gen_Hc,Gen_S,Gen_R,Gen_MOCoef,Gen_eQP_state,chem_pot,sigma, &
                          nfreqs,wcoord,wweight,Gen_U_QP,db_ERI_AO)
  deallocate(Gen_Hc)
  deallocate(Gen_S)
  deallocate(Gen_R)
  deallocate(Gen_MOCoef)
  deallocate(Gen_eQP_state)
  deallocate(Gen_U_QP)
  deallocate(db_ERI_AO)

end subroutine
