subroutine unrestricted_individual_energy(x_rung,x_DFA,c_rung,c_DFA,LDA_centered,nEns,wEns,nGrid,weight,nBas,AO,dAO, &
                                          nO,nV,T,V,ERI,ENuc,eps,Pw,rhow,drhow,J,Fx,FxHF,Fc,P,rho,drho,Ew,E,Om)

! Compute unrestricted individual energies as well as excitation energies

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: x_rung,c_rung
  character(len=12),intent(in)  :: x_DFA,c_DFA
  logical,intent(in)            :: LDA_centered
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  integer,intent(in)            :: nBas
  double precision,intent(in)   :: AO(nBas,nGrid)
  double precision,intent(in)   :: dAO(ncart,nBas,nGrid)

  integer,intent(in)            :: nO(nspin),nV(nspin)
  double precision,intent(in)   :: T(nBas,nBas)
  double precision,intent(in)   :: V(nBas,nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ENuc

  double precision,intent(in)   :: eps(nBas,nspin)     !!!!!
  double precision,intent(in)   :: Pw(nBas,nBas,nspin) !!!!!
  double precision,intent(in)   :: rhow(nGrid,nspin)
  double precision,intent(in)   :: drhow(ncart,nGrid,nspin)

  double precision,intent(in)   :: P(nBas,nBas,nspin,nEns)
  double precision,intent(in)   :: rho(nGrid,nspin,nEns)
  double precision,intent(in)   :: drho(ncart,nGrid,nspin,nEns)

  double precision,intent(in)   :: J(nBas,nBas,nspin)
  double precision,intent(in)   :: Fx(nBas,nBas,nspin)
  double precision,intent(in)   :: FxHF(nBas,nBas,nspin)
  double precision,intent(in)   :: Fc(nBas,nBas,nspin)
  double precision              :: Ew


! Local variables

  double precision              :: ET(nspin,nEns)
  double precision              :: EV(nspin,nEns)
  double precision              :: EJ(nsp,nEns)
  double precision              :: Ex(nspin,nEns)
  double precision              :: Ec(nsp,nEns)
  double precision              :: Exc(nEns)  
  double precision              :: Eaux(nspin,nEns) 

  double precision              :: ExDD(nspin,nEns) 
  double precision              :: EcDD(nsp,nEns) 
  double precision              :: ExcDD(nsp,nEns)

  double precision              :: Omx(nEns),Omc(nEns),Omxc(nEns) 
  double precision              :: Omaux(nEns)                   
  double precision              :: OmxDD(nEns),OmcDD(nEns),OmxcDD(nEns) 

  double precision,external     :: trace_matrix

  integer                       :: ispin,iEns

! Output variables

  double precision,intent(out)  :: E(nEns)
  double precision,intent(out)  :: Om(nEns)

!------------------------------------------------------------------------
! Kinetic energy
!------------------------------------------------------------------------

  do ispin=1,nspin
    do iEns=1,nEns
      ET(ispin,iEns) = trace_matrix(nBas,matmul(P(:,:,ispin,iEns),T(:,:)))
    end do
  end do

!------------------------------------------------------------------------
! Potential energy
!------------------------------------------------------------------------

  do iEns=1,nEns
    do ispin=1,nspin
      EV(ispin,iEns) = trace_matrix(nBas,matmul(P(:,:,ispin,iEns),V(:,:)))
    end do
  end do

!------------------------------------------------------------------------
! Individual Hartree energy
!------------------------------------------------------------------------

  do iEns=1,nEns

    do ispin=1,nspin
      call hartree_coulomb(nBas,Pw(:,:,ispin),ERI,J(:,:,ispin))
    end do

    EJ(1,iEns) =       trace_matrix(nBas,matmul(P(:,:,1,iEns),J(:,:,1))) &
               - 0.5d0*trace_matrix(nBas,matmul(Pw(:,:,1),J(:,:,1)))

    EJ(2,iEns) = 2.0d0*trace_matrix(nBas,matmul(P(:,:,1,iEns),J(:,:,2))) &
               - 0.5d0*trace_matrix(nBas,matmul(Pw(:,:,1),J(:,:,2)))     &
               - 0.5d0*trace_matrix(nBas,matmul(Pw(:,:,2),J(:,:,1)))

    EJ(3,iEns) =       trace_matrix(nBas,matmul(P(:,:,2,iEns),J(:,:,2))) &
               - 0.5d0*trace_matrix(nBas,matmul(Pw(:,:,2),J(:,:,2)))

  end do

!------------------------------------------------------------------------
! Individual exchange energy
!------------------------------------------------------------------------

  do iEns=1,nEns
    do ispin=1,nspin
      call exchange_individual_energy(x_rung,x_DFA,LDA_centered,nEns,wEns,nGrid,weight,nBas,ERI,      &
                                      Pw(:,:,ispin),P(:,:,ispin,iEns),rhow(:,ispin),drhow(:,:,ispin), & 
                                      rho(:,ispin,iEns),drho(:,:,ispin,iEns),Ex(ispin,iEns))
    end do
  end do

!------------------------------------------------------------------------
! Individual correlation energy
!------------------------------------------------------------------------

  do iEns=1,nEns
    call unrestricted_correlation_individual_energy(c_rung,c_DFA,LDA_centered,nEns,wEns,nGrid,weight, &
                                                    rhow,drhow,rho(:,:,iEns),drho(:,:,:,iEns),Ec(:,iEns))
  end do

!------------------------------------------------------------------------
! Compute auxiliary energies
!------------------------------------------------------------------------
 
  call unrestricted_auxiliary_energy(nBas,nEns,nO(:),eps(:,:),Eaux(:,:))

!------------------------------------------------------------------------
! Compute derivative discontinuities
!------------------------------------------------------------------------

  do ispin=1,nspin 

    call exchange_derivative_discontinuity(x_rung,x_DFA,nEns,wEns,nGrid,weight, &
                                           rhow(:,ispin),drhow(:,:,ispin),ExDD(ispin,:))
  end do

  call unrestricted_correlation_derivative_discontinuity(c_rung,c_DFA,nEns,wEns,nGrid,weight,rhow,drhow,EcDD)

  ExcDD(1,:) = ExDD(1,:) + EcDD(1,:)
  ExcDD(2,:) =             EcDD(2,:)
  ExcDD(3,:) = ExDD(2,:) + EcDD(3,:)

!------------------------------------------------------------------------
! Total energy
!------------------------------------------------------------------------

  do iEns=1,nEns
    Exc(iEns) = sum(Ex(:,iEns)) + sum(Ec(:,iEns))
    E(iEns)   = sum(ET(:,iEns)) + sum(EV(:,iEns)) + sum(EJ(:,iEns)) &
              + sum(Ex(:,iEns)) + sum(Ec(:,iEns)) + sum(ExcDD(:,iEns))
  end do

!------------------------------------------------------------------------
! Excitation energies
!------------------------------------------------------------------------

  do iEns=1,nEns
    Om(iEns) = E(iEns) - E(1)

    Omx(iEns)    = sum(Ex(:,iEns)) - sum(Ex(:,1))
    Omc(iEns)    = sum(Ec(:,iEns)) - sum(Ec(:,1))
    Omxc(iEns)   =     Exc(iEns)   -     Exc(1)

    Omaux(iEns)  = sum(Eaux(:,iEns))  - sum(Eaux(:,1))

    OmxDD(iEns)  = sum(ExDD(:,iEns))  - sum(ExDD(:,1))
    OmcDD(iEns)  = sum(EcDD(:,iEns))  - sum(EcDD(:,1))
    OmxcDD(iEns) = sum(ExcDD(:,iEns)) - sum(ExcDD(:,1))


  end do

!------------------------------------------------------------------------
! Dump results
!------------------------------------------------------------------------

  call print_unrestricted_individual_energy(nEns,ENuc,Ew,ET,EV,EJ,Ex,Ec,Exc,Eaux,ExDD,EcDD,ExcDD,E, & 
                                          Om,Omx,Omc,Omxc,Omaux,OmxDD,OmcDD,OmxcDD)

end subroutine unrestricted_individual_energy
