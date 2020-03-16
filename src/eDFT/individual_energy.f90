subroutine individual_energy(x_rung,x_DFA,c_rung,c_DFA,nEns,wEns,nGrid,weight,nBas,AO,dAO, & 
                             nO,nV,T,V,ERI,ENuc,Pw,rhow,drhow,J,Fx,FxHF,Fc,P,rho,drho,E,Om)

! Compute individual energies as well as excitation energies

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: x_rung,c_rung
  character(len=12),intent(in)  :: x_DFA,c_DFA
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

  double precision,intent(in)   :: Pw(nBas,nBas)
  double precision,intent(in)   :: rhow(nGrid,nspin)
  double precision,intent(in)   :: drhow(ncart,nGrid,nspin)

  double precision,intent(in)   :: P(nBas,nBas,nspin,nEns)
  double precision,intent(in)   :: rho(nGrid,nspin,nEns)
  double precision,intent(in)   :: drho(ncart,nGrid,nspin,nEns)

  double precision,intent(in)   :: J(nBas,nBas,nspin)
  double precision,intent(in)   :: Fx(nBas,nBas,nspin)
  double precision,intent(in)   :: FxHF(nBas,nBas,nspin)
  double precision,intent(in)   :: Fc(nBas,nBas,nspin)

! Local variables

  double precision              :: ET(nspin,nEns)
  double precision              :: EV(nspin,nEns)
  double precision              :: EJ(nsp,nEns)
  double precision              :: Ex(nspin,nEns)
  double precision              :: Ec(nsp,nEns)
  double precision              :: EcDD(nsp,nEns) 

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
! Hartree energy
!------------------------------------------------------------------------

  do iEns=1,nEns

    do ispin=1,nspin
      call hartree_coulomb(nBas,P(:,:,ispin,iEns),ERI,J(:,:,ispin))
    end do

    EJ(1,iEns) = 0.5d0*trace_matrix(nBas,matmul(P(:,:,1,iEns),J(:,:,1)))
    EJ(2,iEns) = trace_matrix(nBas,matmul(P(:,:,1,iEns),J(:,:,2)))
    EJ(3,iEns) = 0.5d0*trace_matrix(nBas,matmul(P(:,:,2,iEns),J(:,:,2)))

  end do

!------------------------------------------------------------------------
! Exchange energy
!------------------------------------------------------------------------

  do iEns=1,nEns
    do ispin=1,nspin

      call exchange_potential(x_rung,x_DFA,nEns,wEns(:),nGrid,weight(:),nBas,P(:,:,ispin,iEns),ERI(:,:,:,:), & 
                              AO(:,:),dAO(:,:,:),rho(:,ispin,iEns),drho(:,:,ispin,iEns),Fx(:,:,ispin),FxHF(:,:,ispin))
      call exchange_energy(x_rung,x_DFA,nEns,wEns(:),nGrid,weight(:),nBas,P(:,:,ispin,iEns),FxHF(:,:,ispin), &
                           rho(:,ispin,iEns),drho(:,:,ispin,iEns),Ex(ispin,iEns))

    end do
  end do

!------------------------------------------------------------------------
! Correlation energy
!------------------------------------------------------------------------

  do iEns=1,nEns

    call correlation_individual_energy(c_rung,c_DFA,nEns,wEns(:),nGrid,weight(:),rhow(:,:),drhow(:,:,:), & 
                                       rho(:,:,iEns),drho(:,:,:,iEns),Ec(:,iEns))

  end do

!------------------------------------------------------------------------
! Compute derivative discontinuities
!------------------------------------------------------------------------

  call correlation_derivative_discontinuity(c_rung,c_DFA,nEns,wEns(:),nGrid,weight(:),rhow(:,:),drhow(:,:,:),EcDD(:,:))

!------------------------------------------------------------------------
! Total energy
!------------------------------------------------------------------------

  do iEns=1,nEns
    E(iEns) = ENuc + sum(ET(:,iEns)) + sum(EV(:,iEns)) + sum(EJ(:,iEns)) &
                   + sum(Ex(:,iEns)) + sum(Ec(:,iEns)) + sum(EcDD(:,iEns))
  end do

!------------------------------------------------------------------------
! Excitation energies
!------------------------------------------------------------------------

  do iEns=1,nEns
    Om(iEns) = E(iEns) - E(1)
  end do

!------------------------------------------------------------------------
! Dump results
!------------------------------------------------------------------------

  call print_individual_energy(nEns,EJ,Ex,Ec,EcDD,E,Om)

end subroutine individual_energy
