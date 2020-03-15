subroutine restricted_individual_energy(x_rung,x_DFA,c_rung,c_DFA,nEns,wEns,nGrid,weight,nBas,AO,dAO, & 
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

  integer,intent(in)            :: nO,nV
  double precision,intent(in)   :: T(nBas,nBas)
  double precision,intent(in)   :: V(nBas,nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ENuc

  double precision,intent(in)   :: Pw(nBas,nBas)
  double precision,intent(in)   :: rhow(nGrid)
  double precision,intent(in)   :: drhow(ncart,nGrid)

  double precision,intent(in)   :: P(nBas,nBas,nEns)
  double precision,intent(in)   :: rho(nGrid,nEns)
  double precision,intent(in)   :: drho(ncart,nGrid,nEns)

  double precision,intent(in)   :: J(nBas,nBas)
  double precision,intent(in)   :: Fx(nBas,nBas)
  double precision,intent(in)   :: FxHF(nBas,nBas)
  double precision,intent(in)   :: Fc(nBas,nBas)

! Local variables

  double precision              :: ET(nEns)
  double precision              :: EV(nEns)
  double precision              :: EJ(nEns)
  double precision              :: Ex(nEns),Ec(nEns),Exc(nEns)
  double precision              :: ExLZ,EcLZ,ExcLZ
  double precision              :: ExDD(nEns),EcDD(nEns),ExcDD(nEns)

  double precision,external     :: trace_matrix

  integer                       :: iEns

! Output variables

  double precision,intent(out)  :: E(nEns)
  double precision,intent(out)  :: Om(nEns)

!------------------------------------------------------------------------
! Kinetic energy
!------------------------------------------------------------------------

  do iEns=1,nEns
    ET(iEns) = trace_matrix(nBas,matmul(P(:,:,iEns),T(:,:)))
  end do

!------------------------------------------------------------------------
! Potential energy
!------------------------------------------------------------------------

  do iEns=1,nEns
    EV(iEns) = trace_matrix(nBas,matmul(P(:,:,iEns),V(:,:)))
  end do

!------------------------------------------------------------------------
! Hartree energy
!------------------------------------------------------------------------

  do iEns=1,nEns
    call hartree_coulomb(nBas,P(:,:,iEns),ERI,J(:,:))
    EJ(iEns) = 0.5d0*trace_matrix(nBas,matmul(P(:,:,iEns),J(:,:)))
  end do

!------------------------------------------------------------------------
! Exchange energy
!------------------------------------------------------------------------

  do iEns=1,nEns

      call exchange_potential(x_rung,x_DFA,nEns,wEns(:),nGrid,weight(:),nBas,P(:,:,iEns),ERI(:,:,:,:), & 
                              AO(:,:),dAO(:,:,:),rho(:,iEns),drho(:,:,iEns),Fx(:,:),FxHF(:,:))
      call exchange_energy(x_rung,x_DFA,nEns,wEns(:),nGrid,weight(:),nBas,P(:,:,iEns),FxHF(:,:), &
                           rho(:,iEns),drho(:,:,iEns),Ex(iEns))

  end do

  Ex(:) = 0.5d0*Ex(:)

!------------------------------------------------------------------------
! Correlation energy
!------------------------------------------------------------------------

  do iEns=1,nEns

    call restricted_correlation_individual_energy(c_rung,c_DFA,nEns,wEns(:),nGrid,weight(:),rhow(:),drhow(:,:), & 
                                       rho(:,iEns),drho(:,:,iEns),Ec(iEns))

  end do

!------------------------------------------------------------------------
! Compute Levy-Zahariev shift
!------------------------------------------------------------------------

  call restricted_correlation_Levy_Zahariev_shift(c_rung,c_DFA,nEns,wEns(:),nGrid,weight(:),rho(:,:),drho(:,:,:), & 
                                       ExLZ,EcLZ,ExcLZ)

!------------------------------------------------------------------------
! Compute derivative discontinuities
!------------------------------------------------------------------------

  call restricted_correlation_derivative_discontinuity(c_rung,c_DFA,nEns,wEns(:),nGrid,weight(:),rhow(:),drhow(:,:), & 
                                            ExDD(:),EcDD(:),ExcDD(:))

!------------------------------------------------------------------------
! Total energy
!------------------------------------------------------------------------

  do iEns=1,nEns
    Exc(iEns) = Ex(iEns) + Ec(iEns)
    E(iEns) = ENuc + ET(iEns) + EV(iEns) + EJ(iEns) &
                   + Ex(iEns) + Ec(iEns) + ExcLZ + ExcDD(iEns)
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

  call print_individual_energy(nEns,ET(:),EV(:),EJ(:),Ex(:),Ec(:),Exc(:),ExLZ,EcLZ,ExcLZ, & 
                               ExDD(:),EcDD(:),ExcDD(:),E(:),Om(:))

end subroutine restricted_individual_energy
