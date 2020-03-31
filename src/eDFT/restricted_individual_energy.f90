subroutine restricted_individual_energy(x_rung,x_DFA,c_rung,c_DFA,nEns,wEns,nGrid,weight,nBas,        & 
                                        nO,nV,T,V,ERI,ENuc,eps,Pw,rhow,drhow,J,P,rho,drho,Ew,EwGIC,E,Om)

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

  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  double precision,intent(in)   :: T(nBas,nBas)
  double precision,intent(in)   :: V(nBas,nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ENuc

  double precision,intent(in)   :: eps(nBas)
  double precision,intent(in)   :: Pw(nBas,nBas)
  double precision,intent(in)   :: rhow(nGrid)
  double precision,intent(in)   :: drhow(ncart,nGrid)

  double precision,intent(in)   :: P(nBas,nBas,nEns)
  double precision,intent(in)   :: rho(nGrid,nEns)
  double precision,intent(in)   :: drho(ncart,nGrid,nEns)

  double precision,intent(in)   :: J(nBas,nBas)

  double precision              :: Ew

! Local variables

  double precision              :: ET(nEns)
  double precision              :: EV(nEns)
  double precision              :: EJ(nEns)
  double precision              :: Ex(nEns),   Ec(nEns),   Exc(nEns)
  double precision              :: Eaux(nEns)
  double precision              :: ExDD(nEns), EcDD(nEns), ExcDD(nEns)
  double precision              :: Omx(nEns),  Omc(nEns),  Omxc(nEns)
  double precision              :: Omaux(nEns)
  double precision              :: OmxDD(nEns),OmcDD(nEns),OmxcDD(nEns)

  double precision,external     :: trace_matrix

  integer                       :: iEns

! Output variables

  double precision,intent(out)  :: EwGIC
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
    call hartree_coulomb(nBas,P(:,:,iEns),ERI(:,:,:,:),J(:,:))
    EJ(iEns) = 0.5d0*trace_matrix(nBas,matmul(P(:,:,iEns),J(:,:)))
  end do

!------------------------------------------------------------------------
! Individual exchange energy
!------------------------------------------------------------------------

  do iEns=1,nEns

      call exchange_individual_energy(x_rung,x_DFA,nEns,wEns(:),nGrid,weight(:),nBas,ERI(:,:,:,:), &
                                      P(:,:,iEns),rhow(:),drhow(:,:),rho(:,iEns),drho(:,:,iEns),Ex(iEns))

  end do

!------------------------------------------------------------------------
! Indivudual correlation energy
!------------------------------------------------------------------------

  do iEns=1,nEns

    call restricted_correlation_individual_energy(c_rung,c_DFA,nEns,wEns(:),nGrid,weight(:),rhow(:),drhow(:,:), & 
                                                  rho(:,iEns),drho(:,:,iEns),Ec(iEns))

  end do




!------------------------------------------------------------------------
! Compute auxiliary energies
!------------------------------------------------------------------------

  call restricted_auxiliary_energy(nBas,nEns,nO,eps(:),Eaux(:))

!------------------------------------------------------------------------
! Compute derivative discontinuities
!------------------------------------------------------------------------

  call exchange_derivative_discontinuity(x_rung,x_DFA,nEns,wEns(:),nGrid,weight(:),rhow(:),drhow(:,:),ExDD(:))

  call restricted_correlation_derivative_discontinuity(c_rung,c_DFA,nEns,wEns(:),nGrid,weight(:),rhow(:),drhow(:,:),EcDD(:))

  ExcDD(:) = ExDD(:) + EcDD(:)

!------------------------------------------------------------------------
! Total energy
!------------------------------------------------------------------------

  do iEns=1,nEns
    Exc(iEns) = Ex(iEns) + Ec(iEns)
    E(iEns)   = ET(iEns) + EV(iEns) + EJ(iEns) &
              + Ex(iEns) + Ec(iEns) + ExcDD(iEns)
  end do

!------------------------------------------------------------------------
! Total energy with ghost-interaction correction
!------------------------------------------------------------------------

  EwGIC = 0d0
  do iEns=1,nEns
    EwGIC = EwGIC + wEns(iEns)*E(iEns)
  end do

!------------------------------------------------------------------------
! Excitation energies
!------------------------------------------------------------------------

  do iEns=1,nEns

    Om(iEns)     = E(iEns)     - E(1)

    Omx(iEns)    = Ex(iEns)    - Ex(1)
    Omc(iEns)    = Ec(iEns)    - Ec(1)
    Omxc(iEns)   = Exc(iEns)   - Exc(1)

    Omaux(iEns)  = Eaux(iEns)  - Eaux(1)

    OmxDD(iEns)  = ExDD(iEns)  - ExDD(1)
    OmcDD(iEns)  = EcDD(iEns)  - EcDD(1)
    OmxcDD(iEns) = ExcDD(iEns) - ExcDD(1)

  end do

!------------------------------------------------------------------------
! Dump results
!------------------------------------------------------------------------

  call print_restricted_individual_energy(nEns,ENuc,Ew,EwGIC,ET(:),EV(:),EJ(:),Ex(:),Ec(:),Exc(:), &
                                          Eaux(:),ExDD(:),EcDD(:),ExcDD(:),E(:),                   & 
                                          Om(:),Omx(:),Omc(:),Omxc(:),Omaux,OmxDD(:),OmcDD(:),OmxcDD(:))

end subroutine restricted_individual_energy
