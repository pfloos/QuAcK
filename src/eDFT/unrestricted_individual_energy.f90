subroutine unrestricted_individual_energy(x_rung,x_DFA,c_rung,c_DFA,LDA_centered,nEns,wEns,nCC,aCC,nGrid,weight,nBas,AO,dAO, &
                                          T,V,ERI,ENuc,eKS,Pw,rhow,drhow,J,Fx,FxHF,Fc,P,rho,drho,occnum,Cx_choice,doNcentered,Ew)

! Compute unrestricted individual energies as well as excitation energies

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: x_rung,c_rung
  integer,intent(in)            :: x_DFA,c_DFA
  logical,intent(in)            :: LDA_centered
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nCC
  double precision,intent(in)   :: aCC(nCC,nEns-1)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  integer,intent(in)            :: nBas
  double precision,intent(in)   :: AO(nBas,nGrid)
  double precision,intent(in)   :: dAO(ncart,nBas,nGrid)

  double precision,intent(in)   :: T(nBas,nBas)
  double precision,intent(in)   :: V(nBas,nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ENuc

  double precision,intent(in)   :: eKS(nBas,nspin)     
  double precision,intent(in)   :: Pw(nBas,nBas,nspin) 
  double precision,intent(in)   :: rhow(nGrid,nspin)
  double precision,intent(in)   :: drhow(ncart,nGrid,nspin)

  double precision,intent(in)   :: P(nBas,nBas,nspin,nEns)
  double precision,intent(in)   :: rho(nGrid,nspin,nEns)
  double precision,intent(in)   :: drho(ncart,nGrid,nspin,nEns)

  double precision,intent(inout):: J(nBas,nBas,nspin)
  double precision,intent(inout):: Fx(nBas,nBas,nspin)
  double precision,intent(inout):: FxHF(nBas,nBas,nspin)
  double precision,intent(inout):: Fc(nBas,nBas,nspin)
  double precision,intent(in)   :: Ew
  double precision,intent(in)   :: occnum(nBas,nspin,nEns)
  integer,intent(in)            :: Cx_choice
  logical,intent(in)            :: doNcentered


! Local variables

  double precision              :: ET(nspin,nEns)
  double precision              :: EV(nspin,nEns)
  double precision              :: EH(nsp,nEns)
  double precision              :: Ex(nspin,nEns)
  double precision              :: Ec(nsp,nEns)
  double precision              :: LZH(nsp)
  double precision              :: LZx(nspin)
  double precision              :: LZc(nsp)
  double precision              :: Eaux(nspin,nEns) 

  double precision              :: ExDD(nspin,nEns) 
  double precision              :: EcDD(nsp,nEns) 

  double precision              :: OmH(nEns)
  double precision              :: Omx(nEns)
  double precision              :: Omc(nEns)
  double precision              :: Omaux(nEns)                   
  double precision              :: OmxDD(nEns)
  double precision              :: OmcDD(nEns)

  double precision,external     :: trace_matrix

  integer                       :: ispin,iEns,iBas
  double precision,allocatable  :: nEl(:)
  double precision,allocatable  :: kappa(:)

  double precision              :: E(nEns)
  double precision              :: Om(nEns)

  double precision,external     :: electron_number

! Compute scaling factor for N-centered ensembles

  allocate(nEl(nEns),kappa(nEns))

  nEl(:) = 0d0
  do iEns=1,nEns
    do iBas=1,nBas
      do ispin=1,nspin
        nEl(iEns) = nEl(iEns) + occnum(iBas,ispin,iEns) 
      end do
    end do
    kappa(iEns) = nEl(iEns)/nEl(1)
  end do

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

  LZH(:) = 0d0
  EH(:,:) = 0d0
  call unrestricted_hartree_individual_energy(nBas,nEns,Pw,P,ERI,LZH,EH)

!------------------------------------------------------------------------
! Individual exchange energy
!------------------------------------------------------------------------

  LZx(:) = 0d0
  Ex(:,:) = 0d0
  call unrestricted_exchange_individual_energy(x_rung,x_DFA,LDA_centered,nEns,wEns,nCC,aCC,nGrid,weight,nBas,ERI,  &
                                               Pw,rhow,drhow,P,rho,drho,Cx_choice,doNcentered,LZx,Ex)

!------------------------------------------------------------------------
! Individual correlation energy
!------------------------------------------------------------------------

  LZc(:) = 0d0
  Ec(:,:) = 0d0
  call unrestricted_correlation_individual_energy(c_rung,c_DFA,LDA_centered,nEns,wEns,nGrid,weight, & 
                                                  rhow,drhow,rho,drho,LZc,Ec)

!------------------------------------------------------------------------
! Compute auxiliary energies
!------------------------------------------------------------------------
 
  call unrestricted_auxiliary_energy(nBas,nEns,eKS,occnum,Eaux)

!------------------------------------------------------------------------
! Compute derivative discontinuities
!------------------------------------------------------------------------

  do ispin=1,nspin 
    call unrestricted_exchange_derivative_discontinuity(x_rung,x_DFA,nEns,wEns,nCC,aCC,nGrid,weight, &
                                                        rhow(:,ispin),drhow(:,:,ispin),Cx_choice,doNcentered,ExDD(ispin,:))
  end do

  call unrestricted_correlation_derivative_discontinuity(c_rung,c_DFA,nEns,wEns,nGrid,weight,rhow,drhow,EcDD)

! Scaling derivative discontinuity for N-centered ensembles

  if(doNcentered) then

     do iEns=1,nEns
       ExDD(:,iEns) = (1d0 - kappa(iEns))*ExDD(:,iEns)
       EcDD(:,iEns) = (1d0 - kappa(iEns))*EcDD(:,iEns)
     end do

  end if

!------------------------------------------------------------------------
! Total energy
!------------------------------------------------------------------------

  if(doNcentered) then

    do iEns=1,nEns
      E(iEns) = sum(Eaux(:,iEns))                                     & 
              + kappa(iEns)*(sum(LZH(:)) + sum(LZx(:)) + sum(LZc(:))) &
              + sum(ExDD(:,iEns)) + sum(EcDD(:,iEns))
    end do

  else 

    do iEns=1,nEns
      E(iEns) = sum(Eaux(:,iEns))                       & 
              + sum(LZH(:)) + sum(LZx(:)) + sum(LZc(:)) &
              + sum(ExDD(:,iEns)) + sum(EcDD(:,iEns))
    end do

  end if
  
! do iEns=1,nEns
!   E(iEns)   = sum(ET(:,iEns)) + sum(EV(:,iEns))                     & 
!             + sum(EH(:,iEns)) + sum(Ex(:,iEns)) + sum(Ec(:,iEns))   & 
!             + sum(LZH(:)) + sum(LZx(:)) + sum(LZc(:)) &
!             + sum(ExDD(:,iEns)) + sum(EcDD(:,iEns))
! end do

!------------------------------------------------------------------------
! Excitation energies
!------------------------------------------------------------------------

  do iEns=1,nEns

    Om(iEns) = E(iEns) - E(1)

    OmH(iEns)   = sum(EH(:,iEns)) - sum(EH(:,1)) 
    Omx(iEns)   = sum(Ex(:,iEns)) - sum(Ex(:,1))
    Omc(iEns)   = sum(Ec(:,iEns)) - sum(Ec(:,1))

    Omaux(iEns) = sum(Eaux(:,iEns)) - sum(Eaux(:,1))

    OmxDD(iEns) = sum(ExDD(:,iEns)) - sum(ExDD(:,1))
    OmcDD(iEns) = sum(EcDD(:,iEns)) - sum(EcDD(:,1))

  end do

  if(doNcentered) then

    do iEns=1,nEns
      OmH(iEns)   = OmH(iEns) + (kappa(iEns) - kappa(1))*sum(LZH(:))
      Omx(iEns)   = Omx(iEns) + (kappa(iEns) - kappa(1))*sum(LZx(:))
      Omc(iEns)   = Omc(iEns) + (kappa(iEns) - kappa(1))*sum(LZc(:))

      Omaux(iEns) = Omaux(iEns) &
                  + (kappa(iEns) - kappa(1))*(sum(LZH(:)) + sum(LZx(:)) + sum(LZc(:))) 

    end do

  end if

!------------------------------------------------------------------------
! Dump results
!------------------------------------------------------------------------

  call print_unrestricted_individual_energy(nEns,ENuc,Ew,ET,EV,EH,Ex,Ec,Eaux,LZH,LZx,LZc,ExDD,EcDD,E, & 
                                            Om,OmH,Omx,Omc,Omaux,OmxDD,OmcDD)

end subroutine unrestricted_individual_energy
