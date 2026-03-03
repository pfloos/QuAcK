subroutine CVS_UGW_SRG_self_energy_diag(flow,nBas,nC,nO,nV,nR,nSt,nCVS,nFC,occupations,virtuals,e,Om,rho,EcGM,SigC,Z)

! Compute diagonal of the correlation part of the self-energy

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: flow
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC(nspin)
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nV(nspin)
  integer,intent(in)            :: nR(nspin)
  integer,intent(in)            :: nSt
  double precision,intent(in)   :: e(nBas,nspin)
  double precision,intent(in)   :: Om(nSt)
  double precision,intent(in)   :: rho(nBas,nBas,nSt,nspin)
  
  integer,intent(in)            :: nCVS(nspin),nFC(nspin)
  integer,intent(in)            :: occupations(maxval(nO-nFC),nspin)
  integer,intent(in)            :: virtuals(nBas-minval(nO),nspin)

! Local variables

  integer                       :: i,a,p,m,ispin
  double precision              :: num,eps
  double precision              :: s
  double precision              :: Dpim,Dpam,Diam

! Output variables

  double precision,intent(out)  :: SigC(nBas,nspin)
  double precision,intent(out)  :: Z(nBas,nspin)
  double precision              :: EcGM(nspin)

! Initialize 

  SigC(:,:) = 0d0
  Z(:,:)   = 0d0
  EcGM(:)  = 0d0
  s = flow
  

  do ispin=1,nspin
    ! Occupied part of the correlation self-energy

    do p=1,nBas
      do i=1,nO(ispin)-nFC(ispin)
        do m=1,nSt
           Dpim = e(p,ispin) - e(occupations(i,ispin),ispin) + Om(m)
           SigC(p,ispin) = SigC(p,ispin) + rho(p,occupations(i,ispin),m,ispin)**2&
                                         * (1d0-exp(-2d0*s*Dpim*Dpim))/Dpim
        end do
      end do
    end do

    ! Virtual part of the correlation self-energy

    do p=1,nBas
      do a=1+nCVS(ispin),nBas - nO(ispin)
        do m=1,nSt
          Dpam = e(p,ispin) - e(virtuals(a,ispin),ispin) - Om(m)
          SigC(p,ispin) = SigC(p,ispin) + rho(p,virtuals(a,ispin),m,ispin)**2&
                                        * (1d0-exp(-2d0*s*Dpam*Dpam))/Dpam
        end do
      end do
    end do

    ! GM correlation energy

    do i=1,nO(ispin)-nFC(ispin)
      do a=nCVS(ispin)+1,nBas-nO(ispin)
        do m=1,nSt
          Diam = e(virtuals(a,ispin),ispin) - e(occupations(i,ispin),ispin) + Om(m)
          EcGM(ispin) = EcGM(ispin)& 
                      - rho(virtuals(a,ispin),occupations(i,ispin),m,ispin)&
                      * rho(virtuals(a,ispin),occupations(i,ispin),m,ispin)&
                      * (1d0-exp(-2d0*s*Diam*Diam))/Diam 
        end do
      end do
    end do
  end do

! Compute renormalization factor from derivative 

  Z(:,:) = 1d0/(1d0 - Z(:,:))
end subroutine 
