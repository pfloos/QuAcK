subroutine CVS_UGW_self_energy_diag(eta,nBas,nC,nO,nV,nR,nSt,nCVS,nFC,occupations,virtuals,e,Om,rho,EcGM,Sig,Z)

! Compute diagonal of the correlation part of the self-energy

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: eta
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

  integer                       :: i,a,p,m
  double precision              :: num,eps

! Output variables

  double precision,intent(out)  :: Sig(nBas,nspin)
  double precision,intent(out)  :: Z(nBas,nspin)
  double precision              :: EcGM(nspin)

! Initialize 

  Sig(:,:) = 0d0
  Z(:,:)   = 0d0
  EcGM(:)  = 0d0

!--------------!
! Spin-up part !
!--------------!

  ! Occupied part of the correlation self-energy

  do p=1,nBas
    do i=1,nO(1)-nFC(1)
      do m=1,nSt
        eps = e(p,1) - e(occupations(i,1),1) + Om(m)
        num = rho(p,occupations(i,1),m,1)**2
        Sig(p,1) = Sig(p,1) + num*eps/(eps**2 + eta**2)
        Z(p,1)   = Z(p,1)   - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2
      end do
    end do
  end do

  ! Virtual part of the correlation self-energy

  do p=1,nBas
    do a=1+nCVS(1),nBas - nO(1)
      do m=1,nSt
        eps = e(p,1) - e(virtuals(a,1),1) - Om(m)
        num = rho(p,virtuals(a,1),m,1)**2
        Sig(p,1) = Sig(p,1) + num*eps/(eps**2 + eta**2)
        Z(p,1)   = Z(p,1)   - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2
      end do
    end do
  end do

  ! GM correlation energy

  do i=1,nO(1)-nFC(1)
    do a=nCVS(1)+1,nBas-nO(1)
      do m=1,nSt
        eps = e(virtuals(a,1),1) - e(occupations(i,1),1) + Om(m)  
        num = rho(virtuals(a,1),occupations(i,1),m,1)**2
        EcGM(1) = EcGM(1) - num*eps/(eps**2 + eta**2)
      end do
    end do
  end do

!----------------!
! Spin-down part !
!----------------!

  ! Occupied part of the correlation self-energy

  do p=1,nBas
    do i=1,nO(2)-nFC(2)
      do m=1,nSt
        eps = e(p,2) - e(occupations(i,2),2) + Om(m)
        num = rho(p,occupations(i,2),m,2)**2
        Sig(p,2) = Sig(p,2) + num*eps/(eps**2 + eta**2)
        Z(p,2)   = Z(p,2)   - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2
      end do
    end do
  end do

  ! Virtual part of the correlation self-energy

  do p=1,nBas
    do a=nCVS(2)+1,nBas-nO(2)
      do m=1,nSt
        eps = e(p,2) - e(virtuals(a,2),2) - Om(m)
        num = rho(p,virtuals(a,2),m,2)**2
        Sig(p,2) = Sig(p,2) + num*eps/(eps**2 + eta**2)
        Z(p,2)   = Z(p,2)   - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2
      end do
    end do
  end do

  ! GM correlation energy

  do i=1,nO(2)-nFC(2)
    do a=nCVS(2)+1,nBas-nO(2)
      do m=1,nSt
        eps = e(virtuals(a,2),2) - e(occupations(i,2),2) + Om(m)
        num = rho(virtuals(a,2),occupations(i,2),m,2)**2
        EcGM(2) = EcGM(2) - num*eps/(eps**2 + eta**2)
      end do
    end do
  end do

! Compute renormalization factor from derivative 

  Z(:,:) = 1d0/(1d0 - Z(:,:))

end subroutine 
