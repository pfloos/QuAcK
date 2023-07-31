subroutine UGTpp_self_energy_diag(eta,nBas,nC,nO,nV,nR,nHaa,nHab,nHbb,nPaa,nPab,nPbb,e,Om1aa,Om1ab,Om1bb,&
                                  rho1aa,rho1ab,rho1bb,Om2aa,Om2ab,Om2bb,rho2aa,rho2ab,rho2bb,EcGM,SigT,Z)

! Compute diagonal of the correlation part of the T-matrix self-energy

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC(nspin)
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nV(nspin)
  integer,intent(in)            :: nR(nspin)
  integer,intent(in)            :: nHaa,nHab,nHbb
  integer,intent(in)            :: nPaa,nPab,nPbb
  double precision,intent(in)   :: e(nBas,nspin)
  double precision,intent(in)   :: Om1aa(nPaa),Om1ab(nPab),Om1bb(nPbb)
  double precision,intent(in)   :: rho1aa(nBas,nBas,nPaa),rho1ab(nBas,nBas,nPab)
  double precision,intent(in)   :: rho1bb(nBas,nBas,nPbb)
  double precision,intent(in)   :: Om2aa(nHaa),Om2ab(nHab),Om2bb(nHbb)
  double precision,intent(in)   :: rho2aa(nBas,nBas,nHaa),rho2ab(nBas,nBas,nHab)
  double precision,intent(in)   :: rho2bb(nBas,nBas,nHbb)

! Local variables

  integer                       :: i,j,a,b,p,cd,kl
  double precision              :: num,eps

! Output variables

  double precision,intent(inout)  :: EcGM(nspin)
  double precision,intent(inout)  :: SigT(nBas,nspin)
  double precision,intent(inout)  :: Z(nBas,nspin)

! Initialization

  EcGM(:)   = 0d0
  SigT(:,:) = 0d0
  Z(:,:)    = 0d0

!----------------------------------------------
! Occupied part of the T-matrix self-energy 
!----------------------------------------------

! spin up part 

  do p=nC(1)+1,nBas-nR(1)
    do i=nC(1)+1,nO(1)
      do cd=1,nPaa
        eps = e(p,1) + e(i,1) - Om1aa(cd)
        num = rho1aa(p,i,cd)**2
        SigT(p,1) = SigT(p,1) + num*eps/(eps**2 + eta**2)
        Z(p,1) = Z(p,1) - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2
      end do
    end do

    do i=nC(2)+1,nO(2)
      do cd=1,nPab
        eps = e(p,1) + e(i,1) - Om1ab(cd)
        num = rho1ab(p,i,cd)**2
        SigT(p,1) = SigT(p,1) + num*eps/(eps**2 + eta**2)
        Z(p,1) = Z(p,1) - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2
      end do
    end do  
  end do
 
! spin down part
  
  do p=nC(2)+1,nBas-nR(2)
    do i=nC(2)+1,nO(2)
      do cd=1,nPbb
        eps = e(p,2) + e(i,2) - Om1bb(cd)
        num = rho1bb(p,i,cd)**2
        SigT(p,2) = SigT(p,2) + num*eps/(eps**2 + eta**2)
        Z(p,2) = Z(p,2) - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2
      end do
    end do
 
    do i=nC(2)+1,nO(2)
      do cd=1,nPab
        eps = e(p,2) + e(i,2) - Om1ab(cd)
        num = rho1ab(p,i,cd)**2
        SigT(p,2) = SigT(p,2) + num*eps/(eps**2 + eta**2)
        Z(p,2) = Z(p,2) - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2
      end do 
    end do  
  end do
  
!----------------------------------------------
! Virtual part of the T-matrix self-energy
!----------------------------------------------
  
! spin up part
  
  do p=nC(1)+1,nBas-nR(1)
    do a=nO(1)+1,nBas-nR(1)
      do kl=1,nHaa
        eps = e(p,1) + e(a,1) - Om2aa(kl)
        num = rho2aa(p,a,kl)**2
        SigT(p,1) = SigT(p,1) + num*eps/(eps**2 + eta**2)
        Z(p,1) = Z(p,1) - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2
      end do
    end do 
      
    do a=nO(1)+1,nBas-nR(1)
      do kl=1,nHab
        eps = e(p,1) + e(a,1) - Om2ab(kl)
        num = rho2ab(p,a,kl)**2
        SigT(p,1) = SigT(p,1) + num*eps/(eps**2 + eta**2)
        Z(p,1) = Z(p,1) - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2
      end do
    end do   
  end do

! spin down part

  do p=nC(2)+1,nBas-nR(2)
    do a=nO(2)+1,nBas-nR(2)
      do kl=1,nHbb
        eps = e(p,2) + e(a,2) - Om2bb(kl)
        num = rho2bb(p,a,kl)**2
        SigT(p,2) = SigT(p,2) + num*eps/(eps**2 + eta**2)
        Z(p,2) = Z(p,2) - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2
      end do
    end do

    do a=nO(2)+1,nBas-nR(2)
      do kl=1,nHab
        eps = e(p,2) + e(a,2) - Om2ab(kl)
        num = rho2ab(p,a,kl)**2
        SigT(p,2) = SigT(p,2) + num*eps/(eps**2 + eta**2)
        Z(p,2) = Z(p,2) - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2
      end do
    end do 
  end do

  Z(:,:) = 1d0/(1d0 - Z(:,:))

!----------------------------------------------
! Galitskii-Migdal correlation energy
!----------------------------------------------

! spin up part

  do i=nC(1)+1,nO(1)
    do j=nC(1)+1,nO(1)
      do cd=1,nPaa
        eps = e(i,1) + e(j,1) - Om1aa(cd)
        EcGM(1) = EcGM(1) + rho1aa(i,j,cd)*rho1aa(i,j,cd)*eps/(eps**2 + eta**2)
      end do
    end do
  end do 
  
  do i=nC(1)+1,nO(1)
    do j=nC(2)+1,nO(2)
      do cd=1,nPab
        eps = e(i,1) + e(j,1) - Om1ab(cd)
        EcGM(1) = EcGM(1) + rho1ab(i,j,cd)*rho1ab(i,j,cd)*eps/(eps**2 + eta**2)
      end do
    end do
  end do

  do a=nO(1)+1,nBas-nR(1)
    do b=nO(1)+1,nBas-nR(1)
      do kl=1,nHaa
        eps = e(a,1) + e(b,1) - Om2aa(kl)
        EcGM(1) = EcGM(1) - rho2aa(a,b,kl)*rho2aa(a,b,kl)*eps/(eps**2 + eta**2)
      end do
    end do
  end do

  do a=nO(1)+1,nBas-nR(1)
    do b=nO(1)+1,nBas-nR(1)
      do kl=1,nHab
        eps = e(a,1) + e(b,1) - Om2ab(kl)
        EcGM(1) = EcGM(1) - rho2ab(a,b,kl)*rho2ab(a,b,kl)*eps/(eps**2 + eta**2)
      end do
    end do
  end do

! spin down part 

  do i=nC(2)+1,nO(2)
    do j=nC(2)+1,nO(2)
      do cd=1,nPbb
        eps = e(i,2) + e(j,2) - Om1bb(cd)
        EcGM(2) = EcGM(2) + rho1bb(i,j,cd)*rho1bb(i,j,cd)*eps/(eps**2 + eta**2)
      end do
    end do
  end do  
  
  do i=nC(1)+1,nO(1)
    do j=nC(2)+1,nO(2)
      do cd=1,nPab
        eps = e(i,2) + e(j,2) - Om1ab(cd)
        EcGM(2) = EcGM(2) + rho1ab(i,j,cd)*rho1ab(i,j,cd)*eps/(eps**2 + eta**2)
      end do
    end do
  end do

  do a=nO(1)+1,nBas-nR(1)
    do b=nO(2)+1,nBas-nR(2)
      do kl=1,nHab
        eps = e(a,2) + e(b,2) - Om2ab(kl)
        EcGM(2) = EcGM(2) - rho2ab(a,b,kl)*rho2ab(a,b,kl)*eps/(eps**2 + eta**2)
      end do
    end do
  end do

  do a=nO(2)+1,nBas-nR(2)
    do b=nO(2)+1,nBas-nR(2)
      do kl=1,nHbb
        eps = e(a,2) + e(b,2) - Om2bb(kl)
        EcGM(2) = EcGM(2) - rho2bb(a,b,kl)*rho2bb(a,b,kl)*eps/(eps**2 + eta**2)
      end do
    end do
  end do

end subroutine 
