subroutine UGTpp_renormalization_factor(eta,nBas,nC,nO,nV,nR,nHaa,nHab,nHbb,nPaa,nPab,nPbb,e,Om1aa,Om1ab, &
                                        Om1bb,rho1aa,rho1ab,rho1bb,Om2aa,Om2ab,Om2bb,rho2aa,rho2ab,rho2bb,Z)

! Compute renormalization factor of the T-matrix self-energy

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas,nC(nspin),nO(nspin),nV(nspin),nR(nspin)
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

  integer                       :: i,a,p,cd,kl
  double precision              :: eps

! Output variables

  double precision,intent(out)  :: Z(nBas,nspin)
  
!spin up part

! Occupied part of the T-matrix self-energy 

  do p=nC(1)+1,nBas-nR(1)
    do i=nC(1)+1,nO(1)
      do cd=1,nPaa
        eps  = e(p,1) + e(i,1) - Om1aa(cd)
        Z(p,1) = Z(p,1) - rho1aa(p,i,cd)**2*(eps/(eps**2 + eta**2))**2
      enddo
    enddo
   
    do i=nC(1)+1,nO(1)
      do cd=1,nPab
        eps  = e(p,1) + e(i,1) - Om1ab(cd)
        Z(p,1) = Z(p,1) - rho1ab(p,i,cd)**2*(eps/(eps**2 + eta**2))**2
      end do
    end do
  enddo 

! Virtual part of the T-matrix self-energy

  do p=nC(1)+1,nBas-nR(1)
    do a=nO(1)+1,nBas-nR(1)
      do kl=1,nHaa
        eps  = e(p,1) + e(a,1) - Om2aa(kl)
        Z(p,1) = Z(p,1) - rho2aa(p,a,kl)**2*(eps/(eps**2 + eta**2))**2
      enddo
    enddo

    do a=nO(2)+1,nBas-nR(2)
      do kl=1,nHab
        eps  = e(p,1) + e(a,1) - Om2ab(kl)
        Z(p,1) = Z(p,1) - rho2ab(p,a,kl)**2*(eps/(eps**2 + eta**2))**2
      enddo
    enddo
  enddo

!spin down part

! Occupied part of the T-matrix self-energy

  do p=nC(2)+1,nBas-nR(2)
    do i=nC(2)+1,nO(2)
      do cd=1,nPbb
        eps  = e(p,2) + e(i,2) - Om1bb(cd)
        Z(p,2) = Z(p,2) - rho1bb(p,i,cd)**2*(eps/(eps**2 + eta**2))**2
      enddo
    enddo

    do i=nC(1)+1,nO(1)
      do cd=1,nPab
        eps  = e(p,2) + e(i,2) - Om1ab(cd)
        Z(p,2) = Z(p,2) - rho1ab(p,i,cd)**2*(eps/(eps**2 + eta**2))**2
      enddo
    enddo
  enddo

! Virtual part of the T-matrix self-energy

  do p=nC(2)+1,nBas-nR(2)
    do a=nO(2)+1,nBas-nR(2)
      do kl=1,nHbb
        eps  = e(p,2) + e(a,2) - Om2bb(kl)
        Z(p,2) = Z(p,2) - rho2bb(p,a,kl)**2*(eps/(eps**2 + eta**2))**2
      enddo
    enddo

    do a=nO(1)+1,nBas-nR(1)
      do kl=1,nHab
        eps  = e(p,2) + e(a,2) - Om2ab(kl)
        Z(p,2) = Z(p,2) - rho2ab(p,a,kl)**2*(eps/(eps**2 + eta**2))**2
      enddo
    enddo
  enddo

end subroutine 