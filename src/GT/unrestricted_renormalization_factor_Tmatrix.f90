subroutine unrestricted_renormalization_factor_Tmatrix(ispin,eta,nBas,nC,nO,nV,nR,nH,nP,e,Omega1,rho1,Omega2,rho2,Z)

! Compute renormalization factor of the T-matrix self-energy

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas,nC(nspin),nO(nspin),nV(nspin),nR(nspin)
  integer,intent(in)            :: nH,ispin
  integer,intent(in)            :: nP
  double precision,intent(in)   :: e(nBas,nspin)
  double precision,intent(in)   :: Omega1(nP)
  double precision,intent(in)   :: rho1(nBas,nBas,nP,nspin)
  double precision,intent(in)   :: Omega2(nH)
  double precision,intent(in)   :: rho2(nBas,nBas,nH,nspin)

! Local variables

  integer                       :: i,a,p,cd,kl
  double precision              :: eps

! Output variables

  double precision,intent(out)  :: Z(nBas,nspin)

! Occupied part of the T-matrix self-energy 
  
  if(ispin==1) then

    do p=nC(1)+1,nBas-nR(1)
      do i=nC(1)+1,nO(1)
        do cd=1,nP
          eps  = e(p,1) + e(i,1) - Omega1(cd)
          Z(p,1) = Z(p,1) - rho1(p,i,cd,1)**2*(eps/(eps**2 + eta**2))**2
        enddo
      enddo
    enddo 

! Virtual part of the T-matrix self-energy

    do p=nC(1)+1,nBas-nR(1)
      do a=nO(1)+1,nBas-nR(1)
        do kl=1,nH
          eps  = e(p,1) + e(a,1) - Omega2(kl)
          Z(p,1) = Z(p,1) - rho2(p,a,kl,1)**2*(eps/(eps**2 + eta**2))**2
        enddo
      enddo
    enddo

  end if

! Occupied part of the T-matrix self-energy 
 
  if(ispin==2) then

    do p=nC(2)+1,nBas-nR(2)
      do i=nC(2)+1,nO(2)
        do cd=1,nP
          eps  = e(p,2) + e(i,2) - Omega1(cd)
          Z(p,2) = Z(p,2) - rho1(p,i,cd,2)**2*(eps/(eps**2 + eta**2))**2
        enddo
      enddo
    enddo

! Virtual part of the T-matrix self-energy

    do p=nC(2)+1,nBas-nR(2)
      do a=nO(2)+1,nBas-nR(2)
        do kl=1,nH
          eps  = e(p,2) + e(a,2) - Omega2(kl)
          Z(p,2) = Z(p,2) - rho2(p,a,kl,2)**2*(eps/(eps**2 + eta**2))**2
        enddo
      enddo
    enddo

  end if

end subroutine unrestricted_renormalization_factor_Tmatrix
