subroutine GTpp_self_energy(eta,nBas,nC,nO,nV,nR,nOO,nVV,e,Om1,rho1,Om2,rho2,EcGM,Sig,Z)

! Compute the correlation part of the T-matrix self-energy and the renormalization factor

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nOO
  integer,intent(in)            :: nVV
  double precision,intent(in)   :: e(nBas)
  double precision,intent(in)   :: Om1(nVV)
  double precision,intent(in)   :: rho1(nBas,nBas,nVV)
  double precision,intent(in)   :: Om2(nOO)
  double precision,intent(in)   :: rho2(nBas,nBas,nOO)

! Local variables

  integer                       :: i,j,a,b,p,q,cd,kl
  double precision              :: num,eps

! Output variables

  double precision,intent(inout):: EcGM
  double precision,intent(inout):: Sig(nBas,nBas)
  double precision,intent(inout):: Z(nBas)

!----------------------------------------------
! Occupied part of the T-matrix self-energy 
!----------------------------------------------

  do p=nC+1,nBas-nR
    do q=nC+1,nBas-nR
      do i=nC+1,nO
        do cd=1,nVV

          eps = e(p) + e(i) - Om1(cd)
          num = rho1(p,i,cd)*rho1(q,i,cd)
          Sig(p,q) = Sig(p,q) + num*eps/(eps**2 + eta**2)
          if(p == q) Z(p) = Z(p) - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2

        enddo
      enddo
    enddo
  enddo

!----------------------------------------------
! Virtual part of the T-matrix self-energy
!----------------------------------------------

  do p=nC+1,nBas-nR
    do q=nC+1,nBas-nR
      do a=nO+1,nBas-nR
        do kl=1,nOO

          eps = e(p) + e(a) - Om2(kl)
          num = rho2(p,a,kl)*rho2(q,a,kl)
          Sig(p,q) = Sig(p,q) + num*eps/(eps**2 + eta**2)
          if(p == q) Z(p) = Z(p) - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2

        enddo
      enddo
    enddo
  enddo

!----------------------------------------------
! Galitskii-Migdal correlation energy
!----------------------------------------------

  do i=nC+1,nO
    do j=nC+1,nO
      do cd=1,nVV

        eps = e(i) + e(j) - Om1(cd)
        num = rho1(i,j,cd)*rho1(i,j,cd)
        EcGM = EcGM + num*eps/(eps**2 + eta**2)

      enddo
    enddo
  enddo

  do a=nO+1,nBas-nR
    do b=nO+1,nBas-nR
      do kl=1,nOO

        eps = e(a) + e(b) - Om2(kl)
        num = rho2(a,b,kl)*rho2(a,b,kl)
        EcGM = EcGM - num*eps/(eps**2 + eta**2)

      enddo
    enddo
  enddo

end subroutine 
