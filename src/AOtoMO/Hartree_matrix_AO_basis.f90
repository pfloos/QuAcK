subroutine Hartree_matrix_AO_basis(nBas,P,G,H)

! Compute Hartree matrix in the AO basis

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  double precision,intent(in)   :: P(nBas,nBas)
  double precision,intent(in)   :: G(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: mu,nu,la,si

! Output variables

  double precision,intent(out)  :: H(nBas,nBas)

  H(:,:) = 0d0

  do si=1,nBas
    do nu=1,nBas
      do la=1,nBas
        do mu=1,nBas
          H(mu,nu) = H(mu,nu) + P(la,si)*G(mu,la,nu,si)
        end do
      end do
    end do
  end do

end subroutine 

! ---

subroutine Hartree_matrix_AO_basis_hpc(nBas, ERI_size, P, ERI_chem, H)

  implicit none

  integer,          intent(in)  :: nBas
  integer*8,        intent(in)  :: ERI_size
  double precision, intent(in)  :: P(nBas,nBas)
  double precision, intent(in)  :: ERI_chem(ERI_size)
  double precision, intent(out) :: H(nBas,nBas)

  integer                       :: mu, nu, la, si
  integer                       :: nunu, lala, nula, lasi
  integer*8                     :: nunununu, nunulala, nununula, nunulasi, lalanunu, lasinunu
  integer*8                     :: munu0, munu
  integer*8                     :: sila0, sila
  integer*8                     :: munulasi0, munulasi

  integer*8, external           :: Yoshimine_4ind

!  do nu = 1, nBas
!    do mu = 1, nBas
!      H(mu,nu) = 0.d0
!      do si = 1, nBas
!        do la = 1, nBas
!          munulasi = Yoshimine_4ind(mu, nu, la, si)
!          H(mu,nu) = H(mu,nu) + P(la,si) * ERI_chem(munulasi)
!        enddo
!      enddo
!    enddo
!  enddo

!  do nu = 1, nBas
!    do mu = 1, nu
!      H(mu,nu) = 0.d0
!      do si = 1, nBas
!        munulasi = Yoshimine_4ind(mu, nu, si, si)
!        H(mu,nu) = H(mu,nu) + P(si,si) * ERI_chem(munulasi)
!        do la = 1, si-1
!          munulasi = Yoshimine_4ind(mu, nu, la, si)
!          H(mu,nu) = H(mu,nu) + 2.d0 * P(la,si) * ERI_chem(munulasi)
!        enddo
!      enddo
!    enddo
!  enddo



  do nu = 1, nBas

    nunu = (nu * (nu - 1)) / 2 + nu
    nunununu = (nunu * (nunu - 1)) / 2 + nunu
    !nunununu = Yoshimine_4ind(nu, nu, nu, nu)
    H(nu,nu) = P(nu,nu) * ERI_chem(nunununu)

    do la = 1, nu-1

      ! la < nu
      lala = (la * (la - 1)) / 2 + la
      nunulala = (nunu * (nunu - 1)) / 2 + lala
      !nunulala = Yoshimine_4ind(nu, nu, la, la)
      H(nu,nu) = H(nu,nu) + P(la,la) * ERI_chem(nunulala)

      ! la < nu
      nula = (nu * (nu - 1)) / 2 + la
      nununula = (nunu * (nunu - 1)) / 2 + nula
      !nununula = Yoshimine_4ind(nu, nu, la, nu)
      H(nu,nu) = H(nu,nu) + 2.d0 * P(la,nu) * ERI_chem(nununula)

      do si = 1, la - 1
        ! lasi < nunu
        lasi = (la * (la - 1)) / 2 + si
        nunulasi = (nunu * (nunu - 1)) / 2 + lasi
        !nunulasi = Yoshimine_4ind(nu, nu, si, la)
        H(nu,nu) = H(nu,nu) + 2.d0 * P(si,la) * ERI_chem(nunulasi)
      enddo
    enddo

    do la = nu+1, nBas

      ! nu < la
      lala = (la * (la - 1)) / 2 + la
      lalanunu = (lala * (lala - 1)) / 2 + nunu
      !lalanunu = Yoshimine_4ind(nu, nu, la, la)
      H(nu,nu) = H(nu,nu) + P(la,la) * ERI_chem(lalanunu)

      do si = 1, la - 1
        ! nunu < lasi
        lasi = (la * (la - 1)) / 2 + si
        lasinunu = (lasi * (lasi - 1)) / 2 + nunu
        !lasinunu = Yoshimine_4ind(nu, nu, si, la)
        H(nu,nu) = H(nu,nu) + 2.d0 * P(si,la) * ERI_chem(lasinunu)
      enddo
    enddo

    do mu = 1, nu-1
      H(mu,nu) = 0.d0
      do si = 1, nBas
        munulasi = Yoshimine_4ind(mu, nu, si, si)
        H(mu,nu) = H(mu,nu) + P(si,si) * ERI_chem(munulasi)
        do la = 1, si-1
          munulasi = Yoshimine_4ind(mu, nu, la, si)
          H(mu,nu) = H(mu,nu) + 2.d0 * P(la,si) * ERI_chem(munulasi)
        enddo
      enddo
    enddo
  enddo


  do nu = 1, nBas
    do mu = nu+1, nBas
      H(mu,nu) = H(nu,mu)
    enddo
  enddo

  return
end subroutine 

! ---

