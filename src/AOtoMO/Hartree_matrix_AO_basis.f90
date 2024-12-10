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
  integer                       :: nunu, lala, nula, lasi, numu
  integer*8                     :: nunununu, nunulala, nununula, nunulasi
  integer*8                     :: lalanunu, lasinunu, numulala, lalanumu
  integer*8                     :: numunula, numulasi, lasinumu, nununumu
  integer*8                     :: munu0, munu
  integer*8                     :: sila0, sila
  integer*8                     :: munulasi0, munulasi


  do nu = 1, nBas

    nunu = (nu * (nu - 1)) / 2 + nu
    nunununu = (nunu * (nunu - 1)) / 2 + nunu
    H(nu,nu) = P(nu,nu) * ERI_chem(nunununu)

    do la = 1, nu-1

      lala = (la * (la - 1)) / 2 + la
      nunulala = (nunu * (nunu - 1)) / 2 + lala
      H(nu,nu) = H(nu,nu) + P(la,la) * ERI_chem(nunulala)

      nula = (nu * (nu - 1)) / 2 + la
      nununula = (nunu * (nunu - 1)) / 2 + nula
      H(nu,nu) = H(nu,nu) + 2.d0 * P(la,nu) * ERI_chem(nununula)

      do si = 1, la - 1
        lasi = (la * (la - 1)) / 2 + si
        nunulasi = (nunu * (nunu - 1)) / 2 + lasi
        H(nu,nu) = H(nu,nu) + 2.d0 * P(si,la) * ERI_chem(nunulasi)
      enddo
    enddo

    do la = nu + 1, nBas

      lala = (la * (la - 1)) / 2 + la
      lalanunu = (lala * (lala - 1)) / 2 + nunu
      H(nu,nu) = H(nu,nu) + P(la,la) * ERI_chem(lalanunu)

      do si = 1, la - 1
        lasi = (la * (la - 1)) / 2 + si
        lasinunu = (lasi * (lasi - 1)) / 2 + nunu
        H(nu,nu) = H(nu,nu) + 2.d0 * P(si,la) * ERI_chem(lasinunu)
      enddo
    enddo

    do mu = 1, nu - 1

      numu = (nu * (nu - 1)) / 2 + mu
      nununumu = (nunu * (nunu - 1)) / 2 + numu
      H(mu,nu) = p(nu,nu) * ERI_chem(nununumu)

      do la = 1, nu - 1
        lala = (la * (la - 1)) / 2 + la
        numulala = (numu * (numu - 1)) / 2 + lala
        H(mu,nu) = H(mu,nu) + p(la,la) * ERI_chem(numulala)
      enddo

      do la = nu + 1, nBas
        lala = (la * (la - 1)) / 2 + la
        lalanumu = (lala * (lala - 1)) / 2 + numu
        H(mu,nu) = H(mu,nu) + p(la,la) * ERI_chem(lalanumu)
      enddo

      do la = 1, mu
        nula = (nu * (nu - 1)) / 2 + la
        numunula = (numu * (numu - 1)) / 2 + nula
        H(mu,nu) = H(mu,nu) + 2.d0 * P(la,nu) * ERI_chem(numunula)
      enddo

      do la = mu + 1, nu - 1
        nula = (nu * (nu - 1)) / 2 + la
        numunula = (nula * (nula - 1)) / 2 + numu
        H(mu,nu) = H(mu,nu) + 2.d0 * P(la,nu) * ERI_chem(numunula)
      enddo

      do la = 2, nu - 1
        do si = 1, la - 1
          lasi = (la * (la - 1)) / 2 + si
          numulasi = (numu * (numu - 1)) / 2 + lasi
          H(mu,nu) = H(mu,nu) + 2.d0 * P(si,la) * ERI_chem(numulasi)
        enddo
      enddo

      do la = nu + 1, nBas
        do si = 1, la - 1
          lasi = (la * (la - 1)) / 2 + si
          lasinumu = (lasi * (lasi - 1)) / 2 + numu
          H(mu,nu) = H(mu,nu) + 2.d0 * P(si,la) * ERI_chem(lasinumu)
        enddo
      enddo

    enddo ! mu
  enddo ! nu


  do nu = 1, nBas
    do mu = nu+1, nBas
      H(mu,nu) = H(nu,mu)
    enddo
  enddo

  return
end subroutine 

! ---

