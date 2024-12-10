subroutine exchange_matrix_AO_basis(nBas,P,ERI,K)

! Compute exchange matrix in the AO basis

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  double precision,intent(in)   :: P(nBas,nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: mu,nu,la,si

! Output variables

  double precision,intent(out)  :: K(nBas,nBas)

  K = 0d0
  do nu=1,nBas
    do si=1,nBas
      do la=1,nBas
        do mu=1,nBas
          K(mu,nu) = K(mu,nu) - P(la,si)*ERI(mu,la,si,nu)
        end do
      end do
    end do
  end do

end subroutine 

! ---

subroutine exchange_matrix_AO_basis_hpc(nBas, ERI_size, P, ERI_chem, K)

  implicit none

  integer,          intent(in)  :: nBas
  integer*8,        intent(in)  :: ERI_size
  double precision, intent(in)  :: P(nBas,nBas)
  double precision, intent(in)  :: ERI_chem(ERI_size)
  double precision, intent(out) :: K(nBas,nBas)

  integer                       :: mu, nu, la, si
  integer                       :: nunu, nula, lanu, lasi, nusi, sinu
  integer                       :: numu, mumu, mula, lamu, musi, simu
  integer*8                     :: nunununu, nulanula, lanulanu, nulanusi
  integer*8                     :: munulasi, lanunusi, lanusinu, numumumu 
  integer*8                     :: nulamula, nulalamu, lanulamu, nulamusi
  integer*8                     :: nulasimu, lanumusi, lanusimu


  integer*8, external           :: Yoshimine_4ind


  do nu = 1, nBas

    nunu = (nu * (nu - 1)) / 2 + nu
    nunununu = (nunu * (nunu - 1)) / 2 + nunu
    K(nu,nu) = -P(nu,nu) * ERI_chem(nunununu)

    do la = 1, nu - 1
      nula = (nu * (nu - 1)) / 2 + la
      nulanula = (nula * (nula - 1)) / 2 + nula
      K(nu,nu) = K(nu,nu) - P(la,la) * ERI_chem(nulanula)
    enddo

    do la = nu + 1, nBas
      lanu = (la * (la - 1)) / 2 + nu
      lanulanu = (lanu * (lanu - 1)) / 2 + lanu
      K(nu,nu) = K(nu,nu) - P(la,la) * ERI_chem(lanulanu)
    enddo

    do la = 1, nu
      nula = (nu * (nu - 1)) / 2 + la
      do si = 1, la - 1
        nusi = (nu * (nu - 1)) / 2 + si
        nulanusi = (nula * (nula - 1)) / 2 + nusi
        K(nu,nu) = K(nu,nu) - 2.d0 * P(si,la) * ERI_chem(nulanusi)
      enddo
    enddo

    do la = nu + 1, nBas
      lanu = (la * (la - 1)) / 2 + nu
      do si = 1, nu
        nusi = (nu * (nu - 1)) / 2 + si
        lanunusi = (lanu * (lanu - 1)) / 2 + nusi
        K(nu,nu) = K(nu,nu) - 2.d0 * P(si,la) * ERI_chem(lanunusi)
      enddo
      do si = nu + 1, la - 1
        sinu = (si * (si - 1)) / 2 + nu
        lanusinu = (lanu * (lanu - 1)) / 2 + sinu
        K(nu,nu) = K(nu,nu) - 2.d0 * P(si,la) * ERI_chem(lanusinu)
      enddo
    enddo


    do mu = 1, nu - 1

      numu = (nu * (nu - 1)) / 2 + mu
      mumu = (mu * (mu - 1)) / 2 + mu
      numumumu = (numu * (numu - 1)) / 2 + mumu
      K(mu,nu) = - P(mu,mu) * ERI_chem(numumumu)

      do la = 1, mu - 1
        mula = (mu * (mu - 1)) / 2 + la
        nula = (nu * (nu - 1)) / 2 + la
        nulamula = (nula * (nula - 1)) / 2 + mula
        K(mu,nu) = K(mu,nu) - P(la,la) * ERI_chem(nulamula)
      enddo
      do la = mu + 1, nu
        lamu = (la * (la - 1)) / 2 + mu
        nula = (nu * (nu - 1)) / 2 + la
        nulalamu = (nula * (nula - 1)) / 2 + lamu
        K(mu,nu) = K(mu,nu) - P(la,la) * ERI_chem(nulalamu)
      enddo
      do la = nu + 1, nBas
        lamu = (la * (la - 1)) / 2 + mu
        lanu = (la * (la - 1)) / 2 + nu
        lanulamu = (lanu * (lanu - 1)) / 2 + lamu
        K(mu,nu) = K(mu,nu) - P(la,la) * ERI_chem(lanulamu)
      enddo

      do la = 1, mu
        nula = (nu * (nu - 1)) / 2 + la
        do si = 1, la - 1
          musi = (mu * (mu - 1)) / 2 + si
          nulamusi = (nula * (nula - 1)) / 2 + musi
          K(mu,nu) = K(mu,nu) - P(si,la) * ERI_chem(nulamusi)
        enddo
      enddo
      do la = mu + 1, nu
        nula = (nu * (nu - 1)) / 2 + la
        do si = 1, mu
          musi = (mu * (mu - 1)) / 2 + si
          nulamusi = (nula * (nula - 1)) / 2 + musi
          K(mu,nu) = K(mu,nu) - P(si,la) * ERI_chem(nulamusi)
        enddo
        do si = mu + 1, la - 1
          simu = (si * (si - 1)) / 2 + mu
          nulasimu = (nula * (nula - 1)) / 2 + simu
          K(mu,nu) = K(mu,nu) - P(si,la) * ERI_chem(nulasimu)
        enddo
      enddo
      do la = nu + 1, nBas
        lanu = (la * (la - 1)) / 2 + nu
        do si = 1, mu
          musi = (mu * (mu - 1)) / 2 + si
          lanumusi = (lanu * (lanu - 1)) / 2 + musi
          K(mu,nu) = K(mu,nu) - P(si,la) * ERI_chem(lanumusi)
        enddo
        do si = mu + 1, la - 1
          simu = (si * (si - 1)) / 2 + mu
          lanusimu = (lanu * (lanu - 1)) / 2 + simu
          K(mu,nu) = K(mu,nu) - P(si,la) * ERI_chem(lanusimu)
        enddo
      enddo

!TODO
!      do la = 1, mu
!        nula = (nu * (nu - 1)) / 2 + la
!        do si = la + 1, mu
!          musi = (mu * (mu - 1)) / 2 + si
!          nulamusi = (nula * (nula - 1)) / 2 + musi
!          !nulamusi = Yoshimine_4ind(nu, la, si, mu)
!          K(mu,nu) = K(mu,nu) - P(si,la) * ERI_chem(nulamusi)
!        enddo
!        do si = mu + 1, nBas
!          simu = (si * (si - 1)) / 2 + mu
!          nulasimu = (nula * (nula - 1)) / 2 + simu
!          !nulasimu = Yoshimine_4ind(nu, la, si, mu)
!          K(mu,nu) = K(mu,nu) - P(si,la) * ERI_chem(nulasimu)
!        enddo
!      enddo
!      do la = mu + 1, nu
!        nula = (nu * (nu - 1)) / 2 + la
!        do si = la + 1, nu
!          simu = (si * (si - 1)) / 2 + mu
!          nulasimu = (nula * (nula - 1)) / 2 + simu
!          !nulasimu = Yoshimine_4ind(nu, la, si, mu)
!          K(mu,nu) = K(mu,nu) - P(si,la) * ERI_chem(nulasimu)
!        enddo
!        do si = nu + 1, nBas
!          simu = (si * (si - 1)) / 2 + mu
!          munulasi = Yoshimine_4ind(nu, la, si, mu)
!          K(mu,nu) = K(mu,nu) - P(si,la) * ERI_chem(munulasi)
!        enddo
!      enddo
!      do la = nu + 1, nBas
!        lanu = (la * (la - 1)) / 2 + nu
!        do si = la + 1, mu
!          simu = (si * (si - 1)) / 2 + mu
!          munulasi = Yoshimine_4ind(nu, la, si, mu)
!          K(mu,nu) = K(mu,nu) - P(si,la) * ERI_chem(munulasi)
!        enddo
!        do si = mu + 1, nBas
!          musi = (mu * (mu - 1)) / 2 + si
!          munulasi = Yoshimine_4ind(nu, la, si, mu)
!          K(mu,nu) = K(mu,nu) - P(si,la) * ERI_chem(munulasi)
!        enddo
!      enddo

    enddo ! mu
  enddo ! nu


  do nu = 1, nBas
    do mu = nu+1, nBas
      K(mu,nu) = K(nu,mu)
    enddo
  enddo

  return
end subroutine 

! ---

