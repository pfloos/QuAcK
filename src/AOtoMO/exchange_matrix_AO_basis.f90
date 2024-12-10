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
  integer*8                     :: nulasimu, lanumusi, lanusimu, simunula
  integer*8                     :: simulanu



  do nu = 1, nBas

    nunu = shiftr(nu * (nu - 1), 1) + nu
    nunununu = shiftr(nunu * (nunu - 1), 1) + nunu
    K(nu,nu) = -P(nu,nu) * ERI_chem(nunununu)

    do la = 1, nu - 1
      nula = shiftr(nu * (nu - 1), 1) + la
      nulanula = shiftr(nula * (nula - 1), 1) + nula
      K(nu,nu) = K(nu,nu) - P(la,la) * ERI_chem(nulanula)
    enddo

    do la = nu + 1, nBas
      lanu = shiftr(la * (la - 1), 1) + nu
      lanulanu = shiftr(lanu * (lanu - 1), 1) + lanu
      K(nu,nu) = K(nu,nu) - P(la,la) * ERI_chem(lanulanu)
    enddo

    do la = 1, nu
      nula = shiftr(nu * (nu - 1), 1) + la
      do si = 1, la - 1
        nusi = shiftr(nu * (nu - 1), 1) + si
        nulanusi = shiftr(nula * (nula - 1), 1) + nusi
        K(nu,nu) = K(nu,nu) - 2.d0 * P(si,la) * ERI_chem(nulanusi)
      enddo
    enddo

    do la = nu + 1, nBas
      lanu = shiftr(la * (la - 1), 1) + nu
      do si = 1, nu
        nusi = shiftr(nu * (nu - 1), 1) + si
        lanunusi = shiftr(lanu * (lanu - 1), 1) + nusi
        K(nu,nu) = K(nu,nu) - 2.d0 * P(si,la) * ERI_chem(lanunusi)
      enddo
      do si = nu + 1, la - 1
        sinu = shiftr(si * (si - 1), 1) + nu
        lanusinu = shiftr(lanu * (lanu - 1), 1) + sinu
        K(nu,nu) = K(nu,nu) - 2.d0 * P(si,la) * ERI_chem(lanusinu)
      enddo
    enddo


    do mu = 1, nu - 1

      numu = shiftr(nu * (nu - 1), 1) + mu
      mumu = shiftr(mu * (mu - 1), 1) + mu
      numumumu = shiftr(numu * (numu - 1), 1) + mumu
      K(mu,nu) = - P(mu,mu) * ERI_chem(numumumu)

      do la = 1, mu - 1
        mula = shiftr(mu * (mu - 1), 1) + la
        nula = shiftr(nu * (nu - 1), 1) + la
        nulamula = shiftr(nula * (nula - 1), 1) + mula
        K(mu,nu) = K(mu,nu) - P(la,la) * ERI_chem(nulamula)
      enddo
      do la = mu + 1, nu
        lamu = shiftr(la * (la - 1), 1) + mu
        nula = shiftr(nu * (nu - 1), 1) + la
        nulalamu = shiftr(nula * (nula - 1), 1) + lamu
        K(mu,nu) = K(mu,nu) - P(la,la) * ERI_chem(nulalamu)
      enddo
      do la = nu + 1, nBas
        lamu = shiftr(la * (la - 1), 1) + mu
        lanu = shiftr(la * (la - 1), 1) + nu
        lanulamu = shiftr(lanu * (lanu - 1), 1) + lamu
        K(mu,nu) = K(mu,nu) - P(la,la) * ERI_chem(lanulamu)
      enddo

      do la = 1, mu
        nula = shiftr(nu * (nu - 1), 1) + la
        do si = 1, la - 1
          musi = shiftr(mu * (mu - 1), 1) + si
          nulamusi = shiftr(nula * (nula - 1), 1) + musi
          K(mu,nu) = K(mu,nu) - P(si,la) * ERI_chem(nulamusi)
        enddo
      enddo
      do la = mu + 1, nu
        nula = shiftr(nu * (nu - 1), 1) + la
        do si = 1, mu
          musi = shiftr(mu * (mu - 1), 1) + si
          nulamusi = shiftr(nula * (nula - 1), 1) + musi
          K(mu,nu) = K(mu,nu) - P(si,la) * ERI_chem(nulamusi)
        enddo
        do si = mu + 1, la - 1
          simu = shiftr(si * (si - 1), 1) + mu
          nulasimu = shiftr(nula * (nula - 1), 1) + simu
          K(mu,nu) = K(mu,nu) - P(si,la) * ERI_chem(nulasimu)
        enddo
      enddo
      do la = nu + 1, nBas
        lanu = shiftr(la * (la - 1), 1) + nu
        do si = 1, mu
          musi = shiftr(mu * (mu - 1), 1) + si
          lanumusi = shiftr(lanu * (lanu - 1), 1) + musi
          K(mu,nu) = K(mu,nu) - P(si,la) * ERI_chem(lanumusi)
        enddo
        do si = mu + 1, la - 1
          simu = shiftr(si * (si - 1), 1) + mu
          lanusimu = shiftr(lanu * (lanu - 1), 1) + simu
          K(mu,nu) = K(mu,nu) - P(si,la) * ERI_chem(lanusimu)
        enddo
      enddo

      do la = 1, mu
        nula = shiftr(nu * (nu - 1), 1) + la
        do si = la + 1, mu
          musi = shiftr(mu * (mu - 1) , 1) + si
          nulamusi = shiftr(nula * (nula - 1), 1) + musi
          K(mu,nu) = K(mu,nu) - P(si,la) * ERI_chem(nulamusi)
        enddo
        do si = mu + 1, nu - 1
          simu = shiftr(si * (si - 1), 1) + mu
          nulasimu = shiftr(nula * (nula - 1), 1) + simu
          K(mu,nu) = K(mu,nu) - P(si,la) * ERI_chem(nulasimu)
        enddo
        do si = nu, nBas
          simu = shiftr(si * (si - 1), 1) + mu
          simunula = shiftr(simu * (simu - 1), 1) + nula
          K(mu,nu) = K(mu,nu) - P(si,la) * ERI_chem(simunula)
        enddo
      enddo
      do la = mu + 1, nu
        nula = shiftr(nu * (nu - 1), 1) + la
        do si = la + 1, nu
          simu = shiftr(si * (si - 1), 1) + mu
          nulasimu = shiftr(nula * (nula - 1), 1) + simu
          K(mu,nu) = K(mu,nu) - P(si,la) * ERI_chem(nulasimu)
        enddo
        do si = nu + 1, nBas
          simu = shiftr(si * (si - 1), 1) + mu
          simunula = shiftr(simu * (simu - 1), 1) + nula
          K(mu,nu) = K(mu,nu) - P(si,la) * ERI_chem(simunula)
        enddo
      enddo
      do la = nu + 1, nBas
        lanu = shiftr(la * (la - 1), 1) + nu
        do si = la + 1, nBas
          simu = shiftr(si * (si - 1), 1) + mu
          simulanu = shiftr(simu * (simu - 1), 1) + lanu
          K(mu,nu) = K(mu,nu) - P(si,la) * ERI_chem(simulanu)
        enddo
      enddo

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

