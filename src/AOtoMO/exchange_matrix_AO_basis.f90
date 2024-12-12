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
  integer                       :: nunu0, lala0, mumu0
  integer*8                     :: nunununu, nulanula, lanulanu, nulanusi
  integer*8                     :: munulasi, lanunusi, lanusinu, numumumu 
  integer*8                     :: nulamula, nulalamu, lanulamu, nulamusi
  integer*8                     :: nulasimu, lanumusi, lanusimu, simunula
  integer*8                     :: simulanu, nulanula0, lanulanu0

!  integer*8                     :: munusila
!  integer*8, external           :: Yoshimine_4ind
!
!  do nu = 1, nBas
!    do mu = 1, nu
!      K(mu,nu) = 0.d0
!      do la = 1, nBas
!        do si = 1, nBas
!          munusila = Yoshimine_4ind(int(mu, kind=8), &
!                                    int(si, kind=8), &
!                                    int(la, kind=8), &
!                                    int(nu, kind=8))
!          K(mu,nu) = K(mu,nu) - P(si,la) * ERI_chem(munusila)
!        enddo
!      enddo
!    enddo
!  enddo



!  !$OMP PARALLEL SHARED (NONE)
!  !$OMP PRIVATE ()
!  !$OMP SHARED ()
!  !$OMP DO
  do nu = 1, nBas

    nunu0 = shiftr(nu * (nu - 1), 1)
    nunu = nunu0 + nu
    nunununu = shiftr(nunu * (nunu - 1), 1) + nunu
    K(nu,nu) = -P(nu,nu) * ERI_chem(nunununu)

    do la = 1, nu - 1
      nula = nunu0 + la
      nulanula = shiftr(nula * (nula - 1), 1) + nula
      K(nu,nu) = K(nu,nu) - P(la,la) * ERI_chem(nulanula)
    enddo

    do la = nu + 1, nBas
      lanu = shiftr(la * (la - 1), 1) + nu
      lanulanu = shiftr(lanu * (lanu - 1), 1) + lanu
      K(nu,nu) = K(nu,nu) - P(la,la) * ERI_chem(lanulanu)
    enddo

    do la = 1, nu
      nula = nunu0 + la
      nulanula0 = shiftr(nula * (nula - 1), 1)
      do si = 1, la - 1
        nulanusi = nulanula0 + nunu0 + si
        K(nu,nu) = K(nu,nu) - 2.d0 * P(si,la) * ERI_chem(nulanusi)
      enddo
    enddo

    do la = nu + 1, nBas
      lanu = shiftr(la * (la - 1), 1) + nu
      lanulanu0 = shiftr(lanu * (lanu - 1), 1)
      do si = 1, nu
        lanunusi = lanulanu0 + nunu0 + si
        K(nu,nu) = K(nu,nu) - 2.d0 * P(si,la) * ERI_chem(lanunusi)
      enddo
      do si = nu + 1, la - 1
        lanusinu = lanulanu0 + shiftr(si * (si - 1), 1) + nu
        K(nu,nu) = K(nu,nu) - 2.d0 * P(si,la) * ERI_chem(lanusinu)
      enddo
    enddo


    do mu = 1, nu - 1

      numu = nunu0 + mu
      mumu0 = shiftr(mu * (mu - 1), 1)
      mumu = mumu0 + mu
      numumumu = shiftr(numu * (numu - 1), 1) + mumu
      K(mu,nu) = - P(mu,mu) * ERI_chem(numumumu)

      do la = 1, mu - 1
        nula = nunu0 + la
        nulamula = shiftr(nula * (nula - 1), 1) + mumu0 + la
        K(mu,nu) = K(mu,nu) - P(la,la) * ERI_chem(nulamula)
      enddo
      do la = mu + 1, nu
        nula = nunu0 + la
        nulalamu = shiftr(nula * (nula - 1), 1) + shiftr(la * (la - 1), 1) + mu
        K(mu,nu) = K(mu,nu) - P(la,la) * ERI_chem(nulalamu)
      enddo
      do la = nu + 1, nBas
        lala0 = shiftr(la * (la - 1), 1)
        lanu = lala0 + nu
        lanulamu = shiftr(lanu * (lanu - 1), 1) + lala0 + mu
        K(mu,nu) = K(mu,nu) - P(la,la) * ERI_chem(lanulamu)
      enddo

      do la = 1, mu
        nula = nunu0 + la
        nulanula0 = shiftr(nula * (nula - 1), 1)
        do si = 1, la - 1
          nulamusi = nulanula0 + mumu0 + si
          K(mu,nu) = K(mu,nu) - P(si,la) * ERI_chem(nulamusi)
        enddo
      enddo
      do la = mu + 1, nu
        nula = nunu0 + la
        nulanula0 = shiftr(nula * (nula - 1), 1)
        do si = 1, mu
          nulamusi = nulanula0 + mumu0 + si
          K(mu,nu) = K(mu,nu) - P(si,la) * ERI_chem(nulamusi)
        enddo
        do si = mu + 1, la - 1
          nulasimu = nulanula0 + shiftr(si * (si - 1), 1) + mu
          K(mu,nu) = K(mu,nu) - P(si,la) * ERI_chem(nulasimu)
        enddo
      enddo
      do la = nu + 1, nBas
        lanu = shiftr(la * (la - 1), 1) + nu
        lanulanu0 = shiftr(lanu * (lanu - 1), 1)
        do si = 1, mu
          lanumusi = lanulanu0 + mumu0 + si
          K(mu,nu) = K(mu,nu) - P(si,la) * ERI_chem(lanumusi)
        enddo
        do si = mu + 1, la - 1
          lanusimu = lanulanu0 + shiftr(si * (si - 1), 1) + mu
          K(mu,nu) = K(mu,nu) - P(si,la) * ERI_chem(lanusimu)
        enddo
      enddo

      do la = 1, mu
        nula = nunu0 + la
        nulanula0 = shiftr(nula * (nula - 1), 1)
        do si = la + 1, mu
          nulamusi = nulanula0 + mumu0 + si
          K(mu,nu) = K(mu,nu) - P(si,la) * ERI_chem(nulamusi)
        enddo
        do si = mu + 1, nu - 1
          nulasimu = nulanula0 + shiftr(si * (si - 1), 1) + mu
          K(mu,nu) = K(mu,nu) - P(si,la) * ERI_chem(nulasimu)
        enddo
        do si = nu, nBas
          simu = shiftr(si * (si - 1), 1) + mu
          simunula = shiftr(simu * (simu - 1), 1) + nula
          K(mu,nu) = K(mu,nu) - P(si,la) * ERI_chem(simunula)
        enddo
      enddo
      do la = mu + 1, nu
        nula = nunu0 + la
        nulanula0 = shiftr(nula * (nula - 1), 1)
        do si = la + 1, nu
          nulasimu = nulanula0 + shiftr(si * (si - 1), 1) + mu
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
!  !$OMP END DO
!  !$OMP END PARALLEL


  do nu = 1, nBas
    do mu = nu+1, nBas
      K(mu,nu) = K(nu,mu)
    enddo
  enddo

  return
end subroutine 

! ---

