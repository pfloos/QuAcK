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
  integer                       :: nunu0, lala0
  integer*8                     :: nunununu, nunulala, nununula, nunulasi
  integer*8                     :: lalanunu, lasinunu, numulala, lalanumu
  integer*8                     :: numunula, numulasi, lasinumu, nununumu
  integer*8                     :: nunununu0, numunumu0

!  integer*8                     :: munusila
!  integer*8, external           :: Yoshimine_4ind
!
!  do nu = 1, nBas
!    do mu = 1, nu
!      H(mu,nu) = 0.d0
!      do la = 1, nBas
!        do si = 1, nBas
!          munusila = Yoshimine_4ind(int(mu, kind=8), &
!                                    int(nu, kind=8), &
!                                    int(si, kind=8), &
!                                    int(la, kind=8))
!          H(mu,nu) = H(mu,nu) + P(si,la) * ERI_chem(munusila)
!        enddo
!      enddo
!    enddo
!  enddo

  !$OMP PARALLEL DEFAULT(NONE)                                      &
  !$OMP PRIVATE (nu, la, si, mu,                                    &
  !$OMP          nunu0, nunu, nula, lala0, lala, lasi, numu,        &
  !$OMP          nunununu0, nunununu, nununula, numulala, numunula, &
  !$OMP          nunulala, lalanunu, lalanumu, nunulasi, lasinunu,  &
  !$OMP          numunumu0, nununumu, numulasi, lasinumu)           &
  !$OMP SHARED (nBas, H, P, ERI_chem)
  !$OMP DO
  do nu = 1, nBas

    nunu0 = shiftr(nu * (nu - 1), 1)
    nunu = nunu0 + nu
    nunununu0 = shiftr(nunu * (nunu - 1), 1)

    nunununu = nunununu0 + nunu
    H(nu,nu) = P(nu,nu) * ERI_chem(nunununu)

    do la = 1, nu - 1

      lala0 = shiftr(la * (la - 1), 1)

      lala = lala0 + la
      nunulala = nunununu0 + lala
      H(nu,nu) = H(nu,nu) + P(la,la) * ERI_chem(nunulala)

      nula = nunu0 + la
      nununula = nunununu0 + nula
      H(nu,nu) = H(nu,nu) + 2.d0 * P(la,nu) * ERI_chem(nununula)

      do si = 1, la - 1
        lasi = lala0 + si
        nunulasi = nunununu0 + lasi
        H(nu,nu) = H(nu,nu) + 2.d0 * P(si,la) * ERI_chem(nunulasi)
      enddo
    enddo

    do la = nu + 1, nBas

      lala0 = shiftr(la * (la - 1), 1)

      lala = lala0 + la
      lalanunu = shiftr(lala * (lala - 1), 1) + nunu
      H(nu,nu) = H(nu,nu) + P(la,la) * ERI_chem(lalanunu)

      do si = 1, la - 1
        lasi = lala0 + si
        lasinunu = shiftr(lasi * (lasi - 1), 1) + nunu
        H(nu,nu) = H(nu,nu) + 2.d0 * P(si,la) * ERI_chem(lasinunu)
      enddo
    enddo

    do mu = 1, nu - 1

      numu = nunu0 + mu

      numunumu0 = shiftr(numu * (numu - 1), 1)

      nununumu = nunununu0 + numu
      H(mu,nu) = p(nu,nu) * ERI_chem(nununumu)

      do la = 1, nu - 1
        lala = shiftr(la * (la - 1), 1) + la
        numulala = numunumu0 + lala
        H(mu,nu) = H(mu,nu) + p(la,la) * ERI_chem(numulala)
      enddo

      do la = nu + 1, nBas
        lala = shiftr(la * (la - 1), 1) + la
        lalanumu = shiftr(lala * (lala - 1), 1) + numu
        H(mu,nu) = H(mu,nu) + p(la,la) * ERI_chem(lalanumu)
      enddo

      do la = 1, mu
        nula = nunu0 + la
        numunula = numunumu0 + nula
        H(mu,nu) = H(mu,nu) + 2.d0 * P(la,nu) * ERI_chem(numunula)
      enddo

      do la = mu + 1, nu - 1
        nula = nunu0 + la
        numunula = shiftr(nula * (nula - 1), 1) + numu
        H(mu,nu) = H(mu,nu) + 2.d0 * P(la,nu) * ERI_chem(numunula)
      enddo

      do la = 2, nu - 1
        lala0 = shiftr(la * (la - 1), 1)
        do si = 1, la - 1
          lasi = lala0 + si
          numulasi = numunumu0 + lasi
          H(mu,nu) = H(mu,nu) + 2.d0 * P(si,la) * ERI_chem(numulasi)
        enddo
      enddo

      do la = nu + 1, nBas
        lala0 = shiftr(la * (la - 1), 1)
        do si = 1, la - 1
          lasi = lala0 + si
          lasinumu = shiftr(lasi * (lasi - 1), 1) + numu
          H(mu,nu) = H(mu,nu) + 2.d0 * P(si,la) * ERI_chem(lasinumu)
        enddo
      enddo

    enddo ! mu
  enddo ! nu
  !$OMP END DO
  !$OMP END PARALLEL

  do nu = 1, nBas
    do mu = nu+1, nBas
      H(mu,nu) = H(nu,mu)
    enddo
  enddo

  return
end subroutine 

! ---


