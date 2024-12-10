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
  integer                       :: nunu, lala, nula, lasi, numu
  integer*8                     :: nunununu, nunulala, nununula, nunulasi
  integer*8                     :: lalanunu, lasinunu, numulala, lalanumu
  integer*8                     :: numunula, numulasi, lasinumu, nununumu
  integer*8                     :: munu0, munu
  integer*8                     :: sila0, sila
  integer*8                     :: munulasi0, munulasi

  integer*8, external           :: Yoshimine_4ind


  do nu = 1, nBas

    munulasi = Yoshimine_4ind(nu, nu, nu, nu)
    K(nu,nu) = -P(nu,nu) * ERI_chem(munulasi)

    do la = 1, nu - 1
      munulasi = Yoshimine_4ind(nu, la, nu, la)
      K(nu,nu) = K(nu,nu) - P(la,la) * ERI_chem(munulasi)
    enddo

    do la = nu + 1, nBas
      munulasi = Yoshimine_4ind(nu, la, nu, la)
      K(nu,nu) = K(nu,nu) - P(la,la) * ERI_chem(munulasi)
    enddo

    do la = 1, nu
      do si = 1, la - 1
        munulasi = Yoshimine_4ind(nu, la, nu, si)
        K(nu,nu) = K(nu,nu) - 2.d0 * P(si,la) * ERI_chem(munulasi)
      enddo
    enddo

    do la = nu + 1, nBas
      do si = 1, nu
        munulasi = Yoshimine_4ind(nu, la, nu, si)
        K(nu,nu) = K(nu,nu) - 2.d0 * P(si,la) * ERI_chem(munulasi)
      enddo
    enddo

    do la = nu + 1, nBas
      do si = nu + 1, la - 1
        munulasi = Yoshimine_4ind(nu, la, nu, si)
        K(nu,nu) = K(nu,nu) - 2.d0 * P(si,la) * ERI_chem(munulasi)
      enddo
    enddo


    do mu = 1, nu - 1

      K(mu,nu) = 0.d0
      do la = 1, mu
        munulasi = Yoshimine_4ind(nu, la, mu, la)
        K(mu,nu) = K(mu,nu) - P(la,la) * ERI_chem(munulasi)
      enddo
      do la = mu + 1, nu
        munulasi = Yoshimine_4ind(nu, la, mu, la)
        K(mu,nu) = K(mu,nu) - P(la,la) * ERI_chem(munulasi)
      enddo
      do la = nu + 1, nBas
        munulasi = Yoshimine_4ind(nu, la, mu, la)
        K(mu,nu) = K(mu,nu) - P(la,la) * ERI_chem(munulasi)
      enddo

      do la = 1, mu
        do si = 1, la - 1
          munulasi = Yoshimine_4ind(nu, la, si, mu)
          K(mu,nu) = K(mu,nu) - P(si,la) * ERI_chem(munulasi)
        enddo
      enddo
      do la = mu+1, nu
        do si = 1, mu
          munulasi = Yoshimine_4ind(nu, la, si, mu)
          K(mu,nu) = K(mu,nu) - P(si,la) * ERI_chem(munulasi)
        enddo
        do si = mu + 1, la - 1
          munulasi = Yoshimine_4ind(nu, la, si, mu)
          K(mu,nu) = K(mu,nu) - P(si,la) * ERI_chem(munulasi)
        enddo
      enddo
      do la = nu + 1, nBas
        do si = 1, mu
          munulasi = Yoshimine_4ind(nu, la, si, mu)
          K(mu,nu) = K(mu,nu) - P(si,la) * ERI_chem(munulasi)
        enddo
        do si = mu + 1, la - 1
          munulasi = Yoshimine_4ind(nu, la, si, mu)
          K(mu,nu) = K(mu,nu) - P(si,la) * ERI_chem(munulasi)
        enddo
      enddo

      do la = 1, mu
        do si = la + 1, mu
          munulasi = Yoshimine_4ind(nu, la, si, mu)
          K(mu,nu) = K(mu,nu) - P(si,la) * ERI_chem(munulasi)
        enddo
        do si = mu + 1, nBas
          munulasi = Yoshimine_4ind(nu, la, si, mu)
          K(mu,nu) = K(mu,nu) - P(si,la) * ERI_chem(munulasi)
        enddo
      enddo
      do la = mu + 1, nu
        do si = la + 1, nBas
          munulasi = Yoshimine_4ind(nu, la, si, mu)
          K(mu,nu) = K(mu,nu) - P(si,la) * ERI_chem(munulasi)
        enddo
      enddo
      do la = nu + 1, nBas
        do si = la + 1, nBas
          munulasi = Yoshimine_4ind(nu, la, si, mu)
          K(mu,nu) = K(mu,nu) - P(si,la) * ERI_chem(munulasi)
        enddo
      enddo

      !do la = 1, nBas
      !  do si = la + 1, nBas
      !    munulasi = Yoshimine_4ind(nu, la, mu, si)
      !    K(mu,nu) = K(mu,nu) - P(si,la) * ERI_chem(munulasi)
      !  enddo
      !enddo
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

