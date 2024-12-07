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
  integer*8                     :: munu0, munu
  integer*8                     :: sila0, sila
  integer*8                     :: munulasi0, munulasi

  integer*8, external           :: Yoshimine_ind

  do nu = 1, nBas
    do mu = 1, nBas
      H(mu,nu) = 0.d0
      do si = 1, nBas
        do la = 1, nBas
          munulasi = Yoshimine_ind(mu, nu, la, si)
          H(mu,nu) = H(mu,nu) + P(la,si) * ERI_chem(munulasi)
        enddo
      enddo
    enddo
  enddo


!  do nu = 1, nBas
!    munu0 = (nu * (nu + 1)) / 2
!
!    do mu = 1, nu
!      munu = munu0 + mu
!      munulasi0 = (munu * (munu + 1)) / 2
!
!      H(mu,nu) = 0.d0
!
!      do si = 1, nu
!        sila0 = (si * (si + 1)) / 2
!
!        do la = 1, si
!          sila = sila0 + la
!
!          if(nu == si .and. mu < la) cycle
!
!          munulasi = munulasi0 + sila
!
!          H(mu,nu) = H(mu,nu) + 4.d0 * P(la,si) * ERI_chem(munulasi)
!        enddo
!      enddo
!    enddo
!  enddo
!
!
!  do nu = 1, nBas
!    do mu = nu+1, nBas
!      H(mu,nu) = H(nu,mu)
!    enddo
!  enddo

  return
end subroutine 

! ---

integer*8 function Yoshimine_ind(a, b, c, d)

  implicit none

  integer, intent(in) :: a, b, c, d

  integer*8           :: ab, cd, abcd

  if(a > b) then
    ab = (a * (a - 1)) / 2 + b
  else
    ab = (b * (b - 1)) / 2 + a
  endif

  if(c > d) then
    cd = (c * (c - 1)) / 2 + d
  else
    cd = (d * (d - 1)) / 2 + c
  endif

  if(ab > cd) then
    abcd = (ab * (ab - 1)) / 2 + cd
  else
    abcd = (cd * (cd - 1)) / 2 + ab
  endif

  Yoshimine_ind = abcd

  return
end 

! ---

