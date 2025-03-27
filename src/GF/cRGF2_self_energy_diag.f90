subroutine cRGF2_self_energy_diag(eta,nBas,nC,nO,nV,nR,e,ERI,Re_SigC,Im_SigC,Re_Z,Im_Z,e_cap)

! Compute diagonal part of the GF2 self-energy and its renormalization factor

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  double precision,intent(in)   :: e(nBas)
  double precision,intent(in)   :: e_cap(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: i,j,a,b
  integer                       :: p
  double precision              :: eps
  double precision              :: eta_tilde
  double precision              :: num
  double precision,allocatable  :: Re_DS(:)
  double precision,allocatable  :: Im_DS(:)

! Output variables

  double precision,intent(out)  :: Re_SigC(nBas)
  double precision,intent(out)  :: Im_SigC(nBas)
  double precision,intent(out)  :: Re_Z(nBas)
  double precision,intent(out)  :: Im_Z(nBas)

! Initialize 
  allocate(Re_DS(nBas),Im_DS(nBas))
  Re_SigC(:) = 0d0
  Im_SigC(:) = 0d0
  Re_DS(:)    = 0d0
  Im_DS(:)    = 0d0
 

! Compute GF2 self-energy

  do p=nC+1,nBas-nR
    do i=nC+1,nO
      do j=nC+1,nO
        do a=nO+1,nBas-nR

          eps = e(p) + e(a) - e(i) - e(j)
          eta_tilde = eta - e_cap(p) + e_cap(i) + e_cap(a) - e_cap(j)
          num = (2d0*ERI(p,a,i,j) - ERI(p,a,j,i))*ERI(p,a,i,j)

          Re_SigC(p) = Re_SigC(p) + num*eps/(eps**2 + eta_tilde**2)
          Im_SigC(p) = Im_SigC(p) + num*eta_tilde/(eps**2 + eta_tilde**2)
          Re_DS(p)    = Re_DS(p)   - num*(eps**2 - eta_tilde**2)/(eps**2 + eta_tilde**2)**2
          Im_DS(p)   = Im_DS(p)    + 2*num*eta_tilde*eps/(eps**2 + eta_tilde**2)**2

        end do
      end do
    end do
  end do

  do p=nC+1,nBas-nR
    do i=nC+1,nO
      do a=nO+1,nBas-nR
        do b=nO+1,nBas-nR

          eps = e(p) + e(i) - e(a) - e(b)
          eta_tilde = eta + e_cap(p)  - e_cap(a)  - e_cap(b) + e_cap(i)
          num = (2d0*ERI(p,i,a,b) - ERI(p,i,b,a))*ERI(p,i,a,b)

          Re_SigC(p) = Re_SigC(p) + num*eps/(eps**2 + eta_tilde**2)
          Im_SigC(p) = Im_SigC(p) - num*eta_tilde/(eps**2 + eta_tilde**2)
          Re_DS(p)    = Re_DS(p)   - num*(eps**2 - eta_tilde**2)/(eps**2 + eta_tilde**2)**2
          Im_DS(p)   = Im_DS(p)    + 2*num*eta_tilde*eps/(eps**2 + eta_tilde**2)**2

        end do
      end do
    end do
  end do
  Re_Z(:) = (1d0-Re_DS(:))/((1d0 - Re_DS(:))**2 + Im_DS(:)**2)
  Im_Z(:) = -Im_DS(:)/((1d0 - Re_DS(:))**2 + Im_DS(:)**2)
  deallocate(Re_DS,Im_DS)
end subroutine 
