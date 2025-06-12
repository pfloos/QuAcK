subroutine complex_cRGF2_SRG_self_energy_diag(flow,eta,nBas,nC,nO,nV,nR,Re_e,Im_e,ERI,Re_SigC,Im_SigC,Re_Z,Im_Z)

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
  double precision,intent(in)   :: flow
  double precision,intent(in)   :: Re_e(nBas)
  double precision,intent(in)   :: Im_e(nBas)
  complex*16,intent(in)         :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: i,j,a,b
  integer                       :: p
  double precision              :: eps
  double precision              :: s
  double precision              :: eta_tilde
  complex*16                    :: num
  double precision,allocatable  :: Re_DS(:)
  double precision,allocatable  :: Im_DS(:)
  complex*16                    :: z_dummy

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
  s = flow 

! Compute GF2 self-energy

  do p=nC+1,nBas-nR
    do i=nC+1,nO
      do j=nC+1,nO
        do a=nO+1,nBas-nR

          eps = Re_e(p) + Re_e(a) - Re_e(i) - Re_e(j)
          eta_tilde = eta - Im_e(p) + Im_e(i)  - (Im_e(a) - Im_e(j))
          num = (2d0*ERI(p,a,i,j) - ERI(p,a,j,i))*ERI(p,a,i,j) &
                  *(1d0 - exp(-2d0*s*(eps**2 + eta_tilde**2)))

          z_dummy = num*cmplx(eps/(eps**2 + eta_tilde**2),eta_tilde/(eps**2 + eta_tilde**2),kind=8)
          Re_SigC(p) = Re_SigC(p) + real(z_dummy)
          Im_SigC(p) = Im_SigC(p) + aimag(z_dummy)
          z_dummy = num*cmplx(-(eps**2 - eta_tilde**2)/(eps**2 + eta_tilde**2)**2,&
                  -2*eta_tilde*eps/(eps**2 + eta_tilde**2)**2,kind=8)
          Re_DS(p)    = Re_DS(p)   + real(z_dummy)
          Im_DS(p)   = Im_DS(p)    + aimag(z_dummy) 

        end do
      end do
    end do
  end do

  do p=nC+1,nBas-nR
    do i=nC+1,nO
      do a=nO+1,nBas-nR
        do b=nO+1,nBas-nR

          eps = Re_e(p) + Re_e(i) - Re_e(a) - Re_e(b)
          eta_tilde = eta + Im_e(p)  - Im_e(a)  - Im_e(b) + Im_e(i)
          num = (2d0*ERI(p,i,a,b) - ERI(p,i,b,a))*ERI(p,i,a,b)&
                  *(1d0 - exp(-2d0*s*(eps**2 + eta_tilde**2)))

          z_dummy = num*cmplx(eps/(eps**2 + eta_tilde**2),-eta_tilde/(eps**2 + eta_tilde**2),kind=8)
          Re_SigC(p) = Re_SigC(p) + real(z_dummy)
          Im_SigC(p) = Im_SigC(p) + aimag(z_dummy)
          z_dummy = num*cmplx(-(eps**2 - eta_tilde**2)/(eps**2 + eta_tilde**2)**2,&
                  2*eta_tilde*eps/(eps**2 + eta_tilde**2)**2,kind=8)
          Re_DS(p)    = Re_DS(p)   + real(z_dummy)
          Im_DS(p)   = Im_DS(p)    + aimag(z_dummy) 

        end do
      end do
    end do
  end do
  Re_Z(:) = (1d0-Re_DS(:))/((1d0 - Re_DS(:))**2 + Im_DS(:)**2)
  Im_Z(:) = Im_DS(:)/((1d0 - Re_DS(:))**2 + Im_DS(:)**2)
  deallocate(Re_DS,Im_DS)
end subroutine 
