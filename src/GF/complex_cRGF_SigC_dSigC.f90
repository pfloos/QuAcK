subroutine complex_cRGF_SigC_dSigC(p,eta,nBas,nC,nO,nV,nR,Re_w,Im_w,Re_e,Im_e,ERI,Re_SigC,Im_SigC,Re_DS,Im_DS)


  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: p
  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  double precision,intent(in)   :: Re_e(nBas)
  double precision,intent(in)   :: Im_e(nBas)
  double precision,intent(in)   :: Re_w
  double precision,intent(in)   :: Im_w
  complex*16,intent(in)         :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: i,j,a,b
  double precision              :: eps
  double precision              :: eta_tilde
  complex*16                    :: num
  complex*16                    :: z_dummy

! Output variables

  double precision,intent(out)  :: Re_SigC
  double precision,intent(out)  :: Im_SigC
  double precision,intent(out)  :: Re_DS
  double precision,intent(out)  :: Im_DS

! Initialize 
  Re_SigC = 0d0
  Im_SigC = 0d0
  Re_DS    = 0d0
  Im_DS    = 0d0
 

! Compute GF2 self-energy
write(*,*) "DEbugging change back"

    do i=nC+1,nO
      do j=nC+1,nO
        do a=nO+1,nBas-nR

          eps = Re_w + Re_e(a) - Re_e(i) - Re_e(j)
          eta_tilde = eta - Im_w + Im_e(i) + Im_e(a) - Im_e(j)
          num = (2d0*ERI(p,a,i,j) - 0*ERI(p,a,j,i))*ERI(p,a,i,j)
          z_dummy = num*cmplx(eps/(eps**2 + eta_tilde**2),eta_tilde/(eps**2 + eta_tilde**2),kind=8)
          Re_SigC = Re_SigC + real(z_dummy)
          Im_SigC = Im_SigC + aimag(z_dummy)
          z_dummy = num*cmplx(-(eps**2 - eta_tilde**2)/(eps**2 + eta_tilde**2)**2,&
                  -2*eta_tilde*eps/(eps**2 + eta_tilde**2)**2,kind=8)
          Re_DS    = Re_DS   + real(z_dummy)
          Im_DS   = Im_DS    + aimag(z_dummy) 

        end do
      end do
    end do

    do i=nC+1,nO
      do a=nO+1,nBas-nR
        do b=nO+1,nBas-nR

          eps = Re_w + Re_e(i) - Re_e(a) - Re_e(b)
          eta_tilde = eta + Im_w  - Im_e(a)  - Im_e(b) + Im_e(i)
          num = (2d0*ERI(p,i,a,b) - 0*ERI(p,i,b,a))*ERI(p,i,a,b)

          z_dummy = num*cmplx(eps/(eps**2 + eta_tilde**2),-eta_tilde/(eps**2 + eta_tilde**2),kind=8)
          Re_SigC = Re_SigC + real(z_dummy)
          Im_SigC = Im_SigC + aimag(z_dummy)
          z_dummy = num*cmplx(-(eps**2 - eta_tilde**2)/(eps**2 + eta_tilde**2)**2,&
                  2*eta_tilde*eps/(eps**2 + eta_tilde**2)**2,kind=8)
          Re_DS    = Re_DS   + real(z_dummy)
          Im_DS   = Im_DS    + aimag(z_dummy) 

        end do
      end do
    end do
end subroutine 
