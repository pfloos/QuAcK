subroutine complex_cRGF2_reg_self_energy(flow,eta,nBas,nC,nO,nV,nR,e,ERI,SigC,Z)

! Compute diagonal part of the GF2 self-energy and its renormalization factor

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: eta
  double precision,intent(in)   :: flow
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  complex*16,intent(in)         :: e(nBas)
  complex*16,intent(in)         :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: i,j,a,b
  integer                       :: p,q
  double precision              :: eps_p
  double precision              :: eps_q
  double precision              :: s
  double precision              :: eta_tilde_p
  double precision              :: eta_tilde_q
  complex*16                    :: num
  double precision,allocatable  :: Re_DS(:)
  double precision,allocatable  :: Im_DS(:)
  complex*16                    :: z_dummy
  double precision,allocatable  :: Re_SigC(:,:)
  double precision,allocatable  :: Im_SigC(:,:)
  double precision,allocatable  :: Re_Z(:)
  double precision,allocatable  :: Im_Z(:)

! Output variables

  complex*16,intent(out)        :: SigC(nBas,nBas)
  complex*16,intent(out)        :: Z(nBas)

! Initialize 
  allocate(Re_DS(nBas),Im_DS(nBas),Re_SigC(nBas,nBas),Im_SigC(nBas,nBas),&
          Re_Z(nBas),Im_Z(nBas))
  Re_SigC(:,:) = 0d0
  Im_SigC(:,:) = 0d0
  Re_DS(:)    = 0d0
  Im_DS(:)    = 0d0
  s = flow 

! Compute GF2 self-energy

 !$OMP PARALLEL &
 !$OMP SHARED(Re_DS,Im_DS,Re_SigC,Im_SigC,ERI,eta,nC,nO,nBas,nR,e,s) &
 !$OMP PRIVATE(p,q,i,j,a,eps_p,eps_q,num,eta_tilde_p,eta_tilde_q,z_dummy) &
 !$OMP DEFAULT(NONE)
 !$OMP DO
  do p=nC+1,nBas-nR
    do q=nC+1,nBas-nR
      do i=nC+1,nO
        do j=nC+1,nO
          do a=nO+1,nBas-nR

            eps_p = real(e(p)) + real(e(a)) - real(e(i)) - real(e(j))
            eps_q = real(e(q)) + real(e(a)) - real(e(i)) - real(e(j))
            eta_tilde_p = eta - aimag(e(p)) + aimag(e(i))  - (aimag(e(a)) - aimag(e(j)))
            eta_tilde_q = eta - aimag(e(q)) + aimag(e(i))  - (aimag(e(a)) - aimag(e(j)))
            num = (2d0*ERI(p,a,i,j) - ERI(p,a,j,i))*ERI(q,a,i,j)&
                   *(1d0 - exp(-s*(eps_p**2+eta_tilde_p**2 + eps_q**2 + eta_tilde_q**2)))

            z_dummy = num*cmplx((eps_p + eps_q)/(eps_p**2 + eta_tilde_p**2 +eps_q**2 + eta_tilde_q**2 ),&
                    (eta_tilde_p + eta_tilde_q)/(eps_p**2 + eta_tilde_p**2 +eps_q**2 + eta_tilde_q**2 ),kind=8)
            Re_SigC(p,q) = Re_SigC(p,q) + real(z_dummy)
            Im_SigC(p,q) = Im_SigC(p,q) + aimag(z_dummy)
            if(p==q) then
              z_dummy = num*cmplx(-(eps_p**2 - eta_tilde_p**2)/(eps_p**2 + eta_tilde_p**2)**2,&
                      -2*eta_tilde_p*eps_p/(eps_p**2 + eta_tilde_p**2)**2,kind=8)
              Re_DS(p)    = Re_DS(p)   + real(z_dummy)
              Im_DS(p)   = Im_DS(p)    + aimag(z_dummy)
            end if

          end do
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  

  !$OMP PARALLEL &
  !$OMP SHARED(Re_DS,Im_DS,Re_SigC,Im_SigC,ERI,eta,nC,nO,nBas,nR,e,s) &
  !$OMP PRIVATE(p,q,i,a,b,eps_p,eps_q,num,eta_tilde_p,eta_tilde_q,z_dummy) &
  !$OMP DEFAULT(NONE)
  !$OMP DO  
  do p=nC+1,nBas-nR
    do q=nC+1,nBas-nR
      do i=nC+1,nO
        do a=nO+1,nBas-nR
          do b=nO+1,nBas-nR

            eps_p = real(e(p)) + real(e(i)) - real(e(a)) - real(e(b))
            eps_q = real(e(q)) + real(e(i)) - real(e(a)) - real(e(b))
            eta_tilde_p = eta + aimag(e(p))  - aimag(e(a))  - aimag(e(b)) + aimag(e(i))
            eta_tilde_q = eta + aimag(e(q))  - aimag(e(a))  - aimag(e(b)) + aimag(e(i))
            num = (2d0*ERI(p,i,a,b) - ERI(p,i,b,a))*ERI(q,i,a,b)&
                    *(1d0 - exp(-s*(eps_p**2+eta_tilde_p**2 + eps_q**2 + eta_tilde_q**2)))


            z_dummy = num*cmplx((eps_p + eps_q)/(eps_p**2+eta_tilde_p**2 + eps_q**2 + eta_tilde_q**2),&
                    -(eta_tilde_p + eta_tilde_q)/(eps_p**2+eta_tilde_p**2 + eps_q**2 + eta_tilde_q**2),kind=8)
            Re_SigC(p,q) = Re_SigC(p,q) + real(z_dummy)
            Im_SigC(p,q) = Im_SigC(p,q) + aimag(z_dummy)
            if(p==q) then
              z_dummy = num*cmplx(-(eps_p**2 - eta_tilde_p**2)/(eps_p**2 + eta_tilde_p**2)**2,&
                      2*eta_tilde_p*eps_p/(eps_p**2 + eta_tilde_p**2)**2,kind=8)
              Re_DS(p)    = Re_DS(p)   + real(z_dummy)
              Im_DS(p)   = Im_DS(p)    + aimag(z_dummy) 
            end if

          end do
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  
  Re_Z(:) = (1d0-Re_DS(:))/((1d0 - Re_DS(:))**2 + Im_DS(:)**2)
  Im_Z(:) = Im_DS(:)/((1d0 - Re_DS(:))**2 + Im_DS(:)**2)
  
  Z = cmplx(Re_Z,Im_Z,kind=8)
  SigC = cmplx(Re_SigC,Im_SigC,kind=8)
  

  deallocate(Re_DS,Im_DS,Re_Z,Im_Z,Re_SigC,Im_SigC)
end subroutine 
