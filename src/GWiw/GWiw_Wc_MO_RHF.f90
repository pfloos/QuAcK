subroutine GWiw_Wc_MO_RHF(nOrb,nO,e,omega,vMAT,Wres)

! Correlation part of the RPA screened interaction at the complex frequency omega, restricted case, MO basis:
!
!   W_c(omega) = v [ 1 - X0(omega) v ]^-1 X0(omega) v
!
! (same construction as Sigmac_MO_RHF_GW_w, with the spin factor 2 inside X0), and the matrix elements needed in
! the residue terms of the contour deformation,
!
!   Wres(p,q,r) = < phi_p phi_r | W_c(omega) | phi_r phi_q >
!
! The energies e must be measured from any origin (only differences enter). W_c is even in omega.

  implicit none

! Input variables

  integer,intent(in)            :: nOrb,nO
  double precision,intent(in)   :: e(nOrb)
  complex(kind=8),intent(in)    :: omega
  double precision,intent(in)   :: vMAT(nOrb*nOrb,nOrb*nOrb)

! Output variables

  complex(kind=8),intent(out)   :: Wres(nOrb,nOrb,nOrb)

! Local variables

  integer                       :: n2,i,a,p,q,r,x,info
  integer,allocatable           :: ipiv(:)
  double precision              :: Delta
  complex(kind=8)               :: chi
  complex(kind=8),parameter     :: cone = (1d0,0d0), czero = (0d0,0d0)
  complex(kind=8),allocatable   :: M(:,:),Y(:,:),T(:,:),V(:,:),chi0(:)

  n2 = nOrb*nOrb

  allocate(M(n2,n2),Y(n2,n2),T(n2,n2),V(n2,n2),chi0(n2),ipiv(n2))

  V = cmplx(vMAT,0d0,kind=8)

! Independent-particle response (diagonal in the pair index of the (occupied,virtual) pairs)

  chi0 = czero
  do i=1,nO
    do a=nO+1,nOrb
      Delta = e(a) - e(i)
      chi0(1+(a-1)+(i-1)*nOrb) = 2d0*( 1d0/(omega-Delta) - 1d0/(omega+Delta) )
    end do
  end do

! M = 1 - X0 v ;  Y = M^-1 X0

  do x=1,n2
    M(x,:) = -chi0(x)*V(x,:)
    M(x,x) = M(x,x) + cone
  end do

  Y = czero
  do x=1,n2
    Y(x,x) = chi0(x)
  end do

  call zgesv(n2,n2,M,n2,ipiv,Y,n2,info)
  if(info /= 0) then
    write(*,*) ' Error in the RPA Dyson equation at complex frequency, info = ',info
    error stop
  end if

! W_c = v Y v

  call zgemm('N','N',n2,n2,n2,cone,Y,n2,V,n2,czero,T,n2)
  call zgemm('N','N',n2,n2,n2,cone,V,n2,T,n2,czero,Y,n2)

  do r=1,nOrb
    do q=1,nOrb
      do p=1,nOrb
        Wres(p,q,r) = Y(1+(r-1)+(p-1)*nOrb,1+(q-1)+(r-1)*nOrb)
      end do
    end do
  end do

  deallocate(M,Y,T,V,chi0,ipiv)

end subroutine
