subroutine RGW_COT_SRG(nAct,nAux,flow,eRef,B,D)

! Freeze the dynamical SRG-GW diagonal residue factors at the reference energies.
! Applying their square roots to B gives a separable, Hermitian supermatrix.

  implicit none
  integer,intent(in) :: nAct,nAux
  double precision,intent(in) :: flow,eRef(nAct),D(nAux)
  double precision,intent(inout) :: B(nAct,nAux)

  integer :: p,alpha
  double precision :: x,factor

  if(flow == 0d0) then
    B=0d0
    return
  end if

  do alpha=1,nAux
    do p=1,nAct
      x=2d0*flow*(eRef(p)-D(alpha))**2
      ! Evaluate 1-exp(-x) without cancellation near a zero gap.
      if(x < 1d-4) then
        factor=x*(1d0-x*(0.5d0-x*(1d0/6d0-x/24d0)))
      else
        factor=1d0-exp(-x)
      end if
      B(p,alpha)=B(p,alpha)*sqrt(factor)
    end do
  end do

end subroutine RGW_COT_SRG
