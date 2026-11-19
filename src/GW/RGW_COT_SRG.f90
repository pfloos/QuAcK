subroutine RGW_COT_SRG(nAct,nAux,flow,F,B,D)

! Finite-flow dynamical SRG-GW through second order in the couplings.
! The auxiliary block D is kept at zeroth order in the resolvent.

  implicit none
  integer,intent(in) :: nAct,nAux
  double precision,intent(in) :: flow,D(nAux)
  double precision,intent(inout) :: F(nAct,nAct),B(nAct,nAux)

  integer :: p,q,alpha
  double precision :: x,gapSquared,kernel
  double precision :: eigenvalues(nAct),gaps(nAct),damping(nAct)
  double precision :: rotation(nAct,nAct),correction(nAct,nAct)
  double precision,allocatable :: couplings(:,:)

  if(flow == 0d0) return
  rotation=F
  call diagonalize_matrix(nAct,rotation,eigenvalues)
  allocate(couplings(nAct,nAux))
  couplings=matmul(transpose(rotation),B)
  correction=0d0
  do alpha=1,nAux
    gaps=eigenvalues-D(alpha)
    damping=exp(-flow*gaps**2)
    do q=1,nAct
      do p=1,q
        gapSquared=gaps(p)**2+gaps(q)**2
        x=flow*gapSquared
        ! Evaluate (1-exp(-x))/gapSquared continuously at zero gaps.
        if(x < 1d-4) then
          kernel=flow*(1d0-x*(0.5d0-x*(1d0/6d0-x/24d0)))
        else
          kernel=(1d0-exp(-x))/gapSquared
        end if
        correction(p,q)=correction(p,q)+couplings(p,alpha)*couplings(q,alpha)*(gaps(p)+gaps(q))*kernel
      end do
    end do
    couplings(:,alpha)=couplings(:,alpha)*damping
  end do
  do q=1,nAct
    correction(q,1:q-1)=correction(1:q-1,q)
  end do
  correction=matmul(rotation,matmul(correction,transpose(rotation)))
  F=F+0.5d0*(correction+transpose(correction))
  B=matmul(rotation,couplings)

end subroutine RGW_COT_SRG
