subroutine scCOT(TDA_COT,useDavidson,davidsonThreshold,maxSCF,threshold,maxDIIS, &
                 nBas,nOrb,nC,nO,nV,nR,nS,ENuc,T,V,Hc,ERI_AO,cHF,eHF,info)

! Iterate the Hermitian COT Hamiltonian in the fixed RHF orbital basis.
! The frozen occupied and removed virtual orbitals are held fixed.

  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

  logical,intent(in) :: TDA_COT,useDavidson
  integer,intent(in) :: maxSCF,maxDIIS,nBas,nOrb,nC,nO,nV,nR,nS
  double precision,intent(in) :: davidsonThreshold,threshold,ENuc
  double precision,intent(in) :: T(nBas,nBas),V(nBas,nBas),Hc(nBas,nBas)
  double precision,intent(in) :: ERI_AO(nBas,nBas,nBas,nBas),cHF(nBas,nOrb),eHF(nOrb)
  integer,intent(out) :: info

  integer :: first,last,nAct,iteration,i,nDIIS,historySize
  double precision :: convergence,rcond,totalEnergy
  double precision,allocatable :: c(:,:),rotation(:,:),eRef(:),eCOT(:)
  double precision,allocatable :: P(:,:),J(:,:),K(:,:),F(:,:),fock(:,:),T_MO(:,:),V_MO(:,:)
  double precision,allocatable :: ERI(:,:,:,:),hCOT(:,:),hHerm(:,:),hInput(:,:),hOutput(:,:),error(:,:)
  double precision,allocatable :: errorHistory(:,:),hHistory(:,:),errorVector(:),candidate(:)

  first=nC+1
  last=nOrb-nR
  nAct=last-first+1
  if(nAct < 1 .or. nO < nC .or. nO > last) then
    info=1
    return
  end if
  allocate(c(nBas,nOrb),rotation(nAct,nAct),eRef(nOrb),eCOT(nOrb))
  allocate(P(nBas,nBas),J(nBas,nBas),K(nBas,nBas),F(nBas,nBas))
  allocate(fock(nOrb,nOrb),T_MO(nOrb,nOrb),V_MO(nOrb,nOrb),ERI(nOrb,nOrb,nOrb,nOrb))
  allocate(hCOT(nOrb,nOrb),hHerm(nOrb,nOrb),hInput(nAct,nAct),hOutput(nAct,nAct),error(nAct,nAct))
  historySize=min(maxDIIS,maxSCF)
  if(historySize > 1) then
    allocate(errorHistory(nAct*nAct,historySize),hHistory(nAct*nAct,historySize))
    allocate(errorVector(nAct*nAct),candidate(nAct*nAct))
  end if
  nDIIS=0

  c=cHF
  eRef=eHF
  fock=0d0
  rotation=0d0
  do i=1,nOrb
    fock(i,i)=eHF(i)
  end do
  do i=1,nAct
    rotation(i,i)=1d0
  end do
  hInput=fock(first:last,first:last)
  T_MO=0d0
  V_MO=0d0

  write(*,'(/,A)') 'scCOT self-consistency'
  write(*,'(A)') ' Iteration   Hamiltonian residual (Eh)'
  do iteration=1,maxSCF
    if(iteration > 1) then
      P=2d0*matmul(c(:,1:nO),transpose(c(:,1:nO)))
      call Hartree_matrix_AO_basis(nBas,P,ERI_AO,J)
      call exchange_matrix_AO_basis(nBas,P,ERI_AO,K)
      F=Hc+J+0.5d0*K
      call AOtoMO(nBas,nOrb,c,F,fock)
      fock=0.5d0*(fock+transpose(fock))
    end if
    call AOtoMO_ERI_RHF(nBas,nOrb,c,ERI_AO,ERI)
    call osCOT(TDA_COT,useDavidson,davidsonThreshold,nOrb,nC,nO,nV,nR,nS, &
               ENuc,T_MO,V_MO,fock,eRef,ERI,hCOT,hHerm,eCOT,totalEnergy,info,'')
    if(info /= 0) return

    hOutput=matmul(rotation,matmul(hHerm(first:last,first:last),transpose(rotation)))
    hOutput=0.5d0*(hOutput+transpose(hOutput))
    error=hOutput-hInput
    convergence=maxval(abs(error))
    write(*,'(I10,3X,ES24.12)') iteration,convergence
    if(.not.all(ieee_is_finite(error))) exit
    if(convergence <= threshold) then
      write(*,'(A,I0,A)') 'scCOT converged in ',iteration,' iteration(s).'
      call AOtoMO(nBas,nOrb,c,T,T_MO)
      call AOtoMO(nBas,nOrb,c,V,V_MO)
      call osCOT(TDA_COT,useDavidson,davidsonThreshold,nOrb,nC,nO,nV,nR,nS, &
                 ENuc,T_MO,V_MO,fock,eRef,ERI,hCOT,hHerm,eCOT,totalEnergy,info,'scCOT input')
      return
    end if
    if(iteration == maxSCF) exit

    hInput=hOutput
    if(historySize > 1) then
      nDIIS=min(nDIIS+1,historySize)
      errorVector=reshape(error,[nAct*nAct])
      candidate=reshape(hOutput,[nAct*nAct])
      call DIIS_extrapolation(rcond,nAct*nAct,nAct*nAct,nDIIS,errorHistory,hHistory,errorVector,candidate)
      if(rcond > 1d-13 .and. all(ieee_is_finite(candidate))) then
        hInput=reshape(candidate,[nAct,nAct])
        hInput=0.5d0*(hInput+transpose(hInput))
      else
        ! Discard a linearly dependent DIIS history and take the plain update.
        nDIIS=0
      end if
    end if
    rotation=hInput
    call diagonalize_matrix(nAct,rotation,eRef(first:last))
    c(:,first:last)=matmul(cHF(:,first:last),rotation)
  end do

  info=9
  write(*,'(A,I0,A)') 'scCOT did not converge within ',maxSCF,' iteration(s).'

end subroutine scCOT
