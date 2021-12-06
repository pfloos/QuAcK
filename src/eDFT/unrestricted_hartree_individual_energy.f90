subroutine unrestricted_hartree_individual_energy(nBas,nEns,Pw,P,ERI,doNcentered,kappa,LZH,EH)

! Compute the hartree contribution to the individual energies

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: Pw(nBas,nBas,nspin)
  double precision,intent(in)   :: P(nBas,nBas,nspin,nEns)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: kappa(nEns)
  logical,intent(in)            :: doNcentered


! Local variables

  double precision,allocatable  :: J(:,:,:)
  double precision,external     :: trace_matrix

  integer                       :: iEns
  integer                       :: ispin

! Output variables

  double precision,intent(out)  :: LZH(nsp,nEns)
  double precision,intent(out)  :: EH(nsp,nEns)

! Compute HF exchange matrix

  allocate(J(nBas,nBas,nspin))

  LZH(:,:) = 0.d0
  EH(:,:)  = 0.d0

  do ispin=1,nspin
    call unrestricted_hartree_potential(nBas,Pw(:,:,ispin),ERI,J(:,:,ispin))
  end do

  do iEns=1,nEns

!    if(doNcentered) then
   
!      LZH(1,iEns) = - 0.5d0*kappa(iEns)*kappa(iEns)*trace_matrix(nBas,matmul(Pw(:,:,1),J(:,:,1)))
!      LZH(2,iEns) = - 0.5d0*kappa(iEns)*kappa(iEns)*trace_matrix(nBas,matmul(Pw(:,:,1),J(:,:,2))) &
!                    - 0.5d0*kappa(iEns)*kappa(iEns)*trace_matrix(nBas,matmul(Pw(:,:,2),J(:,:,1)))
!      LZH(3,iEns) = - 0.5d0*kappa(iEns)*trace_matrix(nBas,matmul(Pw(:,:,2),J(:,:,2)))

!      EH(1,iEns) = kappa(iEns)*trace_matrix(nBas,matmul(P(:,:,1,iEns),J(:,:,1)))
!      EH(2,iEns) = kappa(iEns)*trace_matrix(nBas,matmul(P(:,:,1,iEns),J(:,:,2))) &
!                 + kappa(iEns)*trace_matrix(nBas,matmul(P(:,:,2,iEns),J(:,:,1)))
!      EH(3,iEns) = kappa(iEns)*trace_matrix(nBas,matmul(P(:,:,2,iEns),J(:,:,2)))

!    else

      LZH(1,iEns) = - 0.5d0*trace_matrix(nBas,matmul(Pw(:,:,1),J(:,:,1)))
      LZH(2,iEns) = - 0.5d0*trace_matrix(nBas,matmul(Pw(:,:,1),J(:,:,2))) &
             - 0.5d0*trace_matrix(nBas,matmul(Pw(:,:,2),J(:,:,1)))
      LZH(3,iEns) = - 0.5d0*trace_matrix(nBas,matmul(Pw(:,:,2),J(:,:,2)))

      EH(1,iEns) = trace_matrix(nBas,matmul(P(:,:,1,iEns),J(:,:,1)))
      EH(2,iEns) = trace_matrix(nBas,matmul(P(:,:,1,iEns),J(:,:,2))) &
                 + trace_matrix(nBas,matmul(P(:,:,2,iEns),J(:,:,1)))
      EH(3,iEns) = trace_matrix(nBas,matmul(P(:,:,2,iEns),J(:,:,2)))

!    endif

  end do

end subroutine unrestricted_hartree_individual_energy
