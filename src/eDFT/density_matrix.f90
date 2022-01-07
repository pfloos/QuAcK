subroutine density_matrix(nBas,nEns,c,P,occnum)

! Calculate density matrices

  implicit none

  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: c(nBas,nBas,nspin)
  double precision,intent(in)   :: occnum(nBas,nspin,nEns)


! Local variables

  integer                       :: ispin
  integer                       :: iEns
  integer                       :: q
  integer                       :: mu,nu

! Output variables

  double precision,intent(out)  :: P(nBas,nBas,nspin,nEns)

! Compute density matrix for each state of the ensemble based on occupation numbers

  P(:,:,:,:) = 0d0
  do iEns=1,nEns
    do ispin=1,nspin
      do mu=1,nBas
        do nu=1,nBas
          do q=1,nBas

            P(mu,nu,ispin,iEns) = P(mu,nu,ispin,iEns) & 
                                + occnum(q,ispin,iEns)*c(mu,q,ispin)*c(nu,q,ispin)

          end do 
        end do 
      end do 
    end do 
  end do



end subroutine density_matrix
