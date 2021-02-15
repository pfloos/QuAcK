subroutine xc_potential(nBas,c,Fx,Fc,Vxc)

! Compute the exchange-correlation potential in the MO basis

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  double precision,intent(in)   :: c(nBas,nBas,nspin)
  double precision,intent(in)   :: Fx(nBas,nBas,nspin)
  double precision,intent(in)   :: Fc(nBas,nBas,nspin)

! Local variables

  integer                       :: mu,nu
  integer                       :: p
  integer                       :: ispin

! Output variables

  double precision,intent(out)  :: Vxc(nBas,nspin)

! Compute Vxc 

  Vxc(:,:) = 0d0
  do p=1,nBas
    do ispin=1,nspin
      do mu=1,nBas
        do nu=1,nBas
          Vxc(p,ispin) = Vxc(p,ispin) & 
                       + c(mu,p,ispin)*(Fx(mu,nu,ispin) + Fc(mu,nu,ispin))*c(nu,p,ispin)
    
        end do
      end do
    end do
  end do

end subroutine xc_potential
