subroutine RGW_phBSE_static_kernel_fullmat(nBas,nC,nO,nV,nR,nS,ERI,e,Om,rho,W)

! Compute the second-order static BSE kernel for the resonant block (only for singlets!)

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS

  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: e(nBas)
  double precision,intent(in)   :: Om(nS)

  double precision,intent(in)   :: rho(nBas,nBas,nS)


! Local variables

  double precision              :: num,dem,reg
  double precision              :: flow = 1d6
  integer                       :: p,q,r,s
  integer                       :: m

! Output variables

  double precision,intent(out)   :: W(nBas,nBas,nBas,nBas)

!------------------------------------------------
! Compute static screening (physicist's notation)
!------------------------------------------------

  do p=nC+1,nBas-nR
    do q=nC+1,nBas-nR
      do r=nC+1,nBas-nR
        do s=nC+1,nBas-nR

          W(p,s,q,r) = ERI(p,s,q,r)

          do m=1,nS
            
            num = 2d0*rho(p,q,m)*rho(r,s,m)
            dem = e(s) - e(p) + Om(m)
            reg = (1d0 - exp(-2d0*flow*dem*dem))/dem

            W(p,s,q,r) = W(p,s,q,r) + num*reg

            num = 2d0*rho(p,q,m)*rho(r,s,m)
            dem = e(r) - e(q) + Om(m)
            reg = (1d0 - exp(-2d0*flow*dem*dem))/dem

            W(p,s,q,r) = W(p,s,q,r) + num*reg

          end do

        end do
      end do
    end do
  end do

end subroutine 
