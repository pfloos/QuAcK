subroutine self_energy_GF2(eta,nBas,nC,nO,nV,nR,nS,eHF,eGF2,ERI,SigC,Z,Ec)

! Compute GF2 self-energy and its renormalization factor

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas,nC,nO,nV,nR,nS
  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: eGF2(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: i,j,a,b
  integer                       :: p,q
  double precision              :: eps
  double precision              :: num

! Output variables

  double precision,intent(out)  :: SigC(nBas,nBas)
  double precision,intent(out)  :: Z(nBas)
  double precision,intent(out)  :: Ec

! Initialize 

  SigC(:,:) = 0d0
  Z(:)     = 0d0

! Compute GF2 self-energy and renormalization factor

  do p=nC+1,nBas-nR
    do q=nC+1,nBas-nR
      do i=nC+1,nO
        do j=nC+1,nO
          do a=nO+1,nBas-nR

            eps = eGF2(p) + eHF(a) - eHF(i) - eHF(j)
            num = (2d0*ERI(p,a,i,j) - ERI(p,a,j,i))*ERI(q,a,i,j)

            SigC(p,q) = SigC(p,q) + num*eps/(eps**2 + eta**2)
            if(p == q) Z(p)   = Z(p)   - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2

          end do
        end do
      end do
    end do
  end do

  do p=nC+1,nBas-nR
    do q=nC+1,nBas-nR
      do i=nC+1,nO
        do a=nO+1,nBas-nR
          do b=nO+1,nBas-nR

            eps = eGF2(p) + eHF(i) - eHF(a) - eHF(b)
            num = (2d0*ERI(p,i,a,b) - ERI(p,i,b,a))*ERI(q,i,a,b)

            SigC(p,q) = SigC(p,q) + num*eps/(eps**2 + eta**2)
            if(p == q) Z(p)   = Z(p)   - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2

          end do
        end do
      end do
    end do
  end do

  Z(:) = 1d0/(1d0 - Z(:))

! Compute correlaiton energy

  Ec = 0d0

  do j=nC+1,nO
    do i=nC+1,nO
      do a=nO+1,nBas-nR
        do b=nO+1,nBas-nR

          eps = eGF2(j) + eHF(i) - eHF(a) - eHF(b)
          num = (2d0*ERI(j,i,a,b) - ERI(j,i,b,a))*ERI(j,i,a,b)

          Ec = Ec + num*eps/(eps**2 + eta**2)

        end do
      end do
    end do
  end do

end subroutine self_energy_GF2
