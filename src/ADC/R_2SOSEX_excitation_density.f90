subroutine R_2SOSEX_excitation_density(flow,nOrb,nC,nO,nR,nS,e,Om,ERI,XpY,rho)

! Compute excitation densities for 2SOSEX-spd

  implicit none

! Input variables

  double precision,intent(in)   :: flow
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: e(nOrb)
  double precision,intent(in)   :: Om(nS)
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: XpY(nS,nS)

! Local variables

  integer                       :: ia,jb,p,q,j,a,b,k,c,m
  double precision              :: num,dem,reg
  double precision, allocatable :: tmp(:,:,:)
  double precision,allocatable  :: w(:,:,:)

! Output variables

  double precision,intent(out)  :: rho(nOrb,nOrb,nS)

! Memory allocation

  allocate(w(nOrb,nOrb,nS))

  rho(:,:,:) = 0d0   

  do q=nC+1,nOrb-nR
     do p=nC+1,nOrb-nR
        jb = 0
        do j=nC+1,nO
           do b=nO+1,nOrb-nR
              jb = jb + 1
              do ia=1,nS
                 rho(p,q,ia) = rho(p,q,ia) + ERI(p,j,q,b)*XpY(ia,jb)
              end do
           end do
        end do
     end do
  end do

!---------------!
! Intermediates !
!---------------!

  w(:,:,:) = 0d0

  do p=nC+1,nOrb-nR
    do k=nC+1,nO
      do m=1,nS

        do j=nC+1,nO
          do a=nO+1,nOrb-nR

            num = rho(a,j,m)*ERI(p,j,a,k)
            dem = e(a) - e(j) + Om(m)
            reg = (1d0 - exp(-2d0*flow*dem**2))/dem

            w(p,k,m) = w(p,k,m) + num*reg

            num = rho(a,j,m)*ERI(p,a,j,k)
            dem = e(a) - e(j) - Om(m)
            reg = (1d0 - exp(-2d0*flow*dem**2))/dem

            w(p,k,m) = w(p,k,m) + num*reg

          end do
        end do

      end do
    end do
  end do

  do p=nC+1,nOrb-nR
    do c=nO+1,nOrb-nR
      do m=1,nS

        do j=nC+1,nO
          do a=nO+1,nOrb-nR

            num = rho(a,j,m)*ERI(p,j,a,c)
            dem = e(a) - e(j) - Om(m)
            reg = (1d0 - exp(-2d0*flow*dem**2))/dem

            w(p,c,m) = w(p,c,m) + num*reg

            num = rho(a,j,m)*ERI(p,a,j,c)
            dem = e(a) - e(j) + Om(m)
            reg = (1d0 - exp(-2d0*flow*dem**2))/dem

            w(p,c,m) = w(p,c,m) + num*reg

          end do
        end do

      end do
    end do
  end do

! Final excitation densities for 2SOSEX-psd

  do m=1,nS
    do q=nC+1,nOrb-nR
      do p=nC+1,nOrb-nR

        rho(p,q,m) = rho(p,q,m) + w(p,q,m)

      end do
    end do
  end do


end subroutine 
