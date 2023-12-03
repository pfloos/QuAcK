subroutine GTeh_excitation_density(nBas,nC,nO,nR,nS,ERI,XpY,XmY,rhoL,rhoR)

! Compute excitation densities

  implicit none

! Input variables

  integer,intent(in)            :: nBas,nC,nO,nR,nS
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: XpY(nS,nS)
  double precision,intent(in)   :: XmY(nS,nS)

! Local variables

  integer                       :: m,jb,p,q,j,b
  double precision              :: X,Y

! Output variables

  double precision,intent(out)  :: rhoL(nBas,nBas,nS)
  double precision,intent(out)  :: rhoR(nBas,nBas,nS)

! Initialization
 
  rhoL(:,:,:) = 0d0   
  rhoR(:,:,:) = 0d0   

  !$OMP PARALLEL &
  !$OMP SHARED(nC,nBas,nR,nO,nS,rhoL,rhoR,ERI,XpY,XmY) &
  !$OMP PRIVATE(q,p,jb,m,X,Y) &
  !$OMP DEFAULT(NONE)
  !$OMP DO
  do q=nC+1,nBas-nR
     do p=nC+1,nBas-nR
        jb = 0
        do j=nC+1,nO
           do b=nO+1,nBas-nR
              jb = jb + 1
              do m=1,nS

                X = 0.5d0*(XpY(m,jb) + XmY(m,jb))
                Y = 0.5d0*(XpY(m,jb) - XmY(m,jb))

                rhoL(p,q,m) = rhoL(p,q,m) + ERI(p,j,b,q)*X + ERI(p,b,j,q)*Y

                rhoR(p,q,m) = rhoR(p,q,m) + (2d0*ERI(p,j,b,q) - ERI(p,j,q,b))*X + (2d0*ERI(p,b,j,q) - ERI(p,b,q,j))*Y

              end do
           end do
        end do
     end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

end subroutine 
