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

! Output variables

  double precision,intent(out)  :: rhoL(nBas,nBas,nS)
  double precision,intent(out)  :: rhoR(nBas,nBas,nS)

  rhoL(:,:,:) = 0d0   
  rhoR(:,:,:) = 0d0   
  !$OMP PARALLEL &
  !$OMP SHARED(nC,nBas,nR,nO,nS,rhoL,rhoR,ERI,XpY,XmY) &
  !$OMP PRIVATE(q,p,jb,m) &
  !$OMP DEFAULT(NONE)
  !$OMP DO
  do q=nC+1,nBas-nR
     do p=nC+1,nBas-nR
        jb = 0
        do j=nC+1,nO
           do b=nO+1,nBas-nR
              jb = jb + 1
              do m=1,nS

                 rhoL(p,q,m) = rhoL(p,q,m)                                 & 
                              + ERI(p,j,b,q)*0.5d0*(XpY(m,jb) + XmY(m,jb)) &
                              + ERI(p,b,j,q)*0.5d0*(XpY(m,jb) - XmY(m,jb))

                 rhoR(p,q,m) = rhoR(p,q,m)                                                      & 
                              + (2d0*ERI(p,j,b,q) - ERI(p,j,q,b))*0.5d0*(XpY(m,jb) + XmY(m,jb)) &
                              + (2d0*ERI(p,b,j,q) - ERI(p,b,q,j))*0.5d0*(XpY(m,jb) - XmY(m,jb))
              enddo
           enddo
        enddo
     enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

end subroutine 
