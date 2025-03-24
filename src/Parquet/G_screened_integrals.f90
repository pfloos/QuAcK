subroutine G_eh_screened_integral(nOrb,nC,nO,nR,nS,ERI,eh_Gam,XpY,XmY,rho)

! Compute excitation densities
  implicit none

! Input variables
  integer,intent(in)            :: nOrb,nC,nO,nR,nS
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: eh_Gam(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: XpY(nS,nS),XmY(nS,nS)

! Local variables
  integer                       :: ia,jb,p,q,j,b
  double precision              :: X,Y

! Output variables
  double precision,intent(out)  :: rho(nOrb,nOrb,nS)
  
  rho(:,:,:) = 0d0   
!  !$OMP PARALLEL &
!  !$OMP SHARED(nC,nOrb,nR,nO,nS,rho,ERI,XpY,eh_sing_Gam) &
!  !$OMP PRIVATE(q,p,jb,ia) &
!  !$OMP DEFAULT(NONE)
!  !$OMP DO
  do q=nC+1,nOrb-nR
     do p=nC+1,nOrb-nR

        jb = 0
        do j=nC+1,nO
           do b=nO+1,nOrb-nR
              jb = jb + 1
              
              do ia=1,nS

                 X = 0.5d0*(XpY(ia,jb) + XmY(ia,jb))
                 Y = 0.5d0*(XpY(ia,jb) - XmY(ia,jb))

                 rho(p,q,ia) = rho(p,q,ia)                              &
                             !+ (ERI(p,j,q,b) - ERI(p,j,b,q))*XpY(ia,jb) &
                             + (ERI(p,j,q,b) - ERI(p,j,b,q))*X &
                             + (ERI(p,b,q,j) - ERI(p,b,j,q))*Y &
                             + 1d0*eh_Gam(p,j,q,b)*X &
                             + 1d0*eh_Gam(p,b,q,j)*Y
                 

              end do
           end do
        end do
     end do
  end do
!  !$OMP END DO
!  !$OMP END PARALLEL
  
end subroutine G_eh_screened_integral

subroutine G_pp_screened_integral(nOrb,nC,nO,nV,nR,nOO,nVV,ERI,pp_Gam,X1,Y1,rho1,X2,Y2,rho2)

! Compute excitation densities in the singlet pp channel

  implicit none

! Input variables

  integer,intent(in)            :: nOrb,nC,nO,nV,nR
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: pp_Gam(nOrb,nOrb,nOrb,nOrb)
  integer,intent(in)            :: nOO
  integer,intent(in)            :: nVV 
  double precision,intent(in)   :: X1(nVV,nVV)
  double precision,intent(in)   :: Y1(nOO,nVV)
  double precision,intent(in)   :: X2(nVV,nOO)
  double precision,intent(in)   :: Y2(nOO,nOO)

! Local variables

  integer                       :: i,j,k,l
  integer                       :: a,b,c,d
  integer                       :: p,q
  integer                       :: ab,cd,ij,kl
  double precision,external     :: Kronecker_delta

! Output variables

  double precision,intent(out)  :: rho1(nOrb,nOrb,nVV)
  double precision,intent(out)  :: rho2(nOrb,nOrb,nOO)

  integer                       :: dim_1, dim_2

! Initialization

  rho1(:,:,:) = 0d0
  rho2(:,:,:) = 0d0

!  !$OMP PARALLEL DEFAULT(NONE)                                         &
!  !$OMP          PRIVATE(p, q, a, b, ab, c, d, cd, i, j, ij, k, l, kl) &
!  !$OMP          SHARED(nC, nOrb, nR, nO, rho1, rho2, ERI, pp_sing_Gam, X1, Y1, X2, Y2)
!  !$OMP DO COLLAPSE(2)
  do q=nC+1,nOrb-nR
     do p=nC+1,nOrb-nR
           
        ab = 0
        do a=nO+1,nOrb-nR
           do b=a+1,nOrb-nR
              ab = ab + 1

              cd = 0
              do c=nO+1,nOrb-nR
                 do d=c+1,nOrb-nR
                    cd = cd + 1
                    rho1(p,q,ab) = rho1(p,q,ab) & 
                         + (ERI(p,q,c,d) - ERI(p,q,d,c) + 1d0*pp_Gam(p,q,c,d))*X1(cd,ab)
                 end do
              end do
          
              kl = 0
              do k=nC+1,nO
                 do l=k+1,nO
                    kl = kl + 1
                    rho1(p,q,ab) = rho1(p,q,ab) & 
                         + (ERI(p,q,k,l) - ERI(p,q,l,k) + 1d0*pp_Gam(p,q,k,l))*Y1(kl,ab)
                 end do
              end do
              
           end do
        end do

        ij = 0
        do i=nC+1,nO
           do j=i+1,nO
              ij = ij + 1
              cd = 0
              do c=nO+1,nOrb-nR
                 do d=c+1,nOrb-nR
                    cd = cd + 1
                    rho2(p,q,ij) = rho2(p,q,ij) &
                         + (ERI(p,q,c,d) - ERI(p,q,d,c) + 1d0*pp_Gam(p,q,c,d))*X2(cd,ij)
                 end do
              end do
              
              kl = 0
              do k=nC+1,nO
                 do l=k+1,nO
                    kl = kl + 1
                    rho2(p,q,ij) = rho2(p,q,ij) &
                         + (ERI(p,q,k,l) - ERI(p,q,l,k) + 1d0*pp_Gam(p,q,k,l))*Y2(kl,ij) 
                 end do
              end do
           end do
        end do
     end do
  end do
!  !$OMP END DO
!  !$OMP END PARALLEL

end subroutine G_pp_screened_integral
