subroutine G_eh_screened_integral(nOrb,nC,nO,nR,nS,ERI,eh_Phi,pp_Phi,XpY,XmY,rho)

! Compute excitation densities
  implicit none

! Input variables
  integer,intent(in)            :: nOrb,nC,nO,nR,nS
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: eh_Phi(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: pp_Phi(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: XpY(nS,nS),XmY(nS,nS)

! Local variables
  integer                       :: ia,jb,p,q,j,b
  double precision              :: X,Y

! Output variables
  double precision,intent(out)  :: rho(nOrb,nOrb,nS)
  
  rho(:,:,:) = 0d0   
  !$OMP PARALLEL DEFAULT(NONE)                                 &
  !$OMP PRIVATE(q,p,j,b,jb,ia,X,Y)                             &
  !$OMP SHARED(nC,nOrb,nR,nO,nS,rho,ERI,XpY,XmY,eh_Phi,pp_Phi)
  !$OMP DO COLLAPSE(2)
  do q=nC+1,nOrb-nR
     do p=nC+1,nOrb-nR
              
        do ia=1,nS
           
           jb = 0
           do j=nC+1,nO
              do b=nO+1,nOrb-nR
                 jb = jb + 1

                 X = 0.5d0*(XpY(ia,jb) + XmY(ia,jb))
                 Y = 0.5d0*(XpY(ia,jb) - XmY(ia,jb))

                 rho(p,q,ia) = rho(p,q,ia) &
                             + (ERI(p,b,q,j) - ERI(p,b,j,q)) * X &
                             + (-eh_Phi(p,b,j,q) + pp_Phi(p,b,q,j)) * X &
                             + (ERI(p,j,q,b) - ERI(p,j,b,q)) * Y &
                             + (-eh_Phi(p,j,b,q) + pp_Phi(p,j,q,b)) * Y

              end do
           end do
           
        end do
        
     end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  
end subroutine 
