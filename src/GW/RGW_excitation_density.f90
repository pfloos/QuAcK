subroutine RGW_excitation_density(nOrb,nC,nO,nR,nS,ERI,XpY,rho)

! Compute excitation densities

  implicit none

! Input variables

  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: XpY(nS,nS)

! Local variables

  integer                       :: ia,jb,p,q,j,b
  double precision, allocatable :: tmp(:,:,:)

! Output variables

  double precision,intent(out)  :: rho(nOrb,nOrb,nS)

  if(nOrb .lt. 256) then

    allocate(tmp(nOrb,nOrb,nS))

    !$OMP PARALLEL DEFAULT(NONE)             &
    !$OMP PRIVATE(p, q, j, b, jb)            &
    !$OMP SHARED(nOrb, nC, nO, nR, ERI, tmp)
    !$OMP DO COLLAPSE(2)
    do p = 1, nOrb
      do q = 1, nOrb
        jb = 0
        do j = nC+1, nO
          do b = nO+1, nOrb-nR
            jb = jb + 1
            tmp(p,q,jb) = ERI(p,j,q,b)
          enddo
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

    call dgemm("N", "T", nOrb*nOrb, nS, nS, 1.d0,   &
               tmp(1,1,1), nOrb*nOrb, XpY(1,1), nS, &
               0.d0, rho(1,1,1), nOrb*nOrb)

    deallocate(tmp)

  else

    rho(:,:,:) = 0d0   
    !$OMP PARALLEL &
    !$OMP SHARED(nC,nOrb,nR,nO,nS,rho,ERI,XpY) &
    !$OMP PRIVATE(q,p,jb,ia) &
    !$OMP DEFAULT(NONE)
    !$OMP DO
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
    !$OMP END DO
    !$OMP END PARALLEL

  endif

end subroutine 
