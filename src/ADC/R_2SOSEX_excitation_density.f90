subroutine R_2SOSEX_excitation_density(nOrb,nC,nO,nR,nS,e,Om,ERI,XpY,rho)

! Compute excitation densities for 2SOSEX-spd

  implicit none

! Input variables

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
  double precision, allocatable :: tmp(:,:,:)
  double precision,allocatable  :: w(:,:,:)

! Output variables

  double precision,intent(out)  :: rho(nOrb,nOrb,nS)

! Memory allocation

  allocate(w(nOrb,nOrb,nS))

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

!---------------!
! Intermediates !
!---------------!

  w(:,:,:) = 0d0

  do p=nC+1,nOrb-nR
    do k=nC+1,nO
      do m=1,nS

        do j=nC+1,nO
          do a=nO+1,nOrb-nR

            w(p,k,m) = w(p,k,m) + rho(a,j,m)*(ERI(p,j,a,k)/(e(a) - e(j) + Om(m)) + ERI(p,a,j,k)/(e(a) - e(j) - Om(m)))

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

            w(p,c,m) = w(p,c,m) + rho(a,j,m)*(ERI(p,j,a,c)/(e(a) - e(j) - Om(m)) + ERI(p,a,j,c)/(e(a) - e(j) + Om(m)))

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
