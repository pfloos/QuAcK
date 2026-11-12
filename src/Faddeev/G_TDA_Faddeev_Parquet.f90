subroutine G_TDA_Faddeev_Parquet(nOrb,nC,nO,nR,nS,nOO,nVV,             &
                                  eQP,                                      &
                                  eh_Om,eh_X,eh_rho,                       &
                                  ee_Om,X1,ee_rho,                         &
                                  hh_Om,Y2,hh_rho,                         &
                                  nFm,nFp,                                  &
                                  MFm,eFm,VLm,VRm,                         &
                                  MFp,eFp,VLp,VRp)

! ======================================================================
! TDA Faddeev-Parquet intermediate-state matrices
! ======================================================================
!
! Generalized spin-orbital implementation.
!
! Removal sectors:
!
!   1: eh_1^- = (j,mu)    interacting pair (i,a), spectator j
!   2: eh_2^- = (i,mu)    interacting pair (j,a), spectator i
!   3: pp^-   = (a,nu)    hh pair (i,j), spectator a
!
! Addition sectors:
!
!   1: eh_1^+ = (b,mu)    interacting pair (i,a), spectator b
!   2: eh_2^+ = (a,mu)    interacting pair (i,b), spectator a
!   3: pp^+   = (i,lam)   ee pair (a,b), spectator i
!
! The TDA Faddeev matrices are
!
!                  ( K_eh   C_12   C_13 )
!     M_F^- / + =  ( C_21   K_eh   C_23 )
!                  ( C_31   C_32   K_pp )
!
! All directed off-diagonal blocks are constructed with the screened
! parquet transition densities and the TDA pair amplitudes.
!
! IMPORTANT:
!
!   * The amplitudes are used directly.  No division by
!       Omega - Delta epsilon
!     is performed.
!
!   * phGLR stores XpY with the mode index first.  In TDA,
!       XpY = XmY = X, so XpY can be passed directly as eh_X.
!
!   * In ppGLR, the TDA forward amplitudes are
!       ee : X1
!       hh : Y2.
!
!   * The matrices are non-Hermitian in general.
!
!   * This routine stops at the Faddeev eigensystem.  The physical
!     2h1p/2p1h projection and the one-body self-energy should be built
!     in a separate layer after the ADC(3) limiting-case test is passed.
!
! ======================================================================

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nOrb,nC,nO,nR,nS,nOO,nVV
  integer,intent(in)            :: nFm,nFp

  double precision,intent(in)   :: eQP(nOrb)

  double precision,intent(in)   :: eh_Om(nS)
  double precision,intent(in)   :: eh_X(nS,nS)
  double precision,intent(in)   :: eh_rho(nOrb,nOrb,nS)

  double precision,intent(in)   :: ee_Om(nVV)
  double precision,intent(in)   :: X1(nVV,nVV)
  double precision,intent(in)   :: ee_rho(nOrb,nOrb,nVV)

  double precision,intent(in)   :: hh_Om(nOO)
  double precision,intent(in)   :: Y2(nOO,nOO)
  double precision,intent(in)   :: hh_rho(nOrb,nOrb,nOO)

! Output variables

  double precision,intent(out)  :: MFm(nFm,nFm)
  double precision,intent(out)  :: eFm(nFm)
  double precision,intent(out)  :: VLm(nFm,nFm)
  double precision,intent(out)  :: VRm(nFm,nFm)

  double precision,intent(out)  :: MFp(nFp,nFp)
  double precision,intent(out)  :: eFp(nFp)
  double precision,intent(out)  :: VLp(nFp,nFp)
  double precision,intent(out)  :: VRp(nFp,nFp)

! Local variables

  integer                       :: nOcc,nVir
  integer                       :: nEhm,nHhm,nEhp,nEep
  integer                       :: nFm_ref,nFp_ref

  integer                       :: i,j,k,l
  integer                       :: a,b,c,d
  integer                       :: mu,mup,nu,lam
  integer                       :: ia,ij,ab,kl,cd
  integer                       :: row,col
  integer                       :: pEH

  double precision              :: val
  double precision,allocatable  :: eh_contract(:,:)
  double precision,allocatable  :: Am(:,:),Ap(:,:)

  integer,allocatable           :: eh_pair(:,:)
  integer,allocatable           :: hh_pair(:,:)
  integer,allocatable           :: ee_pair(:,:)

! ======================================================================
! Dimensions
! ======================================================================

  nOcc = nO - nC
  nVir = nOrb - nR - nO

  nEhm = nOcc*nS
  nHhm = nVir*nOO

  nEhp = nVir*nS
  nEep = nOcc*nVV

  nFm_ref = 2*nEhm + nHhm
  nFp_ref = 2*nEhp + nEep

  if(nFm /= nFm_ref) then
    write(*,*) 'ERROR in G_TDA_Faddeev_Parquet: wrong removal dimension'
    write(*,*) 'nFm supplied = ',nFm
    write(*,*) 'nFm expected = ',nFm_ref
    stop
  end if

  if(nFp /= nFp_ref) then
    write(*,*) 'ERROR in G_TDA_Faddeev_Parquet: wrong addition dimension'
    write(*,*) 'nFp supplied = ',nFp
    write(*,*) 'nFp expected = ',nFp_ref
    stop
  end if

! ======================================================================
! TDA convention
!
! The caller must provide TDA response amplitudes:
!
!   eh_X(mu,ia) = X^eh_{ia,mu}
!   X1(ab,lam)  = X^ee_{ab,lam}
!   Y2(ij,nu)   = Y^hh_{ij,nu}
!
! In the current G_evParquet implementation, eh_X can be supplied
! directly as XpY when TDAeh=.true.
! ======================================================================

! ======================================================================
! Pair-index lookup tables
! ======================================================================

  allocate(eh_pair(nOrb,nOrb))
  allocate(hh_pair(nOrb,nOrb))
  allocate(ee_pair(nOrb,nOrb))

  eh_pair(:,:) = 0
  hh_pair(:,:) = 0
  ee_pair(:,:) = 0

! eh pair: occupied i, virtual a

  ia = 0
  do i=nC+1,nO
    do a=nO+1,nOrb-nR
      ia = ia + 1
      eh_pair(i,a) = ia
    end do
  end do

  if(ia /= nS) then
    write(*,*) 'ERROR in G_TDA_Faddeev_Parquet: inconsistent eh dimension'
    write(*,*) 'number of active ia pairs = ',ia
    write(*,*) 'nS                        = ',nS
    stop
  end if

! hh pair: i < j

  ij = 0
  do i=nC+1,nO
    do j=i+1,nO
      ij = ij + 1
      hh_pair(i,j) = ij
    end do
  end do

  if(ij /= nOO) then
    write(*,*) 'ERROR in G_TDA_Faddeev_Parquet: inconsistent hh dimension'
    write(*,*) 'number of active ij pairs = ',ij
    write(*,*) 'nOO                       = ',nOO
    stop
  end if

! ee pair: a < b

  ab = 0
  do a=nO+1,nOrb-nR
    do b=a+1,nOrb-nR
      ab = ab + 1
      ee_pair(a,b) = ab
    end do
  end do

  if(ab /= nVV) then
    write(*,*) 'ERROR in G_TDA_Faddeev_Parquet: inconsistent ee dimension'
    write(*,*) 'number of active ab pairs = ',ab
    write(*,*) 'nVV                       = ',nVV
    stop
  end if

! ======================================================================
! Repeated eh contraction
!
!   eh_contract(mu,mup)
!     = sum_{kc} rho^eh_{kc,mu} X^eh_{kc,mup}
!
! This is only an implementation intermediate.
! ======================================================================

  allocate(eh_contract(nS,nS))
  eh_contract(:,:) = 0d0

  do mu=1,nS
    do mup=1,nS

      val = 0d0

      do k=nC+1,nO
        do c=nO+1,nOrb-nR
          pEH = eh_pair(k,c)
          val = val + eh_rho(k,c,mu)*eh_X(mup,pEH)
        end do
      end do

      eh_contract(mu,mup) = val

    end do
  end do

! ======================================================================
! Initialize matrices
! ======================================================================

  MFm(:,:) = 0d0
  MFp(:,:) = 0d0

! ======================================================================
! Removal diagonal blocks
! ======================================================================

  do j=nC+1,nO
    do mu=1,nS
      row = idx_m1(j,mu)
      MFm(row,row) = eQP(j) - eh_Om(mu)
    end do
  end do

  do i=nC+1,nO
    do mu=1,nS
      row = idx_m2(i,mu)
      MFm(row,row) = eQP(i) - eh_Om(mu)
    end do
  end do

  do a=nO+1,nOrb-nR
    do nu=1,nOO
      row = idx_m3(a,nu)
      MFm(row,row) = hh_Om(nu) - eQP(a)
    end do
  end do

! ======================================================================
! Removal C_12^-
! ======================================================================

  do j=nC+1,nO
    do mu=1,nS

      row = idx_m1(j,mu)

      do i=nC+1,nO
        do mup=1,nS

          col = idx_m2(i,mup)
          val = 0d0

          if(i == j) val = val + eh_contract(mu,mup)

          do c=nO+1,nOrb-nR
            val = val - eh_rho(i,c,mu)*eh_X(mup,eh_pair(j,c))
          end do

          MFm(row,col) = val

        end do
      end do

    end do
  end do

! ======================================================================
! Removal C_21^-
! ======================================================================

  do i=nC+1,nO
    do mu=1,nS

      row = idx_m2(i,mu)

      do j=nC+1,nO
        do mup=1,nS

          col = idx_m1(j,mup)
          val = 0d0

          if(i == j) val = val + eh_contract(mu,mup)

          do c=nO+1,nOrb-nR
            val = val - eh_rho(j,c,mu)*eh_X(mup,eh_pair(i,c))
          end do

          MFm(row,col) = val

        end do
      end do

    end do
  end do

! ======================================================================
! Removal C_13^-
! ======================================================================

  do j=nC+1,nO
    do mu=1,nS

      row = idx_m1(j,mu)

      do b=nO+1,nOrb-nR
        do nu=1,nOO

          col = idx_m3(b,nu)
          val = 0d0

          kl = 0
          do k=nC+1,nO
            do l=k+1,nO
              kl = kl + 1

              if(j == k) val = val + eh_rho(l,b,mu)*Y2(kl,nu)
              if(j == l) val = val - eh_rho(k,b,mu)*Y2(kl,nu)

            end do
          end do

          MFm(row,col) = val

        end do
      end do

    end do
  end do

! ======================================================================
! Removal C_23^-
! ======================================================================

  do i=nC+1,nO
    do mu=1,nS

      row = idx_m2(i,mu)

      do b=nO+1,nOrb-nR
        do nu=1,nOO

          col = idx_m3(b,nu)
          val = 0d0

          kl = 0
          do k=nC+1,nO
            do l=k+1,nO
              kl = kl + 1

              if(i == k) val = val - eh_rho(l,b,mu)*Y2(kl,nu)
              if(i == l) val = val + eh_rho(k,b,mu)*Y2(kl,nu)

            end do
          end do

          MFm(row,col) = val

        end do
      end do

    end do
  end do

! ======================================================================
! Removal C_31^-
! ======================================================================

  do a=nO+1,nOrb-nR
    do nu=1,nOO

      row = idx_m3(a,nu)

      do j=nC+1,nO
        do mu=1,nS

          col = idx_m1(j,mu)
          val = 0d0

          do k=nC+1,nO
            val = val - hh_rho(k,j,nu)*eh_X(mu,eh_pair(k,a))
          end do

          MFm(row,col) = val

        end do
      end do

    end do
  end do

! ======================================================================
! Removal C_32^-
! ======================================================================

  do a=nO+1,nOrb-nR
    do nu=1,nOO

      row = idx_m3(a,nu)

      do i=nC+1,nO
        do mu=1,nS

          col = idx_m2(i,mu)
          val = 0d0

          do l=nC+1,nO
            val = val - hh_rho(i,l,nu)*eh_X(mu,eh_pair(l,a))
          end do

          MFm(row,col) = val

        end do
      end do

    end do
  end do

! ======================================================================
! Addition diagonal blocks
! ======================================================================

  do b=nO+1,nOrb-nR
    do mu=1,nS
      row = idx_p1(b,mu)
      MFp(row,row) = eQP(b) + eh_Om(mu)
    end do
  end do

  do a=nO+1,nOrb-nR
    do mu=1,nS
      row = idx_p2(a,mu)
      MFp(row,row) = eQP(a) + eh_Om(mu)
    end do
  end do

  do i=nC+1,nO
    do lam=1,nVV
      row = idx_p3(i,lam)
      MFp(row,row) = ee_Om(lam) - eQP(i)
    end do
  end do

! ======================================================================
! Addition C_12^+
! ======================================================================

  do b=nO+1,nOrb-nR
    do mu=1,nS

      row = idx_p1(b,mu)

      do a=nO+1,nOrb-nR
        do mup=1,nS

          col = idx_p2(a,mup)
          val = 0d0

          do k=nC+1,nO
            val = val + eh_rho(k,a,mu)*eh_X(mup,eh_pair(k,b))
          end do

          if(a == b) val = val - eh_contract(mu,mup)

          MFp(row,col) = val

        end do
      end do

    end do
  end do

! ======================================================================
! Addition C_21^+
! ======================================================================

  do a=nO+1,nOrb-nR
    do mu=1,nS

      row = idx_p2(a,mu)

      do b=nO+1,nOrb-nR
        do mup=1,nS

          col = idx_p1(b,mup)
          val = 0d0

          do k=nC+1,nO
            val = val + eh_rho(k,b,mu)*eh_X(mup,eh_pair(k,a))
          end do

          if(a == b) val = val - eh_contract(mu,mup)

          MFp(row,col) = val

        end do
      end do

    end do
  end do

! ======================================================================
! Addition C_13^+
! ======================================================================

  do b=nO+1,nOrb-nR
    do mu=1,nS

      row = idx_p1(b,mu)

      do i=nC+1,nO
        do lam=1,nVV

          col = idx_p3(i,lam)
          val = 0d0

          cd = 0
          do c=nO+1,nOrb-nR
            do d=c+1,nOrb-nR
              cd = cd + 1

              if(b == d) val = val + eh_rho(i,c,mu)*X1(cd,lam)
              if(b == c) val = val - eh_rho(i,d,mu)*X1(cd,lam)

            end do
          end do

          MFp(row,col) = val

        end do
      end do

    end do
  end do

! ======================================================================
! Addition C_23^+
! ======================================================================

  do a=nO+1,nOrb-nR
    do mu=1,nS

      row = idx_p2(a,mu)

      do i=nC+1,nO
        do lam=1,nVV

          col = idx_p3(i,lam)
          val = 0d0

          cd = 0
          do c=nO+1,nOrb-nR
            do d=c+1,nOrb-nR
              cd = cd + 1

              if(a == c) val = val + eh_rho(i,d,mu)*X1(cd,lam)
              if(a == d) val = val - eh_rho(i,c,mu)*X1(cd,lam)

            end do
          end do

          MFp(row,col) = val

        end do
      end do

    end do
  end do

! ======================================================================
! Addition C_31^+
! ======================================================================

  do i=nC+1,nO
    do lam=1,nVV

      row = idx_p3(i,lam)

      do b=nO+1,nOrb-nR
        do mu=1,nS

          col = idx_p1(b,mu)
          val = 0d0

          do c=nO+1,nOrb-nR
            val = val + ee_rho(c,b,lam)*eh_X(mu,eh_pair(i,c))
          end do

          MFp(row,col) = val

        end do
      end do

    end do
  end do

! ======================================================================
! Addition C_32^+
! ======================================================================

  do i=nC+1,nO
    do lam=1,nVV

      row = idx_p3(i,lam)

      do a=nO+1,nOrb-nR
        do mu=1,nS

          col = idx_p2(a,mu)
          val = 0d0

          do d=nO+1,nOrb-nR
            val = val + ee_rho(a,d,lam)*eh_X(mu,eh_pair(i,d))
          end do

          MFp(row,col) = val

        end do
      end do

    end do
  end do

! ======================================================================
! Non-Hermitian diagonalization
! ======================================================================

  allocate(Am(nFm,nFm))
  allocate(Ap(nFp,nFp))

  Am(:,:) = MFm(:,:)
  Ap(:,:) = MFp(:,:)

  call diagonalize_general_matrix_LR(nFm,Am,eFm,VLm,VRm)
  call diagonalize_general_matrix_LR(nFp,Ap,eFp,VLp,VRp)

! ======================================================================
! Summary
! ======================================================================

  write(*,*)
  write(*,*)'***********************************************'
  write(*,*)'* Generalized TDA Faddeev-Parquet matrices    *'
  write(*,*)'***********************************************'
  write(*,'(1X,A28,I10)') 'Removal Faddeev dimension =',nFm
  write(*,'(1X,A28,I10)') 'Addition Faddeev dimension=',nFp
  write(*,*)
  write(*,*)'TDA amplitudes used directly: eh_X, X1, Y2'
  write(*,*)
  write(*,*)'Next validation step: project to the antisymmetric'
  write(*,*)'2h1p/2p1h spaces and compare element by element'
  write(*,*)'with the ADC(3) K+C blocks.'
  write(*,*)

! ======================================================================
! Clean up
! ======================================================================

  deallocate(eh_contract)
  deallocate(eh_pair,hh_pair,ee_pair)
  deallocate(Am,Ap)

contains

! ======================================================================
! Compact Faddeev indices
! ======================================================================

  integer function idx_m1(j,mu)

    implicit none
    integer,intent(in) :: j,mu

    idx_m1 = (j-nC-1)*nS + mu

  end function idx_m1


  integer function idx_m2(i,mu)

    implicit none
    integer,intent(in) :: i,mu

    idx_m2 = nEhm + (i-nC-1)*nS + mu

  end function idx_m2


  integer function idx_m3(a,nu)

    implicit none
    integer,intent(in) :: a,nu

    idx_m3 = 2*nEhm + (a-nO-1)*nOO + nu

  end function idx_m3


  integer function idx_p1(b,mu)

    implicit none
    integer,intent(in) :: b,mu

    idx_p1 = (b-nO-1)*nS + mu

  end function idx_p1


  integer function idx_p2(a,mu)

    implicit none
    integer,intent(in) :: a,mu

    idx_p2 = nEhp + (a-nO-1)*nS + mu

  end function idx_p2


  integer function idx_p3(i,lam)

    implicit none
    integer,intent(in) :: i,lam

    idx_p3 = 2*nEhp + (i-nC-1)*nVV + lam

  end function idx_p3

end subroutine G_TDA_Faddeev_Parquet
