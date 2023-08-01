subroutine UGF2_reg_self_energy_diag(nBas,nC,nO,nV,nR,eta,ERI_aa,ERI_ab,ERI_bb,eHF,eGF2,SigC,Z)

! Perform unrestricted GF2 self-energy and its renormalization factor

  implicit none
  include 'parameters.h'


! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC(nspin)
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nV(nspin)
  integer,intent(in)            :: nR(nspin)
  double precision,intent(in)   :: eta
  double precision,intent(in)   :: ERI_aa(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_ab(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_bb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: eHF(nBas,nspin)
  double precision,intent(in)   :: eGF2(nBas,nspin)

! Local variables

  integer                       :: p
  integer                       :: i,j,a,b
  double precision              :: eps,num

  double precision              :: s
  double precision              :: kappa

! Output variables

  double precision,intent(out)  :: SigC(nBas,nspin)
  double precision,intent(out)  :: Z(nBas,nspin)

!---------------------!
! Compute self-energy |
!---------------------!

  SigC(:,:) = 0d0

!-----------------------------------------!
! Parameters for regularized calculations !
!-----------------------------------------!

  s = 100d0

  !----------------!
  ! Spin-up sector
  !----------------!

  do p=nC(1)+1,nBas-nR(1)

    ! Addition part: aa

    do i=nC(1)+1,nO(1)
      do a=nO(1)+1,nBas-nR(1)
        do b=nO(1)+1,nBas-nR(1)

          eps = eGF2(p,1) + eHF(i,1) - eHF(a,1) - eHF(b,1) 
          num = ERI_aa(i,p,a,b)*ERI_aa(a,b,i,p) &
              - ERI_aa(i,p,a,b)*ERI_aa(a,b,p,i)
          kappa = 1d0 - exp(-2d0*eps**2*s)
          num = kappa*num

          SigC(p,1) = SigC(p,1) + num*eps/(eps**2 + eta**2)
          Z(p,1)    = Z(p,1)    - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2

        enddo
      enddo
    enddo

    ! Addition part: ab

    do i=nC(2)+1,nO(2)
      do a=nO(2)+1,nBas-nR(2)
        do b=nO(1)+1,nBas-nR(1)

          eps = eGF2(p,1) + eHF(i,2) - eHF(a,2) - eHF(b,1) 
          num = ERI_ab(p,i,b,a)*ERI_ab(b,a,p,i)
          kappa = 1d0 - exp(-2d0*eps**2*s)
          num = kappa*num

          SigC(p,1) = SigC(p,1) + num*eps/(eps**2 + eta**2)
          Z(p,1)    = Z(p,1)    - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2

        enddo
      enddo
    enddo

   ! Removal part: aa

    do i=nC(1)+1,nO(1)
      do a=nO(1)+1,nBas-nR(1)
        do j=nC(1)+1,nO(1)

          eps = eGF2(p,1) + eHF(a,1) - eHF(i,1) - eHF(j,1) 
          num = ERI_aa(a,p,i,j)*ERI_aa(i,j,a,p) &
              - ERI_aa(a,p,i,j)*ERI_aa(i,j,p,a)
          kappa = 1d0 - exp(-2d0*eps**2*s)
          num = kappa*num

          SigC(p,1) = SigC(p,1) + num*eps/(eps**2 + eta**2)
          Z(p,1)    = Z(p,1)    - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2

        enddo
      enddo
    enddo

   ! Removal part: ab

    do i=nC(2)+1,nO(2)
      do a=nO(2)+1,nBas-nR(2)
        do j=nC(1)+1,nO(1)

          eps = eGF2(p,1) + eHF(a,2) - eHF(i,2) - eHF(j,1) 
          num = ERI_ab(p,a,j,i)*ERI_ab(j,i,p,a)
          kappa = 1d0 - exp(-2d0*eps**2*s)
          num = kappa*num

          SigC(p,1) = SigC(p,1) + num*eps/(eps**2 + eta**2)
          Z(p,1)    = Z(p,1)    - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2

        enddo
      enddo
    enddo

  enddo

  !------------------!
  ! Spin-down sector !
  !------------------!

  do p=nC(2)+1,nBas-nR(2)

    ! Addition part: bb

    do i=nC(2)+1,nO(2)
      do a=nO(2)+1,nBas-nR(2)
        do b=nO(2)+1,nBas-nR(2)

          eps = eGF2(p,2) + eHF(i,2) - eHF(a,2) - eHF(b,2) 
          num = ERI_bb(i,p,a,b)*ERI_bb(a,b,i,p) &
              - ERI_bb(i,p,a,b)*ERI_bb(a,b,p,i)
          kappa = 1d0 - exp(-2d0*eps**2*s)
          num = kappa*num
       
          SigC(p,2) = SigC(p,2) + num*eps/(eps**2 + eta**2)
          Z(p,2)    = Z(p,2)    - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2

        enddo
      enddo
    enddo

    ! Addition part: ab

    do i=nC(1)+1,nO(1)
      do a=nO(1)+1,nBas-nR(1)
        do b=nO(2)+1,nBas-nR(2)

          eps = eGF2(p,2) + eHF(i,1) - eHF(a,1) - eHF(b,2) 
          num = ERI_ab(i,p,a,b)*ERI_ab(a,b,i,p)
          kappa = 1d0 - exp(-2d0*eps**2*s)
          num = kappa*num
       
          SigC(p,2) = SigC(p,2) + num*eps/(eps**2 + eta**2)
          Z(p,2)    = Z(p,2)    - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2

        enddo
      enddo
    enddo

   ! Removal part: bb

    do i=nC(2)+1,nO(2)
      do a=nO(2)+1,nBas-nR(2)
        do j=nC(2)+1,nO(2)

          eps = eGF2(p,2) + eHF(a,2) - eHF(i,2) - eHF(j,2) 
          num = ERI_bb(a,p,i,j)*ERI_bb(i,j,a,p) &
              - ERI_bb(a,p,i,j)*ERI_bb(i,j,p,a)
          kappa = 1d0 - exp(-2d0*eps**2*s)
          num = kappa*num

          SigC(p,2) = SigC(p,2) + num*eps/(eps**2 + eta**2)
          Z(p,2) = Z(p,2) - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2

        enddo
      enddo
    enddo

   ! Removal part: ab

    do i=nC(1)+1,nO(1)
      do a=nO(1)+1,nBas-nR(1)
        do j=nC(2)+1,nO(2)

          eps = eGF2(p,2) + eHF(a,1) - eHF(i,1) - eHF(j,2) 
          num = ERI_ab(a,p,i,j)*ERI_ab(i,j,a,p)
          kappa = 1d0 - exp(-2d0*eps**2*s)
          num = kappa*num

          SigC(p,2) = SigC(p,2) + num*eps/(eps**2 + eta**2)
          Z(p,2) = Z(p,2) - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2

        enddo
      enddo
    enddo

  enddo

  Z(:,:) = 1d0/(1d0 - Z(:,:))

end subroutine 
