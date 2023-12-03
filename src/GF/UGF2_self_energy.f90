subroutine UGF2_self_energy(nBas,nC,nO,nV,nR,eta,ERI_aa,ERI_ab,ERI_bb,eHF,eGF2,SigC,Z)

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

  integer                       :: p,q
  integer                       :: i,j,a,b
  double precision              :: eps,num

! Output variables

  double precision,intent(out)  :: SigC(nBas,nBas,nspin)
  double precision,intent(out)  :: Z(nBas,nspin)

!---------------------!
! Compute self-energy |
!---------------------!

  SigC(:,:,:) = 0d0
  Z(:,:)      = 0d0

  !----------------!
  ! Spin-up sector
  !----------------!

  do p=nC(1)+1,nBas-nR(1)
    do q=nC(1)+1,nBas-nR(1)

      ! Addition part: aa

      do i=nC(1)+1,nO(1)
        do a=nO(1)+1,nBas-nR(1)
          do b=nO(1)+1,nBas-nR(1)

            eps = eGF2(p,1) + eHF(i,1) - eHF(a,1) - eHF(b,1) 
            num = ERI_aa(i,q,a,b)*ERI_aa(a,b,i,p) &
                - ERI_aa(i,q,a,b)*ERI_aa(a,b,p,i)
         
            SigC(p,q,1) = SigC(p,q,1) + num*eps/(eps**2 + eta**2)  
            if(p == q) Z(p,1) = Z(p,1) - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2

          end do
        end do
      end do

      ! Addition part: ab

      do i=nC(2)+1,nO(2)
        do a=nO(2)+1,nBas-nR(2)
          do b=nO(1)+1,nBas-nR(1)

            eps = eGF2(p,1) + eHF(i,2) - eHF(a,2) - eHF(b,1) 
            num = ERI_ab(q,i,b,a)*ERI_ab(b,a,p,i)
         
            SigC(p,q,1) = SigC(p,q,1) + num*eps/(eps**2 + eta**2)
            if(p == q) Z(p,1) = Z(p,1) - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2

          end do
        end do
      end do

     ! Removal part: aa

      do i=nC(1)+1,nO(1)
        do a=nO(1)+1,nBas-nR(1)
          do j=nC(1)+1,nO(1)

            eps = eGF2(p,1) + eHF(a,1) - eHF(i,1) - eHF(j,1) 
            num = ERI_aa(a,q,i,j)*ERI_aa(i,j,a,p) &
                - ERI_aa(a,q,i,j)*ERI_aa(i,j,p,a)
         
            SigC(p,q,1) = SigC(p,q,1) + num*eps/(eps**2 + eta**2)
            if(p == q) Z(p,1) = Z(p,1) - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2

          end do
        end do
      end do

     ! Removal part: ab

      do i=nC(2)+1,nO(2)
        do a=nO(2)+1,nBas-nR(2)
          do j=nC(1)+1,nO(1)

            eps = eGF2(p,1) + eHF(a,2) - eHF(i,2) - eHF(j,1) 
            num = ERI_ab(q,a,j,i)*ERI_ab(j,i,p,a)
         
            SigC(p,q,1) = SigC(p,q,1) + num*eps/(eps**2 + eta**2)
            if(p == q) Z(p,1) = Z(p,1) - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2

          end do
        end do
      end do

    end do
  end do

  !------------------!
  ! Spin-down sector !
  !------------------!

  do p=nC(2)+1,nBas-nR(2)
    do q=nC(2)+1,nBas-nR(2)

      ! Addition part: bb

      do i=nC(2)+1,nO(2)
        do a=nO(2)+1,nBas-nR(2)
          do b=nO(2)+1,nBas-nR(2)

            eps = eGF2(p,2) + eHF(i,2) - eHF(a,2) - eHF(b,2) 
            num = ERI_bb(i,q,a,b)*ERI_bb(a,b,i,p) &
                - ERI_bb(i,q,a,b)*ERI_bb(a,b,p,i)
         
            SigC(p,q,2) = SigC(p,q,2) + num*eps/(eps**2 + eta**2)
            if(p == q) Z(p,2) = Z(p,2) - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2

          end do
        end do
      end do

      ! Addition part: ab

      do i=nC(1)+1,nO(1)
        do a=nO(1)+1,nBas-nR(1)
          do b=nO(2)+1,nBas-nR(2)

            eps = eGF2(p,2) + eHF(i,1) - eHF(a,1) - eHF(b,2) 
            num = ERI_ab(i,q,a,b)*ERI_ab(a,b,i,p)
         
            SigC(p,q,2) = SigC(p,q,2) + num*eps/(eps**2 + eta**2)
            if(p == q) Z(p,2) = Z(p,2) - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2

          end do
        end do
      end do

     ! Removal part: bb

      do i=nC(2)+1,nO(2)
        do a=nO(2)+1,nBas-nR(2)
          do j=nC(2)+1,nO(2)

            eps = eGF2(p,2) + eHF(a,2) - eHF(i,2) - eHF(j,2) 
            num = ERI_bb(a,q,i,j)*ERI_bb(i,j,a,p) &
                - ERI_bb(a,q,i,j)*ERI_bb(i,j,p,a)

            SigC(p,q,2) = SigC(p,q,2) + num*eps/(eps**2 + eta**2)
            if(p == q) Z(p,2) = Z(p,2) - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2

          end do
        end do
      end do

     ! Removal part: ab

      do i=nC(1)+1,nO(1)
        do a=nO(1)+1,nBas-nR(1)
          do j=nC(2)+1,nO(2)

            eps = eGF2(p,2) + eHF(a,1) - eHF(i,1) - eHF(j,2) 
            num = ERI_ab(a,q,i,j)*ERI_ab(i,j,a,p)
         
            SigC(p,q,2) = SigC(p,q,2) + num*eps/(eps**2 + eta**2)
            if(p == q) Z(p,2) = Z(p,2) - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2

          end do
        end do
      end do

    end do
  end do

  Z(:,:) = 1d0/(1d0 - Z(:,:))

end subroutine 
