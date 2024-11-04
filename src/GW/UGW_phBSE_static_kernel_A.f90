subroutine UGW_phBSE_static_kernel_A(ispin,eta,nBas,nC,nO,nV,nR,nSa,nSb,nSt,nS_sc,lambda,Om,rho,KA)

! Compute the extra term for Bethe-Salpeter equation for linear response in the unrestricted formalism

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: ispin
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC(nspin)
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nV(nspin)
  integer,intent(in)            :: nR(nspin)
  integer,intent(in)            :: nSa
  integer,intent(in)            :: nSb
  integer,intent(in)            :: nSt
  integer,intent(in)            :: nS_sc
  double precision,intent(in)   :: eta
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: Om(nS_sc)
  double precision,intent(in)   :: rho(nBas,nBas,nS_sc,nspin)
  
! Local variables

  double precision              :: chi
  double precision              :: eps
  integer                       :: i,j,a,b,ia,jb,kc

! Output variables

  double precision,intent(out)  :: KA(nSt,nSt)

! Initialization

  KA(:,:) = 0d0

!--------------------------------------------------!
! Build BSE matrix for spin-conserving transitions !
!--------------------------------------------------!

  if(ispin == 1) then

    ! aaaa block
 
    ia = 0
    do i=nC(1)+1,nO(1)
      do a=nO(1)+1,nBas-nR(1)
        ia = ia + 1
        jb = 0
        do j=nC(1)+1,nO(1)
          do b=nO(1)+1,nBas-nR(1)
            jb = jb + 1
  
            chi = 0d0
            do kc=1,nS_sc
              eps = Om(kc)**2 + eta**2
              chi = chi + rho(i,j,kc,1)*rho(a,b,kc,1)*Om(kc)/eps
            end do
 
            KA(ia,jb) = 2d0*lambda**2*chi
 
          end do
        end do
      end do
    end do

    ! bbbb block
 
    ia = 0
    do i=nC(2)+1,nO(2)
      do a=nO(2)+1,nBas-nR(2)
        ia = ia + 1
        jb = 0
        do j=nC(2)+1,nO(2)
          do b=nO(2)+1,nBas-nR(2)
            jb = jb + 1
  
            chi = 0d0
            do kc=1,nS_sc
              eps = Om(kc)**2 + eta**2
              chi = chi + rho(i,j,kc,2)*rho(a,b,kc,2)*Om(kc)/eps 
            end do
 
            KA(nSa+ia,nSa+jb) = 2d0*lambda**2*chi
 
          end do
        end do
      end do
    end do

  end if

!--------------------------------------------!
! Build BSE matrix for spin-flip transitions !
!--------------------------------------------!

  if(ispin == 2) then

    ! abab block

    ia = 0
    do i=nC(1)+1,nO(1)
      do a=nO(2)+1,nBas-nR(2)
        ia = ia + 1
        jb = 0
        do j=nC(1)+1,nO(1)
          do b=nO(2)+1,nBas-nR(2)
            jb = jb + 1

            chi = 0d0
            do kc=1,nS_sc
              eps = Om(kc)**2 + eta**2
              chi = chi + rho(i,j,kc,1)*rho(a,b,kc,2)*Om(kc)/eps
            end do

            KA(ia,jb) = 2d0*lambda**2*chi

          end  do
        end  do
      end  do
    end  do

    ! baba block

    ia = 0
    do i=nC(2)+1,nO(2)
      do a=nO(1)+1,nBas-nR(1)
        ia = ia + 1
        jb = 0
        do j=nC(2)+1,nO(2)
          do b=nO(1)+1,nBas-nR(1)
            jb = jb + 1

            chi = 0d0
            do kc=1,nS_sc
              eps = Om(kc)**2 + eta**2
              chi = chi + rho(i,j,kc,2)*rho(a,b,kc,1)*Om(kc)/eps
            end do

            KA(nSa+ia,nSa+jb) = 2d0*lambda**2*chi

          end  do
        end  do
      end  do
    end  do

  end if

end subroutine 
