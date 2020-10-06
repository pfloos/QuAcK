subroutine unrestricted_Bethe_Salpeter_A_matrix_dynamic(ispin,eta,nBas,nC,nO,nV,nR,nSa,nSb,nSt,nS_sc,lambda,eGW, & 
                                                        ERI_aaaa,ERI_aabb,ERI_bbbb,OmRPA,rho_RPA,OmBSE,A_dyn)

! Compute the extra term for dynamical Bethe-Salpeter equation for linear response in the unrestricted formalism

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
  double precision,intent(in)   :: eGW(nBas,nspin)
  double precision,intent(in)   :: ERI_aaaa(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_aabb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_bbbb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: OmRPA(nS_sc)
  double precision,intent(in)   :: rho_RPA(nBas,nBas,nS_sc,nspin)
  double precision,intent(in)   :: OmBSE
  
! Local variables

  double precision              :: chi
  double precision              :: eps
  integer                       :: i,j,a,b,ia,jb,kc

! Output variables

  double precision,intent(out)  :: A_dyn(nSt,nSt)

!--------------------------------------------------!
! Build BSE matrix for spin-conserving transitions !
!--------------------------------------------------!

  A_dyn(:,:) = 0d0

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
              chi = chi + rho_RPA(i,j,kc,1)*rho_RPA(a,b,kc,1)*OmRPA(kc)/(OmRPA(kc)**2 + eta**2)
            enddo
 
            A_dyn(ia,jb) = A_dyn(ia,jb) - 2d0*lambda*chi

            chi = 0d0
            do kc=1,nS_sc
           
              eps = + OmBSE - OmRPA(kc) - (eGW(a,1) - eGW(j,1))
              chi = chi + rho_RPA(i,j,kc,1)*rho_RPA(a,b,kc,1)*eps/(eps**2 + eta**2)
           
              eps = + OmBSE - OmRPA(kc) - (eGW(b,1) - eGW(i,1))
              chi = chi + rho_RPA(i,j,kc,1)*rho_RPA(a,b,kc,1)*eps/(eps**2 + eta**2)
           
            enddo

            A_dyn(ia,jb) = A_dyn(ia,jb) - lambda*chi
 
          enddo
        enddo
      enddo
    enddo

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
              chi = chi + rho_RPA(i,j,kc,2)*rho_RPA(a,b,kc,2)*OmRPA(kc)/(OmRPA(kc)**2 + eta**2)
            enddo
 
            A_dyn(nSa+ia,nSa+jb) = A_dyn(nSa+ia,nSa+jb) - 2d0*lambda*chi

            chi = 0d0
            do kc=1,nS_sc

              eps = + OmBSE - OmRPA(kc) - (eGW(a,2) - eGW(j,2))
              chi = chi + rho_RPA(i,j,kc,2)*rho_RPA(a,b,kc,2)*eps/(eps**2 + eta**2)

              eps = + OmBSE - OmRPA(kc) - (eGW(b,2) - eGW(i,2))
              chi = chi + rho_RPA(i,j,kc,2)*rho_RPA(a,b,kc,2)*eps/(eps**2 + eta**2)

            enddo

            A_dyn(nSa+ia,nSa+jb) = A_dyn(nSa+ia,nSa+jb) - lambda*chi
 
          enddo
        enddo
      enddo
    enddo

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
              chi = chi + rho_RPA(i,j,kc,1)*rho_RPA(a,b,kc,2)*OmRPA(kc)/(OmRPA(kc)**2 + eta**2)
            enddo

            A_dyn(ia,jb) = A_dyn(ia,jb) - 2d0*lambda*chi

            chi = 0d0
            do kc=1,nS_sc

              eps = + OmBSE - OmRPA(kc) - (eGW(a,2) - eGW(j,1))
              chi = chi + rho_RPA(i,j,kc,1)*rho_RPA(a,b,kc,2)*eps/(eps**2 + eta**2)

              eps = + OmBSE - OmRPA(kc) - (eGW(b,2) - eGW(i,1))
              chi = chi + rho_RPA(i,j,kc,1)*rho_RPA(a,b,kc,2)*eps/(eps**2 + eta**2)

            enddo

            A_dyn(ia,jb) = A_dyn(ia,jb) - lambda*chi

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
              chi = chi + rho_RPA(i,j,kc,2)*rho_RPA(a,b,kc,1)*OmRPA(kc)/(OmRPA(kc)**2 + eta**2)
            enddo

            A_dyn(nSa+ia,nSa+jb) = A_dyn(nSa+ia,nSa+jb) - 2d0*lambda*chi

            chi = 0d0
            do kc=1,nS_sc

              eps = + OmBSE - OmRPA(kc) - (eGW(a,1) - eGW(j,2))
              chi = chi + rho_RPA(i,j,kc,2)*rho_RPA(a,b,kc,1)*eps/(eps**2 + eta**2)

              eps = + OmBSE - OmRPA(kc) - (eGW(b,1) - eGW(i,2))
              chi = chi + rho_RPA(i,j,kc,2)*rho_RPA(a,b,kc,1)*eps/(eps**2 + eta**2)

            enddo

            A_dyn(nSa+ia,nSa+jb) = A_dyn(nSa+ia,nSa+jb) - lambda*chi

          end  do
        end  do
      end  do
    end  do

  end if

end subroutine unrestricted_Bethe_Salpeter_A_matrix_dynamic
