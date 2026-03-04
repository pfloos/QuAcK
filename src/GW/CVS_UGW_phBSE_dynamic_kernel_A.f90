subroutine CVS_UGW_phBSE_dynamic_kernel_A(ispin,eta,nBas,nC,nO,nV,nR,nSa,nSb,nSt,nS_sc,lambda,eGW,  & 
                                      nCVS,nFC,occupations,virtuals,                            &
                                      ERI_aaaa,ERI_aabb,ERI_bbbb,OmRPA,rho_RPA,OmBSE,A_dyn,ZA_dyn)

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
  integer,intent(in)            :: nCVS(nspin)
  integer,intent(in)            :: nFC(nspin)
  integer,intent(in)            :: occupations(maxval(nO-nFC),nspin)
  integer,intent(in)            :: virtuals(nBas-minval(nO),nspin)
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
  double precision,intent(out)  :: ZA_dyn(nSt,nSt)

!--------------------------------------------------!
! Build BSE matrix for spin-conserving transitions !
!--------------------------------------------------!

  A_dyn(:,:) = 0d0

  if(ispin == 1) then

    ! aaaa block
 
    ia = 0
    do i=1,nO(1)-nFC(1)
      do a=nCVS(1)+1,nBas-nO(1)
        ia = ia + 1
        jb = 0
        do j=1,nO(1)-nFC(1)
          do b=nCVS(1)+1,nBas-nO(1)
            jb = jb + 1
  
            chi = 0d0
            do kc=1,nS_sc
              chi = chi + rho_RPA(occupations(i,1),occupations(j,1),kc,1)&
                        * rho_RPA(virtuals(a,1),virtuals(b,1),kc,1)      &
                        * OmRPA(kc)/(OmRPA(kc)**2 + eta**2)
            end do
 
            A_dyn(ia,jb) = A_dyn(ia,jb) - 2d0*lambda*chi

            chi = 0d0
            do kc=1,nS_sc
           
              eps = + OmBSE - OmRPA(kc) - (eGW(virtuals(a,1),1) - eGW(occupations(j,1),1))
              chi = chi + rho_RPA(occupations(i,1),occupations(j,1),kc,1)*rho_RPA(virtuals(a,1),virtuals(b,1),kc,1)*eps/(eps**2 + eta**2)
           
              eps = + OmBSE - OmRPA(kc) - (eGW(virtuals(b,1),1) - eGW(occupations(i,1),1))
              chi = chi + rho_RPA(occupations(i,1),occupations(j,1),kc,1)*rho_RPA(virtuals(a,1),virtuals(b,1),kc,1)*eps/(eps**2 + eta**2)
           
            end do

            A_dyn(ia,jb) = A_dyn(ia,jb) - lambda*chi

            chi = 0d0
            do kc=1,nS_sc
  
              eps = + OmBSE - OmRPA(kc) - (eGW(virtuals(a,1),1) - eGW(occupations(j,1),1))
              chi = chi + rho_RPA(occupations(i,1),occupations(j,1),kc,1)*rho_RPA(virtuals(a,1),virtuals(b,1),kc,1)*(eps**2 - eta**2)/(eps**2 + eta**2)**2
  
              eps = + OmBSE - OmRPA(kc) - (eGW(virtuals(b,1),1) - eGW(occupations(i,1),1))
              chi = chi + rho_RPA(i,j,kc,1)*rho_RPA(a,b,kc,1)*(eps**2 - eta**2)/(eps**2 + eta**2)**2
  
            end do

            ZA_dyn(ia,jb) = ZA_dyn(ia,jb) + lambda*chi
 
          end do
        end do
      end do
    end do

    ! bbbb block
 
    ia = 0
    do i=1,nO(2)-nFC(2)
      do a=nCVS(2)+1,nBas-nO(2)
        ia = ia + 1
        jb = 0
        do j=1,nO(2)-nFC(2)
          do b=nCVS(2)+1,nBas-nO(2)
            jb = jb + 1
  
            chi = 0d0
            do kc=1,nS_sc
              chi = chi + rho_RPA(occupations(i,2),occupations(j,2),kc,2)&
                        * rho_RPA(virtuals(a,2),virtuals(b,2),kc,2)      &
                        * OmRPA(kc)/(OmRPA(kc)**2 + eta**2)
            end do
 
            A_dyn(nSa+ia,nSa+jb) = A_dyn(nSa+ia,nSa+jb) - 2d0*lambda*chi

            chi = 0d0
            do kc=1,nS_sc

              eps = + OmBSE - OmRPA(kc) - (eGW(virtuals(a,2),2) - eGW(occupations(j,2),2))
              chi = chi + rho_RPA(occupations(i,2),occupations(j,2),kc,2)*rho_RPA(virtuals(a,2),virtuals(b,2),kc,2)*eps/(eps**2 + eta**2)

              eps = + OmBSE - OmRPA(kc) - (eGW(virtuals(b,2),2) - eGW(occupations(i,2),2))
              chi = chi + rho_RPA(occupations(i,2),occupations(j,2),kc,2)*rho_RPA(virtuals(a,2),virtuals(b,2),kc,2)*eps/(eps**2 + eta**2)

            end do

            A_dyn(nSa+ia,nSa+jb) = A_dyn(nSa+ia,nSa+jb) - lambda*chi

            chi = 0d0
            do kc=1,nS_sc

              eps = + OmBSE - OmRPA(kc) - (eGW(virtuals(a,2),2) - eGW(occupations(j,2),2))
              chi = chi + rho_RPA(occupations(i,2),occupations(j,2),kc,2)*rho_RPA(virtuals(a,2),virtuals(b,2),kc,2)*(eps**2 - eta**2)/(eps**2 + eta**2)**2

              eps = + OmBSE - OmRPA(kc) - (eGW(virtuals(b,2),2) - eGW(occupations(i,2),2))
              chi = chi + rho_RPA(occupations(i,2),occupations(j,2),kc,2)*rho_RPA(virtuals(a,2),virtuals(b,2),kc,2)*(eps**2 - eta**2)/(eps**2 + eta**2)**2

            end do

            ZA_dyn(nSa+ia,nSa+jb) = ZA_dyn(nSa+ia,nSa+jb) + lambda*chi
 
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
    do i=1,nO(1)-nFC(1)
      do a=nCVS(2)+1,nBas-nO(2)
        ia = ia + 1
        jb = 0
        do j=1,nO(1)-nFC(1)
          do b=nCVS(2)+1,nBas-nO(2)
            jb = jb + 1

            chi = 0d0
            do kc=1,nS_sc
              chi = chi + rho_RPA(occupations(i,1),occupations(j,1),kc,1)*rho_RPA(virtuals(a,2),virtuals(b,2),kc,2)*OmRPA(kc)/(OmRPA(kc)**2 + eta**2)
            end do

            A_dyn(ia,jb) = A_dyn(ia,jb) - 2d0*lambda*chi

            chi = 0d0
            do kc=1,nS_sc

              eps = + OmBSE - OmRPA(kc) - (eGW(virtuals(a,2),2) - eGW(occupations(j,1),1))
              chi = chi + rho_RPA(occupations(i,1),occupations(j,1),kc,1)*rho_RPA(virtuals(a,2),virtuals(b,2),kc,2)*eps/(eps**2 + eta**2)

              eps = + OmBSE - OmRPA(kc) - (eGW(virtuals(b,2),2) - eGW(occupations(i,1),1))
              chi = chi + rho_RPA(occupations(i,1),occupations(j,1),kc,1)*rho_RPA(virtuals(a,2),virtuals(b,2),kc,2)*eps/(eps**2 + eta**2)

            end do

            A_dyn(ia,jb) = A_dyn(ia,jb) - lambda*chi

            chi = 0d0
            do kc=1,nS_sc

              eps = + OmBSE - OmRPA(kc) - (eGW(virtuals(a,2),2) - eGW(occupations(j,1),1))
              chi = chi + rho_RPA(occupations(i,1),occupations(j,1),kc,1)*rho_RPA(virtuals(a,2),virtuals(b,2),kc,2)*(eps**2 - eta**2)/(eps**2 + eta**2)**2

              eps = + OmBSE - OmRPA(kc) - (eGW(virtuals(b,2),2) - eGW(occupations(i,1),1))
              chi = chi + rho_RPA(occupations(i,1),occupations(j,1),kc,1)*rho_RPA(virtuals(a,2),virtuals(b,2),kc,2)*(eps**2 - eta**2)/(eps**2 + eta**2)**2

            end do

            ZA_dyn(ia,jb) = ZA_dyn(ia,jb) + lambda*chi

          end  do
        end  do
      end  do
    end  do

    ! baba block

    ia = 0
    do i=1,nO(2)-nFC(2)
      do a=nCVS(1)+1,nBas-nO(1)
        ia = ia + 1
        jb = 0
        do j=1,nO(2)-nFC(2)
          do b=nCVS(1)+1,nBas-nO(1)
            jb = jb + 1

            chi = 0d0
            do kc=1,nS_sc
              chi = chi + rho_RPA(occupations(i,2),occupations(j,2),kc,2)*rho_RPA(virtuals(a,1),virtuals(b,1),kc,1)*OmRPA(kc)/(OmRPA(kc)**2 + eta**2)
            end do

            A_dyn(nSa+ia,nSa+jb) = A_dyn(nSa+ia,nSa+jb) - 2d0*lambda*chi

            chi = 0d0
            do kc=1,nS_sc

              eps = + OmBSE - OmRPA(kc) - (eGW(virtuals(a,1),1) - eGW(occupations(j,2),2))
              chi = chi + rho_RPA(occupations(i,2),occupations(j,2),kc,2)*rho_RPA(virtuals(a,1),virtuals(b,1),kc,1)*eps/(eps**2 + eta**2)

              eps = + OmBSE - OmRPA(kc) - (eGW(virtuals(b,1),1) - eGW(occupations(i,2),2))
              chi = chi + rho_RPA(occupations(i,2),occupations(j,2),kc,2)*rho_RPA(virtuals(a,1),virtuals(b,1),kc,1)*eps/(eps**2 + eta**2)

            end do

            A_dyn(nSa+ia,nSa+jb) = A_dyn(nSa+ia,nSa+jb) - lambda*chi

            chi = 0d0
            do kc=1,nS_sc

              eps = + OmBSE - OmRPA(kc) - (eGW(virtuals(a,1),1) - eGW(occupations(j,2),2))
              chi = chi + rho_RPA(occupations(i,2),occupations(j,2),kc,2)*rho_RPA(virtuals(a,1),virtuals(b,1),kc,1)*(eps**2 - eta**2)/(eps**2 + eta**2)**2

              eps = + OmBSE - OmRPA(kc) - (eGW(virtuals(b,1),1) - eGW(occupations(i,2),2))
              chi = chi + rho_RPA(occupations(i,2),occupations(j,2),kc,2)*rho_RPA(virtuals(a,1),virtuals(b,1),kc,1)*(eps**2 - eta**2)/(eps**2 + eta**2)**2

            end do

            ZA_dyn(nSa+ia,nSa+jb) = ZA_dyn(nSa+ia,nSa+jb) + lambda*chi

          end  do
        end  do
      end  do
    end  do

  end if

end subroutine 
