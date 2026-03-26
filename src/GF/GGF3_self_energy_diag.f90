subroutine GGF3_self_energy_diag(eta,nBas,nC,nO,nV,nR,e,ERI,SigC,Z)

! Compute diagonal part of the GF3 self-energy and its renormalization factor

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  double precision,intent(in)   :: e(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: i,j,k,l,a,b,c,d
  integer                       :: p,q
  double precision              :: eps,eps1,eps2
  double precision              :: num
  double precision,allocatable  :: SigInf(:)
  double precision,allocatable  :: App(:,:),Cpp(:,:),Dpp(:,:)
  double precision,allocatable  :: ZCpp(:,:),ZDpp(:,:)

! Output variables

  double precision,intent(out)  :: SigC(nBas)
  double precision,intent(out)  :: Z(nBas)

! Allocate

  allocate(SigInf(nBas),App(nBas,6),Cpp(nBas,6),Dpp(nBas,6))
  allocate(ZCpp(nBas,6),ZDpp(nBas,6))
  
! Initialize 

  SigInf(:) = 0d0

  App(:,:) = 0d0
  Cpp(:,:) = 0d0
  Dpp(:,:) = 0d0

  ZCpp(:,:) = 0d0
  ZDpp(:,:) = 0d0
  
  SigC(:) = 0d0
  Z(:)    = 0d0


  
!---------------------------!
! Second-order contribution !
!           2h1p            !
!---------------------------!
  do p=nC+1,nBas-nR
    do i=nC+1,nO
      do j=nC+1,nO
        do a=nO+1,nBas-nR

          eps = e(p) + e(a) - e(i) - e(j)
          num = 0.5d0*(ERI(p,a,i,j) - ERI(p,a,j,i))**2

          SigC(p) = SigC(p) + num*eps/(eps**2 + eta**2)
          Z(p)    = Z(p)   - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2

        end do
      end do
    end do
  end do
!---------------------------!
! Second-order contribution !
!           2p1h            !
!---------------------------!
  do p=nC+1,nBas-nR
    do i=nC+1,nO
      do a=nO+1,nBas-nR
        do b=nO+1,nBas-nR

          eps = e(p) + e(i) - e(a) - e(b)
          num = 0.5d0*(ERI(p,i,a,b) - ERI(p,i,b,a))**2

          SigC(p) = SigC(p) + num*eps/(eps**2 + eta**2)
          Z(p)    = Z(p)   - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2

        end do
      end do
    end do
  end do
!--------------------------------------------------------!
! Compute third-order frequency-independent contribution !
!                          3h2p                          !  
!--------------------------------------------------------!
  do p=nC+1,nBas-nR
     do i=nC+1,nO
        do j=nC+1,nO
           do k=nC+1,nO
              do a=nO+1,nBas-nR
                 do b=nO+1,nBas-nR
                    
                    eps1 = e(j) + e(i) - e(a) - e(b)
                    eps2 = e(k) + e(i) - e(a) - e(b)
                    num = - 0.5d0 * (ERI(p,k,p,j) - ERI(p,k,j,p)) * (ERI(j,i,a,b) - ERI(j,i,b,a)) * (ERI(a,b,k,i) - ERI(a,b,i,k))

                    App(p,1) = App(p,1) + num / (eps1*eps2) 
                    
                    eps1 = e(j) + e(i) - e(a) - e(b)
                    eps2 = e(k)        - e(b)
                    num  = - 0.5d0 * (ERI(p,b,p,k) - ERI(p,b,k,p)) * (ERI(j,i,a,b) - ERI(j,i,b,a)) * (ERI(i,j,k,a) - ERI(i,j,a,k))
                    
                    App(p,5) = App(p,5) + num / (eps1*eps2)
                    
                 end do
              end do
           end do
        end do
     end do
  end do
  App(:,6) = App(:,5)
!--------------------------------------------------------!
! Compute third-order frequency-independent contribution !
!                          3p2h                          !  
!--------------------------------------------------------!
  do p=nC+1,nBas-nR
     do i=nC+1,nO
        do j=nC+1,nO
           do a=nO+1,nBas-nR
              do b=nO+1,nBas-nR
                 do c=nO+1,nBas-nR

                    eps1 = e(j) + e(i) - e(a) - e(b)
                    eps2 = e(j) + e(i) - e(a) - e(c)
                    num = + 0.5d0 * (ERI(p,c,p,b) - ERI(p,c,b,p)) * (ERI(j,i,a,b) - ERI(j,i,b,a)) * (ERI(i,j,c,a) - ERI(i,j,a,c))

                    App(p,2) = App(p,2) + num / (eps1*eps2)
                  
                    eps1 = e(j) + e(i) - e(a) - e(b)
                    eps2 = e(j)        - e(c)
                    num = + 0.5d0 * (ERI(p,c,p,j) - ERI(p,c,j,p)) * (ERI(j,i,a,b) - ERI(j,i,b,a)) * (ERI(a,b,c,i) - ERI(a,b,i,c))

                    App(p,3) = App(p,3) + num / (eps1*eps2) 
                  
                 end do
              end do
           end do
        end do
     end do
  end do
  App(:,4) = App(:,3)
!------------------------------------------------------!
! Compute third-order frequency-dependent contribution !
!                     4h1p C terms                     !  
!------------------------------------------------------!
  do p=nC+1,nBas-nR
     do i=nC+1,nO
        do j=nC+1,nO
           do k=nC+1,nO
              do l=nC+1,nO
                 do a=nO+1,nBas-nR
                  
                    eps1 = e(p) + e(a) - e(i) - e(j)
                    eps2 = e(p) + e(a) - e(k) - e(l)
                    num = - 0.25d0 * (ERI(p,a,k,l) - ERI(p,a,l,k)) * (ERI(k,l,i,j) - ERI(k,l,j,i)) * (ERI(p,a,i,j) - ERI(p,a,j,i))
                  
                    Cpp(p,6)  = Cpp(p,6)  + num / (eps1*eps2)
                    ZCpp(p,6) = ZCpp(p,6) - num*(eps1**2 - eta**2)/(eps1**2 + eta**2)**2/eps2
                    ZCpp(p,6) = ZCpp(p,6) - num*(eps2**2 - eta**2)/(eps2**2 + eta**2)**2/eps1
                  
                 end do
              end do
           end do
        end do
     end do
  end do
!------------------------------------------------------!
! Compute third-order frequency-dependent contribution !
!                     4p1h C terms                     !  
!------------------------------------------------------!
  do p=nC+1,nBas-nR
    do i=nC+1,nO
       do a=nO+1,nBas-nR
          do b=nO+1,nBas-nR
             do c=nO+1,nBas-nR
                do d=nO+1,nBas-nR

                   eps1 = e(p) + e(i) - e(a) - e(b)
                   eps2 = e(p) + e(i) - e(c) - e(d)
                   num = + 0.25d0 * (ERI(p,i,a,b) - ERI(p,i,b,a)) * (ERI(a,b,c,d) - ERI(a,b,d,c)) * (ERI(p,i,c,d) - ERI(p,i,d,c))

                   Cpp(p,1)  = Cpp(p,1)  + num / (eps1*eps2)
                   ZCpp(p,1) = ZCpp(p,1) - num*(eps1**2 - eta**2)/(eps1**2 + eta**2)**2/eps2
                   ZCpp(p,1) = ZCpp(p,1) - num*(eps2**2 - eta**2)/(eps2**2 + eta**2)**2/eps1
                 
                end do
             end do
          end do
       end do
    end do
  end do
!------------------------------------------------------!
! Compute third-order frequency-dependent contribution !
!                     3h2p C terms                     !  
!------------------------------------------------------!
  do p=nC+1,nBas-nR
     do i=nC+1,nO
        do j=nC+1,nO
           do k=nC+1,nO
              do a=nO+1,nBas-nR
                 do b=nO+1,nBas-nR
                  
                    eps1 = e(p) + e(i) - e(a) - e(b)
                    eps2 = e(j) + e(k) - e(a) - e(b)
                    num = + 0.25d0 * (ERI(p,i,a,b) - ERI(p,i,b,a)) * (ERI(a,b,j,k) - ERI(a,b,k,j)) * (ERI(p,i,j,k) - ERI(p,i,k,j))
                  
                    Cpp(p,2)  = Cpp(p,2)  + num / (eps1*eps2)
                    ZCpp(p,2) = ZCpp(p,2) - num*(eps1**2 - eta**2)/(eps1**2 + eta**2)**2/eps2
                  
                 end do
              end do
           end do
        end do
     end do
  end do
  Cpp(:,3)  = Cpp(:,2)
  ZCpp(:,3) = ZCpp(:,2)
!------------------------------------------------------!
! Compute third-order frequency-dependent contribution !
!                     3p2h C terms                     !  
!------------------------------------------------------!
  do p=nC+1,nBas-nR
    do i=nC+1,nO
       do j=nC+1,nO
          do a=nO+1,nBas-nR
             do b=nO+1,nBas-nR
                do c=nO+1,nBas-nR
                 
                   eps1 = e(p) + e(a) - e(i) - e(j)
                   eps2 = e(i) + e(j) - e(b) - e(c)
                   num = + 0.25d0 * (ERI(p,a,i,j) - ERI(p,a,j,i)) * (ERI(i,j,b,c) - ERI(i,j,c,b)) * (ERI(p,a,b,c) - ERI(p,a,c,b))
                 
                   Cpp(p,4)  = Cpp(p,4)  + num / (eps1*eps2)
                   ZCpp(p,4) = ZCpp(p,4) - num*(eps1**2 - eta**2)/(eps1**2 + eta**2)**2/eps2
                 
                end do
             end do
          end do
       end do
    end do
  end do
  Cpp(:,5)  = Cpp(:,4)
  ZCpp(:,5) = ZCpp(:,4)
!------------------------------!
! Collecting self-energy terms !  
!------------------------------!
  
  SigInf(:) = App(:,1) + App(:,2) + App(:,3) + App(:,4) + App(:,5) + App(:,6)

  SigC(:) = SigC(:) + SigInf(:)

  SigC(:) = SigC(:) + Cpp(:,1) + Cpp(:,2) + Cpp(:,3) + Cpp(:,4) + Cpp(:,5) + Cpp(:,6)
  
  SigC(:) = SigC(:) + Dpp(:,1) + Dpp(:,2) + Dpp(:,3) + Dpp(:,4) + Dpp(:,5) + Dpp(:,6)

!----------------------------------!
! Collecting renormalization terms !  
!----------------------------------!
  
  Z(:) = Z(:) + ZCpp(:,1) + ZCpp(:,2) + ZCpp(:,3) + ZCpp(:,4) + ZCpp(:,5) + ZCpp(:,6)
  
  Z(:) = Z(:) + ZDpp(:,1) + ZDpp(:,2) + ZDpp(:,3) + ZDpp(:,4) + ZDpp(:,5) + ZDpp(:,6)
  
  Z(:) = 1d0/(1d0 - Z(:))

end subroutine 
