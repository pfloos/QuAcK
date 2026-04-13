double precision function GGF3_SigC(p,w,eta,nBas,nC,nO,nV,nR,e,ERI)

! Compute diagonal of the correlation part of the self-energy

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: p
  double precision,intent(in)   :: w
  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas,nC,nO,nV,nR
  double precision,intent(in)   :: e(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)


! Local variables

  integer                       :: i,j,k,l,a,b,c,d
  double precision              :: eps,eps1,eps2
  double precision              :: num
  double precision              :: SigInf
  double precision,allocatable  :: App(:),Cpp(:),Dpp(:)

! Output variable
  GGF3_SigC = 0d0

! Allocate
  
  allocate(App(6),Cpp(6),Dpp(6))

! Initialize 

  SigInf = 0d0

  App(:) = 0d0
  Cpp(:) = 0d0
  Dpp(:) = 0d0
  
!---------------------------!
! Second-order contribution !
!           2h1p            !
!---------------------------!
  do i=nC+1,nO
    do j=nC+1,nO
      do a=nO+1,nBas-nR

        eps = w + e(a) - e(i) - e(j)
        GGF3_SigC = GGF3_SigC + 0.5d0*(ERI(p,a,i,j) - ERI(p,a,j,i))**2*eps/(eps**2 + eta**2)

        end do
     end do
  end do
!---------------------------!
! Second-order contribution !
!           2p1h            !
!---------------------------!
  do i=nC+1,nO
    do a=nO+1,nBas-nR
      do b=nO+1,nBas-nR

        eps = w + e(i) - e(a) - e(b)
        GGF3_SigC = GGF3_SigC + 0.5d0*(ERI(p,i,a,b) - ERI(p,i,b,a))**2*eps/(eps**2 + eta**2)

      end do
    end do
  end do
! !--------------------------------------------------------!
! ! Compute third-order frequency-independent contribution !
! !                          3h2p                          !  
! !--------------------------------------------------------!
!   do i=nC+1,nO
!      do j=nC+1,nO
!         do k=nC+1,nO
!            do a=nO+1,nBas-nR
!               do b=nO+1,nBas-nR
                 
!                     eps1 = e(j) + e(i) - e(a) - e(b)
!                     eps2 = e(k) + e(i) - e(a) - e(b)
!                     num = - 0.5d0 * (ERI(p,k,p,j) - ERI(p,k,j,p)) * (ERI(j,i,a,b) - ERI(j,i,b,a)) * (ERI(a,b,k,i) - ERI(a,b,i,k))

!                     App(1) = App(1) + num / (eps1*eps2) 
                    
!                     eps1 = e(j) + e(i) - e(a) - e(b)
!                     eps2 = e(k)        - e(b)
!                     num  = - 0.5d0 * (ERI(p,b,p,k) - ERI(p,b,k,p)) * (ERI(j,i,a,b) - ERI(j,i,b,a)) * (ERI(i,j,k,a) - ERI(i,j,a,k))
                    
!                     App(5) = App(5) + num / (eps1*eps2)
                    
!                  end do
!               end do
!            end do
!         end do
!      end do
!   App(6) = App(5)
! !--------------------------------------------------------!
! ! Compute third-order frequency-independent contribution !
! !                          3p2h                          !  
! !--------------------------------------------------------!
!   do i=nC+1,nO
!      do j=nC+1,nO
!         do a=nO+1,nBas-nR
!            do b=nO+1,nBas-nR
!               do c=nO+1,nBas-nR
                 
!                  eps1 = e(j) + e(i) - e(a) - e(b)
!                  eps2 = e(j) + e(i) - e(a) - e(c)
!                  num = + 0.5d0 * (ERI(p,c,p,b) - ERI(p,c,b,p)) * (ERI(j,i,a,b) - ERI(j,i,b,a)) * (ERI(i,j,c,a) - ERI(i,j,a,c))
                 
!                  App(2) = App(2) + num / (eps1*eps2)
                 
!                  eps1 = e(j) + e(i) - e(a) - e(b)
!                  eps2 = e(j)        - e(c)
!                  num = + 0.5d0 * (ERI(p,c,p,j) - ERI(p,c,j,p)) * (ERI(j,i,a,b) - ERI(j,i,b,a)) * (ERI(a,b,c,i) - ERI(a,b,i,c))

!                  App(3) = App(3) + num / (eps1*eps2) 
                    
!               end do
!            end do
!         end do
!      end do
!   end do
!   App(4) = App(3)
!------------------------------------------------------!
! Compute third-order frequency-dependent contribution !
!                     4h1p C terms                     !  
!------------------------------------------------------!
  do i=nC+1,nO
     do j=nC+1,nO
        do k=nC+1,nO
           do l=nC+1,nO
              do a=nO+1,nBas-nR
                 
                 eps1 = w + e(a) - e(i) - e(j)
                 eps2 = w + e(a) - e(k) - e(l)
                 num = - 0.25d0 * (ERI(p,a,k,l) - ERI(p,a,l,k)) * (ERI(k,l,i,j) - ERI(k,l,j,i)) * (ERI(p,a,i,j) - ERI(p,a,j,i))
                 
                 Cpp(6) = Cpp(6) + num / (eps1*eps2)
                 
              end do
           end do
        end do
     end do
  end do
!------------------------------------------------------!
! Compute third-order frequency-dependent contribution !
!                     4p1h C terms                     !  
!------------------------------------------------------!
  do i=nC+1,nO
     do a=nO+1,nBas-nR
        do b=nO+1,nBas-nR
           do c=nO+1,nBas-nR
              do d=nO+1,nBas-nR
                 
                 eps1 = w + e(i) - e(a) - e(b)
                 eps2 = w + e(i) - e(c) - e(d)
                 num = + 0.25d0 * (ERI(p,i,a,b) - ERI(p,i,b,a)) * (ERI(a,b,c,d) - ERI(a,b,d,c)) * (ERI(p,i,c,d) - ERI(p,i,d,c))
                 
                 Cpp(1) = Cpp(1) + num / (eps1*eps2)
                 
              end do
           end do
        end do
     end do
  end do
!------------------------------------------------------!
! Compute third-order frequency-dependent contribution !
!                     3h2p C terms                     !  
!------------------------------------------------------!
  do i=nC+1,nO
     do j=nC+1,nO
        do k=nC+1,nO
           do a=nO+1,nBas-nR
              do b=nO+1,nBas-nR
                 
                 eps1 = w + e(i) - e(a) - e(b)
                 eps2 = e(j) + e(k) - e(a) - e(b)
                 num = + 0.25d0 * (ERI(p,i,a,b) - ERI(p,i,b,a)) * (ERI(a,b,j,k) - ERI(a,b,k,j)) * (ERI(p,i,j,k) - ERI(p,i,k,j))
                 
                 Cpp(2) = Cpp(2) + num / (eps1*eps2)
                 
              end do
           end do
        end do
     end do
  end do
  Cpp(3)  = Cpp(2)
!------------------------------------------------------!
! Compute third-order frequency-dependent contribution !
!                     3p2h C terms                     !  
!------------------------------------------------------!
  do i=nC+1,nO
     do j=nC+1,nO
        do a=nO+1,nBas-nR
           do b=nO+1,nBas-nR
              do c=nO+1,nBas-nR
                 
                 eps1 = w + e(a) - e(i) - e(j)
                 eps2 = e(i) + e(j) - e(b) - e(c)
                 num = + 0.25d0 * (ERI(p,a,i,j) - ERI(p,a,j,i)) * (ERI(i,j,b,c) - ERI(i,j,c,b)) * (ERI(p,a,b,c) - ERI(p,a,c,b))
                 
                 Cpp(4)  = Cpp(4) + num / (eps1*eps2)
                 
              end do
           end do
        end do
     end do
  end do
  Cpp(5)  = Cpp(4)
!------------------------------------------------------!
! Compute third-order frequency-dependent contribution !
!                     3h2p D terms                     !  
!------------------------------------------------------!
  do i=nC+1,nO
     do j=nC+1,nO
        do k=nC+1,nO
           do a=nO+1,nBas-nR
              do b=nO+1,nBas-nR
                 
                 eps1 = w + e(a) - e(i) - e(k)
                 eps2 = w + e(b) - e(j) - e(k)
                 num = - (ERI(p,a,k,i) - ERI(p,a,i,k)) * (ERI(i,b,a,j) - ERI(i,b,j,a)) * (ERI(p,b,k,j) - ERI(p,b,j,k))
                 
                 Dpp(6) = Dpp(6) + num / (eps1*eps2)
                 
                 
                 eps1 = w + e(a) - e(j) - e(k)
                 eps2 = e(i) + e(j) - e(a) - e(b)
                 num = + (ERI(p,a,k,j) - ERI(p,a,j,k)) * (ERI(j,i,a,b) - ERI(j,i,b,a))*(ERI(p,i,k,b) - ERI(p,i,b,k))
                  
                 Dpp(4) = Dpp(4) + num / (eps1*eps2)

                 
              end do
           end do
        end do
     end do
  end do
  Dpp(5) = Dpp(4)
!------------------------------------------------------!
! Compute third-order frequency-dependent contribution !
!                     3p2h D terms                     !  
!------------------------------------------------------!
  do i=nC+1,nO
     do j=nC+1,nO
        do a=nO+1,nBas-nR
           do b=nO+1,nBas-nR
              do c=nO+1,nBas-nR
                 
                 eps1 = w + e(i) - e(a) - e(b)
                 eps2 = w + e(j) - e(b) - e(c)
                 num = + (ERI(p,i,a,b) - ERI(p,i,b,a)) * (ERI(a,j,i,c) - ERI(a,j,c,i)) * (ERI(p,j,c,b) - ERI(p,j,b,c))
                  
                 Dpp(1) = Dpp(1) + num / (eps1*eps2)

                 
                 eps1 = w + e(i) - e(a) - e(c)
                 eps2 = e(i) + e(j) - e(a) - e(b)
                 num = + (ERI(p,i,c,a) - ERI(p,i,a,c)) * (ERI(a,b,i,j) - ERI(a,b,j,i)) * (ERI(p,b,c,j) - ERI(p,b,j,c))
                  
                 Dpp(2) = Dpp(2)  + num / (eps1*eps2)
                  
              end do
           end do
        end do
     end do
  end do
  Dpp(3) = Dpp(2)

!------------------------------!
! Collecting self-energy terms !  
!------------------------------!
  
  SigInf = App(1) + App(2) + App(3) + App(4) + App(5) + App(6)

  GGF3_SigC = GGF3_SigC + SigInf

  GGF3_SigC = GGF3_SigC + Cpp(1) + Cpp(2) + Cpp(3) + Cpp(4) + Cpp(5) + Cpp(6)
  
  GGF3_SigC = GGF3_SigC + Dpp(1) + Dpp(2) + Dpp(3) + Dpp(4) + Dpp(5) + Dpp(6)

  
end function
