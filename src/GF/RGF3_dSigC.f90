double precision function RGF3_dSigC(p,w,eta,nBas,nC,nO,nV,nR,e,ERI)

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
  double precision,allocatable  :: ZCpp(:),ZDpp(:)

! Output variable
  RGF3_dSigC = 0d0

! Allocate
  
  allocate(ZCpp(6),ZDpp(6))

! Initialize 

  SigInf = 0d0

  ZCpp(:) = 0d0
  ZDpp(:) = 0d0

!---------------------------!
! Second-order contribution !
!           2h1p            !
!---------------------------!
  do i=nC+1,nO
     do j=nC+1,nO
        do a=nO+1,nBas-nR
           
           eps = w + e(a) - e(i) - e(j)
           num = (2d0*ERI(p,a,i,j) - ERI(p,a,j,i))*ERI(p,a,i,j)
           RGF3_dSigC = RGF3_dSigC - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2
           
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
           num = (2d0*ERI(p,i,a,b) - ERI(p,i,b,a))*ERI(p,i,a,b)
           RGF3_dSigC = RGF3_dSigC - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2
           
        end do
     end do
  end do
!------------------------------------------------------!
! There is not frequency-independent contribution in Z !  
!------------------------------------------------------!
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
                 num = - (2d0*ERI(p,a,k,l) - ERI(p,a,l,k)) * ERI(k,l,i,j) * ERI(p,a,i,j)
                 
                 ZCpp(6) = ZCpp(6) - num*(eps1**2 - eta**2)/(eps1**2 + eta**2)**2/eps2
                 ZCpp(6) = ZCpp(6) - num*(eps2**2 - eta**2)/(eps2**2 + eta**2)**2/eps1
                 
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
                 num = + (2d0*ERI(p,i,a,b) - ERI(p,i,b,a)) * ERI(a,b,c,d) * ERI(p,i,c,d)
                 
                 ZCpp(1) = ZCpp(1) - num*(eps1**2 - eta**2)/(eps1**2 + eta**2)**2/eps2
                 ZCpp(1) = ZCpp(1) - num*(eps2**2 - eta**2)/(eps2**2 + eta**2)**2/eps1
                 
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
                 num = + (2d0*ERI(p,i,a,b) - ERI(p,i,b,a)) * ERI(a,b,j,k) * ERI(p,i,j,k)
                 
                 ZCpp(2) = ZCpp(2) - num*(eps1**2 - eta**2)/(eps1**2 + eta**2)**2/eps2
                 
              end do
           end do
        end do
     end do
  end do
  ZCpp(3) = ZCpp(2)
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
                 num = + (2d0*ERI(p,a,i,j) - ERI(p,a,j,i)) * ERI(i,j,b,c) * ERI(p,a,b,c)
                 
                 ZCpp(4) = ZCpp(4) - num*(eps1**2 - eta**2)/(eps1**2 + eta**2)**2/eps2
                 
              end do
           end do
        end do
     end do
  end do
  ZCpp(5) = ZCpp(4)
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
                 num = - ERI(p,a,k,i) * ( ERI(i,b,a,j)*(4d0*ERI(p,b,k,j) - 2d0*ERI(p,b,j,k)) &
                                        + ERI(i,b,j,a) * (    ERI(p,b,j,k) - 2d0*ERI(p,b,k,j)) )
                 
                 ZDpp(6) = ZDpp(6) - num*(eps1**2 - eta**2)/(eps1**2 + eta**2)**2/eps2
                 ZDpp(6) = ZDpp(6) - num*(eps2**2 - eta**2)/(eps2**2 + eta**2)**2/eps1
                 
                 num = - ERI(p,a,i,k) * ( ERI(i,b,a,j)*(    ERI(p,b,j,k) - 2d0*ERI(p,b,k,j)) &
                                        + ERI(i,b,j,a)*(    ERI(p,b,k,j) - 2d0*ERI(p,b,j,k)) )

                 ZDpp(6) = ZDpp(6) - num*(eps1**2 - eta**2)/(eps1**2 + eta**2)**2/eps2
                 ZDpp(6) = ZDpp(6) - num*(eps2**2 - eta**2)/(eps2**2 + eta**2)**2/eps1
                 
                 eps1 = w + e(a) - e(j) - e(k)
                 eps2 = e(i) + e(j) - e(a) - e(b)
                 num = + ERI(p,a,k,j) * ( ERI(j,i,a,b)*(4d0*ERI(p,i,k,b) - 2d0*ERI(p,i,b,k)) &
                                        + ERI(j,i,b,a)*(    ERI(p,i,b,k) - 2d0*ERI(p,i,k,b)) )
                  
                 ZDpp(4) = ZDpp(4) - num*(eps1**2 - eta**2)/(eps1**2 + eta**2)**2/eps2

                 num = + ERI(p,a,j,k) * ( ERI(j,i,a,b)*(    ERI(p,i,b,k) - 2d0*ERI(p,i,k,b)) &
                                        + ERI(j,i,b,a)*(    ERI(p,i,k,b) - 2d0*ERI(p,i,b,k)) )
                  
                 ZDpp(4) = ZDpp(4) - num*(eps1**2 - eta**2)/(eps1**2 + eta**2)**2/eps2
                 
              end do
           end do
        end do
     end do
  end do
  ZDpp(5) = ZDpp(4)
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
                 num = + ERI(p,i,a,b) * ( ERI(a,j,i,c)*(    ERI(p,j,c,b) - 2d0*ERI(p,j,b,c)) &
                                        + ERI(a,j,c,i)*(    ERI(p,j,b,c) - 2d0*ERI(p,j,c,b)) )
                  
                 ZDpp(1) = ZDpp(1) - num*(eps1**2 - eta**2)/(eps1**2 + eta**2)**2/eps2
                 ZDpp(1) = ZDpp(1) - num*(eps2**2 - eta**2)/(eps2**2 + eta**2)**2/eps1
                  
                 num = + ERI(p,i,b,a) * ( ERI(a,j,i,c)*(4d0*ERI(p,j,b,c) - 2d0*ERI(p,j,c,b)) &
                                        + ERI(a,j,c,i)*(    ERI(p,j,c,b) - 2d0*ERI(p,j,b,c)) )
                  
                 ZDpp(1) = ZDpp(1) - num*(eps1**2 - eta**2)/(eps1**2 + eta**2)**2/eps2
                 ZDpp(1) = ZDpp(1) - num*(eps2**2 - eta**2)/(eps2**2 + eta**2)**2/eps1
                 
                 eps1 = w + e(i) - e(a) - e(c)
                 eps2 = e(i) + e(j) - e(a) - e(b)
                 num = + ERI(p,i,c,a) * ( ERI(a,b,i,j)*(4d0*ERI(p,b,c,j) - 2d0*ERI(p,b,j,c)) &
                                        + ERI(a,b,j,i)*(    ERI(p,b,j,c) - 2d0*ERI(p,b,c,j)) )
                  
                 ZDpp(2) = ZDpp(2)  - num*(eps1**2 - eta**2)/(eps1**2 + eta**2)**2/eps2
                  
                 num = + ERI(p,i,a,c) * ( ERI(a,b,i,j)*(    ERI(p,b,j,c) - 2d0*ERI(p,b,c,j)) &
                                        + ERI(a,b,j,i)*(    ERI(p,b,c,j) - 2d0*ERI(p,b,j,c)) )

                 ZDpp(2) = ZDpp(2) - num*(eps1**2 - eta**2)/(eps1**2 + eta**2)**2/eps2
                  
              end do
           end do
        end do
     end do
  end do
  ZDpp(3) = ZDpp(2)

!------------------------------!
! Collecting self-energy terms !  
!------------------------------!

  RGF3_dSigC = RGF3_dSigC + ZCpp(1) + ZCpp(2) + ZCpp(3) + ZCpp(4) + ZCpp(5) + ZCpp(6)
  
  RGF3_dSigC = RGF3_dSigC + ZDpp(1) + ZDpp(2) + ZDpp(3) + ZDpp(4) + ZDpp(5) + ZDpp(6)

end function 
