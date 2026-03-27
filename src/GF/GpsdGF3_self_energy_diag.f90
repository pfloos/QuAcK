subroutine GpsdGF3_self_energy_diag(eta,nBas,nC,nO,nV,nR,e,ERI,SigC,Z)

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
  integer                       :: klc,kcd,ija,ijb,iab,jab
  integer                       :: p,q
  integer                       :: n2h1p,n2p1h
  double precision              :: eps,eps1,eps2
  double precision              :: num
  double precision,allocatable  :: SigInf(:)

  double precision,allocatable  :: C1_2h1p(:,:)
  double precision,allocatable  :: C1_2p1h(:,:)
  double precision,allocatable  :: K_2h1p(:,:),dK_2h1p(:,:)
  double precision,allocatable  :: K_2p1h(:,:),dK_2p1h(:,:)
  double precision,allocatable  :: U1_2h1p(:,:),U2_2h1p(:,:)
  double precision,allocatable  :: U1_2p1h(:,:),U2_2p1h(:,:)
  
! Output variables

  double precision,intent(out)  :: SigC(nBas)
  double precision,intent(out)  :: Z(nBas)

! Dimension of the 2h1p and 2p1h subspaces

  n2h1p = nO*nO*nV
  n2p1h = nV*nV*nO
  
! Memory allocation

  allocate(C1_2h1p(n2h1p,n2h1p))
  allocate(C1_2p1h(n2p1h,n2p1h))
  
  allocate(K_2h1p(nBas,n2h1p),dK_2h1p(nBas,n2h1p))
  allocate(K_2p1h(nBas,n2p1h),dK_2p1h(nBas,n2p1h))
  
  allocate(U1_2h1p(nBas,n2h1p),U2_2h1p(nBas,n2h1p))
  allocate(U1_2p1h(nBas,n2p1h),U2_2p1h(nBas,n2p1h))
  
  allocate(SigInf(nBas))

! Initialize 

  C1_2h1p(:,:) = 0d0
  C1_2p1h(:,:) = 0d0
  U1_2h1p(:,:) = 0d0
  U1_2p1h(:,:) = 0d0
  U2_2h1p(:,:) = 0d0
  U2_2p1h(:,:) = 0d0
  K_2h1p(:,:)  = 0d0
  K_2p1h(:,:)  = 0d0
  dK_2h1p(:,:) = 0d0
  dK_2p1h(:,:) = 0d0

  SigInf(:) = 0d0

  SigC(:) = 0d0
  Z(:)    = 0d0

!-----------------------------!
!    Diagonal K 2h1p block    !
!-----------------------------!
  do p=nC+1,nBas-nR
     ija = 0
     do i=nC+1,nO
        do j=nC+1,nO
           do a=nO+1,nBas-nR
              ija = ija + 1
           
              !---------------------------!
              ! Zeroth-order contribution !
              !---------------------------!
              K_2h1p(p,ija)  =   1d0 / (e(p) + e(a) - e(i) - e(j))
              dK_2h1p(p,ija) = - 1d0 / (e(p) + e(a) - e(i) - e(j))**2
           
           end do
        end do
     end do
  end do
!-----------------------------!
!    Diagonal C 2h1p block    !
!-----------------------------!
  ija = 0
  do i=nC+1,nO
     do j=nC+1,nO
        do a=nO+1,nBas-nR
           ija = ija + 1

           klc = 0
           do k=nC+1,nO
              do l=nC+1,nO
                 do c=nO+1,nBas-nR
                    klc = klc + 1
           
                    !---------------------------!
                    !  First-order contribution !
                    !---------------------------!

                    !---------------------------!
                    ! Second-order contribution !
                    !---------------------------!

                 end do
              end do
           end do
           
        end do
     end do
  end do
!-----------------------------!
!    Diagonal K 2p1h block    !
!-----------------------------!
  do p=nC+1,nBas-nR
     iab = 0
     do i=nC+1,nO
        do a=nO+1,nBas-nR
           do b=nO+1,nBas-nR
              iab = iab + 1
           
              !---------------------------!
              ! Zeroth-order contribution !
              !---------------------------!
              K_2p1h(p,iab)  =   1d0 / (e(p) + e(i) - e(a) - e(b))
              dK_2p1h(p,iab) = - 1d0 / (e(p) + e(i) - e(a) - e(b))**2
           
           end do
        end do
     end do
  end do
!-----------------------------!
!    Diagonal C 2p1h block    !
!-----------------------------!
  iab = 0
  do i=nC+1,nO
     do a=nO+1,nBas-nR
        do b=nO+1,nBas-nR
           iab = 0
           
           kcd = 0
           do k=nC+1,nO
              do c=nO+1,nBas-nR
                 do d=nO+1,nBas-nR
                    kcd = kcd + 1
                    !---------------------------!
                    !  First-order contribution !
                    !---------------------------!

                    !---------------------------!
                    ! Second-order contribution !
                    !---------------------------!
           
                 end do
              end do
           end do
           
        end do
     end do
  end do
!---------------------------!
!    2h1p coupling block    !
!---------------------------!
  do p=nC+1,nBas-nR
     ija = 0
     do i=nC+1,nO
        do j=nC+1,nO
           do a=nO+1,nBas-nR
              ija = ija + 1
              !---------------------------!
              !  First-order contribution !
              !---------------------------!
              U1_2h1p(p,ija) = sqrt(0.5d0) * (ERI(p,a,i,j) - ERI(p,a,j,i))

        end do
      end do
    end do
  end do
!---------------------------!
!    2p1h coupling block    !
!---------------------------!
  do p=nC+1,nBas-nR
     iab = 0
     do i=nC+1,nO
        do a=nO+1,nBas-nR
           do b=nO+1,nBas-nR
              iab = iab + 1
              !---------------------------!
              !  First-order contribution !
              !---------------------------!
              U1_2p1h(p,iab) = sqrt(0.5d0) * (ERI(p,i,a,b) - ERI(p,i,b,a))

        end do
      end do
    end do
  end do

!----------------------------!
!    Building self-energy    !
!----------------------------!
  do p=nC+1,nBas-nR
     ija = 0
     do i=nC+1,nO
        do j=nC+1,nO
           do a=nO+1,nBas-nR
              ija = ija + 1
              
              SigC(p) = SigC(p) + U1_2h1p(p,ija) *  K_2h1p(p,ija) * U1_2h1p(p,ija)
              Z(p)    = Z(p)    + U1_2h1p(p,ija) * dK_2h1p(p,ija) * U1_2h1p(p,ija)
              
           end do
        end do
     end do
     iab = 0
     do i=nC+1,nO
        do a=nO+1,nBas-nR
           do b=nO+1,nBas-nR
              iab = iab + 1
              
              SigC(p) = SigC(p) + U1_2p1h(p,iab) *  K_2p1h(p,iab) * U1_2p1h(p,iab)
              Z(p)    = Z(p)    + U1_2p1h(p,iab) * dK_2p1h(p,iab) * U1_2p1h(p,iab)
              
           end do
        end do
     end do

  end do
  
  Z(:) = 1d0/(1d0 - Z(:))

end subroutine 
