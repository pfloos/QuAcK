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

  double precision,external     :: Kronecker_delta
  
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
                    C1_2h1p(ija,klc) = - Kronecker_Delta(a,c) * (ERI(i,j,k,l) - ERI(i,j,l,k)) / 2d0  &
                                       + Kronecker_Delta(i,k) * (ERI(c,j,a,l) - ERI(c,j,l,a)) / 2d0 + Kronecker_Delta(j,l) * (ERI(c,i,a,k) - ERI(c,i,k,a)) / 2d0 &
                                       - Kronecker_Delta(i,l) * (ERI(c,j,a,k) - ERI(c,j,k,a)) / 2d0 - Kronecker_Delta(j,k) * (ERI(c,i,a,l) - ERI(c,i,l,a)) / 2d0
                    
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
           iab = iab + 1
           
           kcd = 0
           do k=nC+1,nO
              do c=nO+1,nBas-nR
                 do d=nO+1,nBas-nR
                    kcd = kcd + 1
                    !---------------------------!
                    !  First-order contribution !
                    !---------------------------!
                    C1_2p1h(iab,kcd) = + Kronecker_Delta(i,k) * (ERI(a,b,c,d) - ERI(a,b,d,c)) / 2d0 &
                                       - Kronecker_Delta(a,c) * (ERI(k,b,i,d) - ERI(k,b,d,i)) / 2d0 - Kronecker_Delta(b,d) * (ERI(k,a,i,c) - ERI(k,a,c,i)) / 2d0 &
                                       + Kronecker_Delta(a,d) * (ERI(k,b,i,c) - ERI(k,b,c,i)) / 2d0 + Kronecker_Delta(b,c) * (ERI(k,a,i,d) - ERI(k,a,d,i)) / 2d0
           
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
              !---------------------------!
              ! Second-order contribution !
              !---------------------------!
              do c=nO+1,nBas-nR
                 do d=nO+1,nBas-nR
                    U2_2h1p(p,ija) = U2_2h1p(p,ija) - 0.25d0 * (ERI(i,j,c,d) - ERI(i,j,d,c)) * (ERI(c,d,p,a) - ERI(c,d,a,p)) / (e(c) + e(d) - e(i) - e(j)) / sqrt(0.5d0)
                 end do
              end do
              do k=nC+1,nO
                 do c=nO+1,nBas-nR   
                    U2_2h1p(p,ija) = U2_2h1p(p,ija) - (ERI(i,k,c,a) - ERI(i,k,a,c)) * (ERI(c,j,p,k) - ERI(c,j,k,p)) / (e(i) + e(k) - e(a) - e(c)) / sqrt(0.5d0)
                 end do
              end do
              
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
              !---------------------------!
              ! Second-order contribution !
              !---------------------------!
              do k=nC+1,nO
                 do l=nC+1,nO
                    U2_2p1h(p,iab) = U2_2p1h(p,iab) + 0.25d0 * (ERI(a,b,k,l) - ERI(a,b,l,k)) * (ERI(k,l,p,i) - ERI(k,l,i,p)) / (e(k) + e(l) - e(a) - e(b)) / sqrt(0.5d0)
                 end do
              end do
              do k=nC+1,nO
                 do c=nO+1,nBas-nR   
                    U2_2p1h(p,iab) = U2_2p1h(p,iab) + (ERI(a,c,k,i) - ERI(a,c,i,k)) * (ERI(k,b,p,c) - ERI(k,b,c,p)) / (e(a) + e(c) - e(i) - e(k)) / sqrt(0.5d0)
                 end do
              end do
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

              SigC(p) = SigC(p) + U2_2h1p(p,ija) *  K_2h1p(p,ija) * U1_2h1p(p,ija)
              Z(p)    = Z(p)    + U2_2h1p(p,ija) * dK_2h1p(p,ija) * U1_2h1p(p,ija)

              SigC(p) = SigC(p) + U1_2h1p(p,ija) *  K_2h1p(p,ija) * U2_2h1p(p,ija)
              Z(p)    = Z(p)    + U1_2h1p(p,ija) * dK_2h1p(p,ija) * U2_2h1p(p,ija)

              klc = 0
              do k=nC+1,nO
                 do l=nC+1,nO
                    do c=nO+1,nBas-nR
                       klc = klc + 1
              
                       SigC(p) = SigC(p) + U1_2h1p(p,ija) *  K_2h1p(p,ija) * C1_2h1p(ija,klc) *  K_2h1p(p,klc) * U1_2h1p(p,klc)
                       Z(p)    = Z(p)    + U1_2h1p(p,ija) * dK_2h1p(p,ija) * C1_2h1p(ija,klc) *  K_2h1p(p,klc) * U1_2h1p(p,klc)
                       Z(p)    = Z(p)    + U1_2h1p(p,ija) *  K_2h1p(p,ija) * C1_2h1p(ija,klc) * dK_2h1p(p,klc) * U1_2h1p(p,klc)

                    end do
                 end do
              end do
              
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

              SigC(p) = SigC(p) + U2_2p1h(p,iab) *  K_2p1h(p,iab) * U1_2p1h(p,iab)
              Z(p)    = Z(p)    + U2_2p1h(p,iab) * dK_2p1h(p,iab) * U1_2p1h(p,iab)

              SigC(p) = SigC(p) + U1_2p1h(p,iab) *  K_2p1h(p,iab) * U2_2p1h(p,iab)
              Z(p)    = Z(p)    + U1_2p1h(p,iab) * dK_2p1h(p,iab) * U2_2p1h(p,iab)

              kcd = 0
              do k=nC+1,nO
                 do c=nO+1,nBas-nR
                    do d=nO+1,nBas-nR
                       kcd = kcd + 1

                       SigC(p) = SigC(p) + U1_2p1h(p,iab) *  K_2p1h(p,iab) * C1_2p1h(iab,kcd) *  K_2p1h(p,kcd) * U1_2p1h(p,kcd)
                       Z(p)    = Z(p)    + U1_2p1h(p,iab) * dK_2p1h(p,iab) * C1_2p1h(iab,kcd) *  K_2p1h(p,kcd) * U1_2p1h(p,kcd)
                       Z(p)    = Z(p)    + U1_2p1h(p,iab) *  K_2p1h(p,iab) * C1_2p1h(iab,kcd) * dK_2p1h(p,kcd) * U1_2p1h(p,kcd)
                       
                    end do
                 end do
              end do
                       
           end do
        end do
     end do

  end do
  
  Z(:) = 1d0/(1d0 - Z(:))

end subroutine 
