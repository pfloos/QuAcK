subroutine G0F2(BSE,TDA,singlet_manifold,triplet_manifold,linearize,eta,nBas,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,eHF)

! Perform a one-shot second-order Green function calculation

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: BSE
  logical,intent(in)            :: TDA
  logical,intent(in)            :: singlet_manifold
  logical,intent(in)            :: triplet_manifold
  logical,intent(in)            :: linearize
  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nO
  integer,intent(in)            :: nC
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: ERHF
  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  double precision              :: eps
  double precision              :: V
  double precision              :: EcBSE(nspin)
  double precision,allocatable  :: eGF2(:)
  double precision,allocatable  :: Sig(:)
  double precision,allocatable  :: Z(:)

  integer                       :: i,j,a,b,p

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'|     One-shot second-order Green function     |'
  write(*,*)'************************************************'
  write(*,*)

! Memory allocation

  allocate(Sig(nBas),Z(nBas),eGF2(nBas))

  if(linearize) then 
  
     write(*,*) '*** Quasiparticle equation will be linearized ***'
     write(*,*)

  end  if

! Frequency-dependent second-order contribution

  Sig(:) = 0d0
  Z(:)   = 0d0

  do p=nC+1,nBas-nR
    do i=nC+1,nO
      do j=nC+1,nO
        do a=nO+1,nBas-nR

          eps = eHF(p) + eHF(a) - eHF(i) - eHF(j)
          V  = (2d0*ERI(p,a,i,j) - ERI(p,a,j,i))*ERI(p,a,i,j)
          Sig(p) = Sig(p) + V*eps/(eps**2 + eta**2)
          Z(p)   = Z(p)   - V*(eps**2 - eta**2)/(eps**2 + eta**2)**2

        end do
      end do
    end do
  end do

  do p=nC+1,nBas-nR
    do i=nC+1,nO
      do a=nO+1,nBas-nR
        do b=nO+1,nBas-nR

          eps = eHF(p) + eHF(i) - eHF(a) - eHF(b)
          V  = (2d0*ERI(p,i,a,b) - ERI(p,i,b,a))*ERI(p,i,a,b)
          Sig(p) = Sig(p) + V*eps/(eps**2 + eta**2)
          Z(p)   = Z(p)   - V*(eps**2 - eta**2)/(eps**2 + eta**2)**2

        end do
      end do
    end do
  end do

  Z(:) = 1d0/(1d0 - Z(:))

  if(linearize) then

    eGF2(:) = eHF(:) + Z(:)*Sig(:)

  else

    eGF2(:) = eHF(:) + Sig(:)

  end if

  ! Print results

  call print_G0F2(nBas,nO,eHF,Sig,eGF2,Z)

! Perform BSE2 calculation

  if(BSE) then

    call BSE2(TDA,singlet_manifold,triplet_manifold,eta,nBas,nC,nO,nV,nR,nS,ERI,eHF,eGF2,EcBSE)

    write(*,*)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(2X,A50,F20.10)') 'Tr@BSE2@G0F correlation energy (singlet) =',EcBSE(1)
    write(*,'(2X,A50,F20.10)') 'Tr@BSE2@G0F correlation energy (triplet) =',EcBSE(2)
    write(*,'(2X,A50,F20.10)') 'Tr@BSE2@G0F correlation energy           =',EcBSE(1) + EcBSE(2)
    write(*,'(2X,A50,F20.10)') 'Tr@BSE2@G0F total energy                 =',ENuc + ERHF + EcBSE(1) + EcBSE(2)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)

!   Compute the BSE correlation energy via the adiabatic connection 

!   if(doACFDT) then

!     write(*,*) '------------------------------------------------------'
!     write(*,*) 'Adiabatic connection version of BSE correlation energy'
!     write(*,*) '------------------------------------------------------'
!     write(*,*) 

!     if(doXBS) then 

!       write(*,*) '*** scaled screening version (XBS) ***'
!       write(*,*)

!     end if

!     call ACFDT(exchange_kernel,doXBS,.true.,TDA,BSE,singlet_manifold,triplet_manifold,eta, & 
!                nBas,nC,nO,nV,nR,nS,ERI,eHF,eGW,Omega,XpY,XmY,rho,EcAC)

!     if(exchange_kernel) then
! 
!       EcAC(1) = 0.5d0*EcAC(1)
!       EcAC(2) = 1.5d0*EcAC(1)
! 
!     end if

!     write(*,*)
!     write(*,*)'-------------------------------------------------------------------------------'
!     write(*,'(2X,A50,F20.10)') 'AC@BSE@G0W0 correlation energy (singlet) =',EcAC(1)
!     write(*,'(2X,A50,F20.10)') 'AC@BSE@G0W0 correlation energy (triplet) =',EcAC(2)
!     write(*,'(2X,A50,F20.10)') 'AC@BSE@G0W0 correlation energy           =',EcAC(1) + EcAC(2)
!     write(*,'(2X,A50,F20.10)') 'AC@BSE@G0W0 total energy                 =',ENuc + ERHF + EcAC(1) + EcAC(2)
!     write(*,*)'-------------------------------------------------------------------------------'
!     write(*,*)

!   end if

  end if
end subroutine G0F2
