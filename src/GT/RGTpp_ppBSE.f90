subroutine RGTpp_ppBSE(TDA_T,TDA,dBSE,dTDA,singlet,triplet,eta,nBas,nC,nO,nV,nR,nOOs,nVVs,nOOt,nVVt, &
                      ERI,dipole_int,eT,eGT,EcBSE)

! Compute the Bethe-Salpeter excitation energies with the T-matrix kernel

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: TDA_T
  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA
  logical,intent(in)            :: singlet
  logical,intent(in)            :: triplet

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR

  integer,intent(in)            :: nOOs
  integer,intent(in)            :: nOOt
  integer,intent(in)            :: nVVs
  integer,intent(in)            :: nVVt

  double precision,intent(in)   :: eT(nBas)
  double precision,intent(in)   :: eGT(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)

! Local variables

  integer                       :: ispin

  double precision              :: EcRPA(nspin)
  double precision,allocatable  :: Bpp(:,:),Cpp(:,:),Dpp(:,:)
  double precision,allocatable  :: Om1s(:),Om1t(:)
  double precision,allocatable  :: X1s(:,:),X1t(:,:)
  double precision,allocatable  :: Y1s(:,:),Y1t(:,:)
  double precision,allocatable  :: rho1s(:,:,:),rho1t(:,:,:)
  double precision,allocatable  :: Om2s(:),Om2t(:)
  double precision,allocatable  :: X2s(:,:),X2t(:,:)
  double precision,allocatable  :: Y2s(:,:),Y2t(:,:)
  double precision,allocatable  :: rho2s(:,:,:),rho2t(:,:,:)
  double precision,allocatable  :: TBs(:,:),TCs(:,:),TDs(:,:)
  double precision,allocatable  :: TBt(:,:),TCt(:,:),TDt(:,:)

! Output variables

  double precision,intent(out)  :: EcBSE(nspin)

!------------------------------------!
! Compute T-matrix for singlet block !
!------------------------------------!
 
  ispin = 1

  allocate(Om1s(nVVs),X1s(nVVs,nVVs),Y1s(nOOs,nVVs),Om2s(nOOs),X2s(nVVs,nOOs),Y2s(nOOs,nOOs), &
           Bpp(nVVs,nOOs),Cpp(nVVs,nVVs),Dpp(nOOs,nOOs))
 
  if(.not.TDA_T) call ppLR_B(ispin,nBas,nC,nO,nV,nR,nOOs,nVVs,1d0,ERI,Bpp)
                 call ppLR_C(ispin,nBas,nC,nO,nV,nR,nVVs,1d0,eT,ERI,Cpp)
                 call ppLR_D(ispin,nBas,nC,nO,nV,nR,nOOs,1d0,eT,ERI,Dpp)
 
  call ppLR(TDA_T,nOOs,nVVs,Bpp,Cpp,Dpp,Om1s,X1s,Y1s,Om2s,X2s,Y2s,EcRPA(ispin))
 
  deallocate(Bpp,Cpp,Dpp)
  allocate(TBs(nVVs,nOOs),TCs(nVVs,nVVs),TDs(nOOs,nOOs))
 
  if(.not.TDA_T) call RGTpp_ppBSE_static_kernel_B(ispin,eta,nBas,nC,nO,nV,nR,nOOs,nVVs,nOOs,nVVs,1d0, & 
                                                  Om1s,rho1s,Om2s,rho2s,TBs)
                 call RGTpp_ppBSE_static_kernel_C(ispin,eta,nBas,nC,nO,nV,nR,nOOs,nVVs,nOOs,nVVs,1d0, & 
                                                  Om1s,rho1s,Om2s,rho2s,TCs)
                 call RGTpp_ppBSE_static_kernel_D(ispin,eta,nBas,nC,nO,nV,nR,nOOs,nVVs,nOOs,nVVs,1d0, & 
                                                  Om1s,rho1s,Om2s,rho2s,TDs)

  deallocate(Om1s,X1s,Y1s,Om2s,X2s,Y2s)

!------------------------------------!
! Compute T-matrix for triplet block !
!------------------------------------!

  ispin = 2

  allocate(Om1t(nVVt),X1t(nVVt,nVVt),Y1t(nOOt,nVVt),Om2t(nOOt),X2t(nVVt,nOOt),Y2t(nOOt,nOOt), &
           Bpp(nVVt,nOOt),Cpp(nVVt,nVVt),Dpp(nOOt,nOOt))

  if(.not.TDA_T) call ppLR_B(ispin,nBas,nC,nO,nV,nR,nOOt,nVVt,1d0,ERI,Bpp)
                 call ppLR_C(ispin,nBas,nC,nO,nV,nR,nVVt,1d0,eT,ERI,Cpp)
                 call ppLR_D(ispin,nBas,nC,nO,nV,nR,nOOt,1d0,eT,ERI,Dpp)

  call ppLR(TDA_T,nOOt,nVVt,Bpp,Cpp,Dpp,Om1t,X1t,Y1t,Om2t,X2t,Y2t,EcRPA(ispin))

  deallocate(Bpp,Cpp,Dpp)
  allocate(TBt(nVVt,nOOt),TCt(nVVt,nVVt),TDt(nOOt,nOOt))

  if(.not.TDA_T) call RGTpp_ppBSE_static_kernel_B(ispin,eta,nBas,nC,nO,nV,nR,nOOt,nVVt,nOOs,nVVs,1d0, & 
                                                 Om1t,rho1t,Om2t,rho2t,TBt)
                 call RGTpp_ppBSE_static_kernel_C(ispin,eta,nBas,nC,nO,nV,nR,nOOt,nVVt,nOOs,nVVs,1d0, & 
                                                 Om1t,rho1t,Om2t,rho2t,TCt)
                 call RGTpp_ppBSE_static_kernel_D(ispin,eta,nBas,nC,nO,nV,nR,nOOt,nVVt,nOOs,nVVs,1d0, & 
                                                 Om1t,rho1t,Om2t,rho2t,TDt)

  deallocate(Om1t,X1t,Y1t,Om2t,X2t,Y2t)

!------------------!
! Singlet manifold !
!------------------!

  if(singlet) then

    ispin = 1

    allocate(Om1s(nVVs),X1s(nVVs,nVVs),Y1s(nOOs,nVVs),Om2s(nOOs),X2s(nVVs,nOOs),Y2s(nOOs,nOOs), &
             Bpp(nVVs,nOOs),Cpp(nVVs,nVVs),Dpp(nOOs,nOOs))

    if(.not.TDA) call ppLR_B(ispin,nBas,nC,nO,nV,nR,nOOs,nVVs,1d0,ERI,Bpp)
                 call ppLR_C(ispin,nBas,nC,nO,nV,nR,nVVs,1d0,eGT,ERI,Cpp)
                 call ppLR_D(ispin,nBas,nC,nO,nV,nR,nOOs,1d0,eGT,ERI,Dpp)

    Bpp(:,:) = Bpp(:,:) - TBs(:,:) - TBt(:,:) 
    Cpp(:,:) = Cpp(:,:) - TCs(:,:) - TCt(:,:) 
    Dpp(:,:) = Dpp(:,:) - TDs(:,:) - TDt(:,:) 

    call ppLR(TDA,nOOs,nVVs,Bpp,Cpp,Dpp,Om1s,X1s,Y1s,Om2s,X2s,Y2s,EcBSE(ispin))

    call ppLR_transition_vectors(.true.,nBas,nC,nO,nV,nR,nOOs,nVVs,dipole_int,Om1s,X1s,Y1s,Om2s,X2s,Y2s)

    deallocate(Om1s,X1s,Y1s,Om2s,X2s,Y2s,Bpp,Cpp,Dpp) 

  end if

!------------------!
! Triplet manifold !
!------------------!

  if(triplet) then

    ispin  = 2

    EcBSE(ispin) = 0d0

    allocate(Om1t(nVVt),X1t(nVVt,nVVt),Y1t(nOOt,nVVt),Om2t(nOOt),X2t(nVVt,nOOt),Y2t(nOOt,nOOt), &
             Bpp(nVVt,nOOt),Cpp(nVVt,nVVt),Dpp(nOOt,nOOt))
  
    if(.not.TDA) call ppLR_B(ispin,nBas,nC,nO,nV,nR,nOOs,nVVs,1d0,ERI,Bpp)
                 call ppLR_C(ispin,nBas,nC,nO,nV,nR,nVVs,1d0,eGT,ERI,Cpp)
                 call ppLR_D(ispin,nBas,nC,nO,nV,nR,nOOs,1d0,eGT,ERI,Dpp)

    Bpp(:,:) = Bpp(:,:) + TBs(:,:) - TBt(:,:)
    Cpp(:,:) = Cpp(:,:) + TCs(:,:) - TCt(:,:)
    Dpp(:,:) = Dpp(:,:) + TDs(:,:) - TDt(:,:)

    call ppLR(TDA,nOOs,nVVs,Bpp,Cpp,Dpp,Om1t,X1t,Y1t,Om2t,X2t,Y2t,EcBSE(ispin))

    call ppLR_transition_vectors(.false.,nBas,nC,nO,nV,nR,nOOt,nVVt,dipole_int,Om1t,X1t,Y1t,Om2t,X2t,Y2t)

    deallocate(Om1t,X1t,Y1t,Om2t,X2t,Y2t,Bpp,Cpp,Dpp)

  end if

end subroutine 
