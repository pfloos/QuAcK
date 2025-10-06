subroutine Generalized_Fock_RHFB(nBas,nBas_twice,nOrb,nOrb_twice,ENuc,sigma,c,Hc,H_HFB_ao,Occ,ERI)

! Building the Genealized Fock operator

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nBas_twice
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nOrb_twice
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: sigma
  double precision,intent(in)   :: Occ(nOrb) 
  double precision,intent(in)   :: c(nBas,nOrb) 
  double precision,intent(in)   :: Hc(nBas,nBas) 
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: H_HFB_ao(nBas_twice,nBas_twice)

! Local variables

  integer                        :: iorb,jorb

  double precision               :: Ehcore,Vee,Eelec_lambda,trace_1rdm
  double precision               :: sqrt_hole1,sqrt_hole2

  double precision,allocatable   :: eigVAL(:)
  double precision,allocatable   :: eigVALmo(:)
  double precision,allocatable   :: eigVEC(:,:)
  double precision,allocatable   :: eigVECmo(:,:)
  double precision,allocatable   :: Lambdas(:,:)
  double precision,allocatable   :: sqrt_occ(:)
  double precision,allocatable   :: DM2_iiii(:)
  double precision,allocatable   :: DM2_J(:,:)
  double precision,allocatable   :: DM2_K(:,:)
  double precision,allocatable   :: DM2_L(:,:)
  double precision,allocatable   :: Ptmp(:,:)
  double precision,allocatable   :: hCORE(:,:)
  double precision,allocatable   :: c_tmp(:,:)
  double precision,allocatable   :: c_ao(:,:)
  double precision,allocatable   :: R_loc(:,:)
  double precision,allocatable   :: ERImol(:,:,:,:)

! Output variables


  write(*,*)
  write(*,*)'**************************************************'
  write(*,*)'* Building the generalized Fock operator for HFB *'
  write(*,*)'**************************************************'
  write(*,*)

  allocate(eigVAL(nOrb_twice))
  allocate(eigVEC(nOrb_twice,nOrb_twice))

  allocate(eigVALmo(nOrb))
  allocate(eigVECmo(nOrb,nOrb))
  allocate(Ptmp(nBas,nBas))
  allocate(sqrt_occ(nOrb))
  allocate(DM2_iiii(nOrb))
  allocate(Lambdas(nOrb,nOrb))
  allocate(DM2_J(nOrb,nOrb))
  allocate(DM2_K(nOrb,nOrb))
  allocate(DM2_L(nOrb,nOrb))
  allocate(hCORE(nOrb,nOrb))
  allocate(c_tmp(nBas,nOrb))
  allocate(ERImol(nOrb,nOrb,nOrb,nOrb))

  allocate(R_loc(nOrb_twice,nOrb_twice))
  allocate(c_ao(nBas_twice,nOrb_twice))

  c_tmp=c
  hCORE=matmul(transpose(c_tmp),matmul(Hc,c_tmp))
  call AOtoMO_ERI_RHF(nBas,nOrb,c,ERI,ERImol)

  DM2_J=0d0; DM2_K=0d0; DM2_L=0d0; DM2_iiii=0d0; Lambdas=0d0;
  sqrt_occ(:) = sqrt(abs(occ(:)))
  
  do iorb=1,nOrb
   sqrt_hole1=sqrt(abs(1d0-Occ(iorb)))
   do jorb=1,nOrb
    sqrt_hole2=sqrt(abs(1d0-Occ(jorb)))
    DM2_J(iorb,jorb) = 2d0*Occ(iorb)*Occ(jorb)
    DM2_K(iorb,jorb) = -Occ(iorb)*Occ(jorb)
    DM2_L(iorb,jorb) = sigma*sqrt_occ(iorb)*sqrt_occ(jorb)*sqrt_hole1*sqrt_hole2
   enddo
  enddo

  do iorb=1,nOrb
   DM2_iiii(iorb)=Occ(iorb)*Occ(iorb)+sigma*(Occ(iorb)-Occ(iorb)*Occ(iorb))
   DM2_J(iorb,iorb)=0d0
   DM2_K(iorb,iorb)=0d0
   DM2_L(iorb,iorb)=0d0
  enddo

  Ehcore=0d0; Vee=0d0; Eelec_lambda=0d0;
  do iorb=1,nOrb
   Eelec_lambda=Eelec_lambda+Occ(iorb)*hCORE(iorb,iorb)
   Ehcore=Ehcore+2d0*Occ(iorb)*hCORE(iorb,iorb)
   Vee=Vee+DM2_iiii(iorb)*ERImol(iorb,iorb,iorb,iorb)
   Lambdas(iorb,:)=Occ(iorb)*hCORE(:,iorb)                                   ! Init: Lambda_pq = n_p hCORE_qp
   Lambdas(iorb,:)=Lambdas(iorb,:)+DM2_iiii(iorb)*ERImol(:,iorb,iorb,iorb)   ! any->iorb,iorb->iorb
   do jorb=1,nOrb
    if(iorb/=jorb) then
     Lambdas(iorb,:)=Lambdas(iorb,:)+DM2_J(iorb,jorb)*ERImol(:,jorb,iorb,jorb) ! any->iorb,jorb->jorb
     Lambdas(iorb,:)=Lambdas(iorb,:)+DM2_K(iorb,jorb)*ERImol(:,jorb,jorb,iorb) ! any->jorb,jorb->iorb
     Lambdas(iorb,:)=Lambdas(iorb,:)+DM2_L(iorb,jorb)*ERImol(:,iorb,jorb,jorb) ! any->jorb,iorb->jorb
    endif
    Vee=Vee+DM2_J(iorb,jorb)*ERImol(iorb,jorb,iorb,jorb)+DM2_K(iorb,jorb)*ERImol(iorb,jorb,jorb,iorb) &
       +DM2_L(iorb,jorb)*ERImol(iorb,iorb,jorb,jorb)
   enddo
   Eelec_lambda=Eelec_lambda+Lambdas(iorb,iorb) 
  enddo

  write(*,'(a)')       '  Energy contributions (a.u.) '
  write(*,'(a)')       '  --------------------------- '
  write(*,'(a,f15.8)') '  Hcore                       ', Ehcore
  write(*,'(a,f15.8)') '  Vee                         ', Vee
  write(*,'(a,f15.8)') '  Eelectronic(np,hpp,lambda)  ', Eelec_lambda
  write(*,'(a,f15.8)') '  Etot                        ', Eelec_lambda+ENuc
  write(*,*)
  write(*,*) '  Gen. Fock matrix'
  do iorb=1,nOrb
   write(*,'(*(f10.5))') Lambdas(iorb,1:nOrb)
  enddo

  eigVECmo=Lambdas
  call diagonalize_matrix(nOrb,eigVECmo,eigVALmo)

  write(*,*) '  Gen. Fock matrix eigenvectors'
  do iorb=1,nOrb
   write(*,'(*(f10.5))') eigVECmo(iorb,1:nOrb)
  enddo
  write(*,*) '  Gen. Fock matrix eigenvalues (a.u.)'
  write(*,'(*(f10.5))') eigVALmo(:)

  c_tmp = matmul(c_tmp,eigVECmo)

  c_ao(:,:) = 0d0
  c_ao(1:nBas      ,1:nOrb      )           = c_tmp(1:nBas,1:nOrb)
  c_ao(nBas+1:nBas_twice,nOrb+1:nOrb_twice) = c_tmp(1:nBas,1:nOrb)
  eigVEC = matmul(transpose(c_ao),matmul(H_HFB_ao,c_ao)) ! This is H_HFB

  call diagonalize_matrix(nOrb_twice,eigVEC,eigVAL)

  ! Build R (as R^no) and save the eigenvectors
    
  trace_1rdm = 0d0 
  R_loc(:,:)     = 0d0
  do iorb=1,nOrb
   R_loc(:,:) = R_loc(:,:) + matmul(eigVEC(:,iorb:iorb),transpose(eigVEC(:,iorb:iorb))) 
  enddo

  write(*,*) '  R in the basis making the Gen. Fock diagonal'
  do iorb=1,nOrb
   write(*,'(*(f10.5))') R_loc(iorb,1:nOrb)
  enddo
  Ptmp=matmul(matmul(c_tmp,R_loc(1:nOrb,1:nOrb)),transpose(c_tmp))
  write(*,*) '  P^ao from the basis making Gen. Fock diagonal'
  do iorb=1,nBas
   write(*,'(*(f10.5))') Ptmp(iorb,1:nOrb)
  enddo
  
  deallocate(DM2_J,DM2_K,DM2_L,sqrt_occ,DM2_iiii,hCORE,ERImol,Lambdas)
  deallocate(Ptmp,c_tmp,c_ao,R_loc)
  deallocate(eigVAL)
  deallocate(eigVEC)
  deallocate(eigVALmo)
  deallocate(eigVECmo)

end subroutine 

