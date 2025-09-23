subroutine linDyson_GTpp_analytic_Sig_RHF(nBas,nOrb,nO,c,eHF,nfreqs,wweight,wcoord,ERI,ERI_MO,&
           Enuc,EcGM,eta,T,V,S,P,Pcorr)


! Use the restricted Sigma_c(E) to compute the linnearized approximation to G

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nfreqs
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nO

  double precision,intent(in)   :: eta
  double precision,intent(in)   :: Enuc
  double precision,intent(inout):: eHF(nOrb)
  double precision,intent(in)   :: wweight(nfreqs)
  double precision,intent(in)   :: wcoord(nfreqs)
  double precision,intent(in)   :: T(nBas,nBas)
  double precision,intent(in)   :: V(nBas,nBas)
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: P(nBas,nBas)
  double precision,intent(in)   :: c(nBas,nOrb)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_MO(nOrb,nOrb,nOrb,nOrb)

! Local variables

  integer                       :: kind_int
  integer                       :: nV
  integer                       :: nOOs,nOOt
  integer                       :: nVVs,nVVt
  integer                       :: ifreq
  integer                       :: iorb,jorb
  integer                       :: ispin
  integer                       :: nOrb2
  integer                       :: nfreqs2

  double precision              :: lim_inf,lim_sup
  double precision              :: alpha,beta
  double precision              :: ET,EV,EJ,EK,ElinG,trace_occ
  double precision,external     :: trace_matrix
  double precision              :: EcRPA
  double precision,allocatable  :: c_inv(:,:)
  double precision,allocatable  :: J(:,:)
  double precision,allocatable  :: K(:,:)
  double precision,allocatable  :: Occ(:)
  double precision,allocatable  :: wweight2(:)
  double precision,allocatable  :: wcoord2(:)
  double precision,allocatable  :: Pcorr_mo(:,:)
  double precision,allocatable  :: Bpp(:,:)
  double precision,allocatable  :: Cpp(:,:)
  double precision,allocatable  :: Dpp(:,:)
  double precision,allocatable  :: Om1s(:),Om1t(:)
  double precision,allocatable  :: X1s(:,:),X1t(:,:)
  double precision,allocatable  :: Y1s(:,:),Y1t(:,:)
  double precision,allocatable  :: rho1s(:,:,:),rho1t(:,:,:)
  double precision,allocatable  :: Om2s(:),Om2t(:)
  double precision,allocatable  :: X2s(:,:),X2t(:,:)
  double precision,allocatable  :: Y2s(:,:),Y2t(:,:)
  double precision,allocatable  :: rho2s(:,:,:),rho2t(:,:,:)

  complex *16,allocatable       :: wcoord2_cpx(:)
  complex *16,allocatable       :: Sigma_c(:,:)
  complex *16,allocatable       :: Tmp_mo(:,:)

! Ouput variables

  double precision,intent(out)  :: EcGM
  double precision,intent(out)  :: Pcorr(nBas,nBas)

! Allocate and initialize arrays and variables

  write(*,*)
  write(*,*) '********************************************'
  write(*,*) '*          G = Go + Go Sigma Go            *'
  write(*,*) '* linearized-Dyson equation approximation  *'
  write(*,*) '*   [ using the analytic Sigma = GT(pp) ]  *'
  write(*,*) '********************************************'
  write(*,*)

  nfreqs2=10*nfreqs
  allocate(Sigma_c(nOrb,nOrb))
  allocate(Tmp_mo(nOrb,nOrb),Pcorr_mo(nOrb,nOrb))

! Prepare second quadrature

  kind_int = 1
  lim_inf = 0d0; lim_sup = 1d0;
  alpha = 0d0;   beta  = 0d0;
  allocate(wweight2(nfreqs2),wcoord2(nfreqs2),wcoord2_cpx(nfreqs2))
  call cgqf(nfreqs2,kind_int,alpha,beta,lim_inf,lim_sup,wcoord2,wweight2)
  wweight2(:)=wweight2(:)/((1d0-wcoord2(:))**2d0)
  wcoord2(:)=wcoord2(:)/(1d0-wcoord2(:))
  wcoord2_cpx(:)=wcoord2(:)*im

! Compute linear response

  nV=nOrb-nO
  nOOs=nO*(nO+1)/2
  nVVs=nV*(nV+1)/2
  nOOt=nO*(nO-1)/2
  nVVt=nV*(nV-1)/2

  allocate(Om1s(nVVs),X1s(nVVs,nVVs),Y1s(nOOs,nVVs),rho1s(nOrb,nOrb,nVVs))
  allocate(Om2s(nOOs),X2s(nVVs,nOOs),Y2s(nOOs,nOOs),rho2s(nOrb,nOrb,nOOs))
  allocate(Om1t(nVVt),X1t(nVVt,nVVt),Y1t(nOOt,nVVt),rho1t(nOrb,nOrb,nVVt))
  allocate(Om2t(nOOt),X2t(nVVt,nOOt),Y2t(nOOt,nOOt),rho2t(nOrb,nOrb,nOOt))

  ispin  = 1
  allocate(Bpp(nVVs,nOOs),Cpp(nVVs,nVVs),Dpp(nOOs,nOOs))
  call ppRLR_B(ispin,nOrb,0,nO,nV,0,nOOs,nVVs,1d0,ERI_MO,Bpp)
  call ppRLR_C(ispin,nOrb,0,nO,nV,0,nVVs,1d0,eHF,ERI_MO,Cpp)
  call ppRLR_D(ispin,nOrb,0,nO,nV,0,nOOs,1d0,eHF,ERI_MO,Dpp)
  call ppRLR(.false.,nOOs,nVVs,Bpp,Cpp,Dpp,Om1s,X1s,Y1s,Om2s,X2s,Y2s,EcRPA)
  deallocate(Bpp,Cpp,Dpp)

  ispin  = 2
  allocate(Bpp(nVVt,nOOt),Cpp(nVVt,nVVt),Dpp(nOOt,nOOt))
  call ppRLR_B(ispin,nOrb,0,nO,nV,0,nOOt,nVVt,1d0,ERI_MO,Bpp)
  call ppRLR_C(ispin,nOrb,0,nO,nV,0,nVVt,1d0,eHF,ERI_MO,Cpp)
  call ppRLR_D(ispin,nOrb,0,nO,nV,0,nOOt,1d0,eHF,ERI_MO,Dpp)
  call ppRLR(.false.,nOOt,nVVt,Bpp,Cpp,Dpp,Om1t,X1t,Y1t,Om2t,X2t,Y2t,EcRPA)
  deallocate(Bpp,Cpp,Dpp)

  ispin = 1
  call RGTpp_excitation_density(ispin,nOrb,0,nO,nV,0,nOOs,nVVs,ERI_MO,X1s,Y1s,rho1s,X2s,Y2s,rho2s)

  ispin = 2
  call RGTpp_excitation_density(ispin,nOrb,0,nO,nV,0,nOOt,nVVt,ERI_MO,X1t,Y1t,rho1t,X2t,Y2t,rho2t)


! Integration along imag. freqs contributions

  Pcorr_mo(:,:)=0d0
  do ifreq=1,nfreqs2
   
   call G_MO_RHF(nOrb,nO,0d0,eHF,wcoord2_cpx(ifreq),Tmp_mo)                                                              ! This is G(iw2)
   call RGTpp_self_energy_iomega(eta,wcoord2_cpx(ifreq),nOrb,0,nO,nV,0,nOOs,nVVs,nOOt,nVVt,eHF,Om1s,rho1s,Om2s,rho2s, &  ! This is Sigma_c(iw2)
   Om1t,rho1t,Om2t,rho2t,EcGM,Sigma_c)
   Tmp_mo(:,:)=matmul(Tmp_mo(:,:),matmul(Sigma_c(:,:),Tmp_mo(:,:)))                                                      ! This is G(iw2) Sigma_c(iw2) G(iw2)
 
   Pcorr_mo(:,:) = Pcorr_mo(:,:) + wweight2(ifreq)*real(Tmp_mo(:,:)+conjg(Tmp_mo(:,:))) ! Integrate along iw2
  
  enddo
  Pcorr_mo(:,:) = Pcorr_mo(:,:)/pi
  
  Pcorr(:,:) = P(:,:) + matmul(c,matmul(Pcorr_mo(:,:),transpose(c)))

! Compute new total energy and Occ numbers

  allocate(J(nBas,nBas))
  allocate(K(nBas,nBas))
  allocate(Occ(nOrb))
  allocate(c_inv(nOrb,nBas))

  c_inv(:,:) = matmul(transpose(c),S)
  Pcorr_mo(:,:) = Pcorr_mo(:,:) + matmul(matmul(c_inv,P),transpose(c_inv)) 
  call diagonalize_matrix(nOrb,Pcorr_mo,Occ)
  call Hartree_matrix_AO_basis(nBas,Pcorr,ERI,J)
  call exchange_matrix_AO_basis(nBas,Pcorr,ERI,K)

  ! Kinetic energy

  ET = trace_matrix(nBas,matmul(Pcorr,T))

  ! Potential energy

  EV = trace_matrix(nBas,matmul(Pcorr,V))

  ! Hartree energy

  EJ = 0.5d0*trace_matrix(nBas,matmul(Pcorr,J))

  ! Exchange energy

  EK = 0.25d0*trace_matrix(nBas,matmul(Pcorr,K))

  ! Total energy (incl. the Galitkii-Migdal contribution)

  ElinG = ET + EV + EJ + EK + EcGM

! Print results

  write(*,*)
  write(*,'(A50)')           '---------------------------------------'
  write(*,'(A33)')           ' Summary               '
  write(*,'(A50)')           '---------------------------------------'
  write(*,'(A33,1X,F16.10,A3)') ' One-electron energy = ',ET + EV,' au'
  write(*,'(A33,1X,F16.10,A3)') ' Kinetic      energy = ',ET,' au'
  write(*,'(A33,1X,F16.10,A3)') ' Potential    energy = ',EV,' au'
  write(*,'(A50)')           '---------------------------------------'
  write(*,'(A33,1X,F16.10,A3)') ' Two-electron energy = ',EJ + EK + EcGM,' au'
  write(*,'(A33,1X,F16.10,A3)') ' Hartree      energy = ',EJ,' au'
  write(*,'(A33,1X,F16.10,A3)') ' Exchange     energy = ',EK,' au'
  write(*,'(A33,1X,F16.10,A3)') ' GM           energy = ',EcGM,' au'
  write(*,'(A50)')           '---------------------------------------'
  write(*,'(A33,1X,F16.10,A3)') ' Electronic   energy = ',ElinG,' au'
  write(*,'(A33,1X,F16.10,A3)') ' Nuclear   repulsion = ',ENuc,' au'
  write(*,'(A33,1X,F16.10,A3)') ' lin-Dys GTpp energy = ',ElinG + ENuc,' au'
  write(*,'(A50)')           '---------------------------------------'
  write(*,'(A50)') '---------------------------------------'
  write(*,'(A50)') ' lin-Dyson GTpp occupation numbers '
  write(*,'(A50)') '---------------------------------------'
  trace_occ=0d0
  do iorb=1,nOrb
   if(abs(Occ(iorb))>1d-8) then
    write(*,'(I7,10F15.8)') iorb,Occ(nOrb-(iorb-1))
   endif
   trace_occ=trace_occ+Occ(iorb)
  enddo
  write(*,*)
  write(*,'(A33,1X,F16.10,A3)') ' Trace [ 1D ]        = ',trace_occ,'   '
  write(*,*)

  ! Deallocate arrays

  deallocate(J,K,Occ,c_inv)
  deallocate(Om1s,X1s,Y1s,rho1s)
  deallocate(Om2s,X2s,Y2s,rho2s)
  deallocate(Om1t,X1t,Y1t,rho1t)
  deallocate(Om2t,X2t,Y2t,rho2t)
  deallocate(Sigma_c,Tmp_mo,Pcorr_mo,wweight2,wcoord2,wcoord2_cpx)

end subroutine

