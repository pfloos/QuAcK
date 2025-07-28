subroutine linDyson_GW_RHF(nBas,nOrb,nO,c,eHF,nfreqs,wweight,wcoord,ERI,vMAT,&
                         Enuc,EcGM,T,V,P,Pcorr)


! Use the restricted Sigma_c(E) to compute the linnearized approximation to G

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nfreqs
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nO

  double precision,intent(in)   :: Enuc
  double precision,intent(in)   :: EcGM
  double precision,intent(inout):: eHF(nOrb)
  double precision,intent(in)   :: wweight(nfreqs)
  double precision,intent(in)   :: wcoord(nfreqs)
  double precision,intent(in)   :: T(nBas,nBas)
  double precision,intent(in)   :: V(nBas,nBas)
  double precision,intent(in)   :: P(nBas,nBas)
  double precision,intent(in)   :: c(nBas,nOrb)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: vMAT(nOrb*nOrb,nOrb*nOrb)

! Local variables

  integer                       :: kind_int
  integer                       :: ifreq
  integer                       :: iorb,jorb
  integer                       :: nOrb2
  integer                       :: nfreqs2

  double precision              :: eta
  double precision              :: lim_inf,lim_sup
  double precision              :: alpha,beta
  double precision              :: ET,EV,EJ,EK,ElinG,trace_occ
  double precision,external     :: trace_matrix
  double precision,allocatable  :: J(:,:)
  double precision,allocatable  :: K(:,:)
  double precision,allocatable  :: Occ(:)
  double precision,allocatable  :: wweight2(:)
  double precision,allocatable  :: wcoord2(:)
  double precision,allocatable  :: Pcorr_mo(:,:)

  complex *16,allocatable       :: wcoord2_cpx(:)
  complex *16,allocatable       :: Sigma_c(:,:,:)
  complex *16,allocatable       :: Tmp_mo(:,:)

! Ouput variables

  double precision,intent(out)  :: Pcorr(nBas,nBas)

! Allocate and initialize arrays and variables

  write(*,*)
  write(*,*) '********************************************'
  write(*,*) '*          G = Go + Go Sigma Go            *'
  write(*,*) '* linearized-Dyson equation approximation  *'
  write(*,*) '********************************************'
  write(*,*)

  eta=0d0
  nfreqs2=10*nfreqs
  allocate(Sigma_c(nfreqs2,nOrb,nOrb))
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

! Build Sigma_c(iw)

  call build_Sigmac_w_RHF(nOrb,nO,nfreqs2,eta,0,wcoord2_cpx,eHF,nfreqs,0,&
                          wweight,wcoord,vMAT,Sigma_c)

! Integration along imag. freqs contributions

  Pcorr_mo(:,:)=0d0
  do ifreq=1,nfreqs2
   
   call G_MO_RHF(nOrb,nO,0d0,eHF,wcoord2_cpx(ifreq),Tmp_mo)                ! This is G(iw2)
   Tmp_mo(:,:)=matmul(Tmp_mo(:,:),matmul(Sigma_c(ifreq,:,:),Tmp_mo(:,:)))  ! This is G(iw2) Sigma_c(iw2) G(iw2)
 
   Pcorr_mo(:,:) = Pcorr_mo(:,:) + wweight2(ifreq)*real(Tmp_mo(:,:)+conjg(Tmp_mo(:,:))) ! Integrate along iw2
  
  enddo
  Pcorr_mo(:,:) = Pcorr_mo(:,:)/pi
  
  Pcorr(:,:) = P(:,:) + matmul(c,matmul(Pcorr_mo(:,:),transpose(c)))

! Compute new total energy and Occ numbers

  allocate(J(nBas,nBas))
  allocate(K(nBas,nBas))
  allocate(Occ(nOrb))

  Pcorr_mo(:,:) = Pcorr_mo(:,:) + matmul(matmul(transpose(c),P),c) 
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
  write(*,'(A33,1X,F16.10,A3)') ' lin-Dyson    energy = ',ElinG + ENuc,' au'
  write(*,'(A50)')           '---------------------------------------'
  write(*,'(A50)') '---------------------------------------'
  write(*,'(A50)') ' lin-Dyson occupation numbers '
  write(*,'(A50)') '---------------------------------------'
  trace_occ=0d0
  do iorb=1,nOrb
   if(abs(Occ(iorb))>1d-8) then
    write(*,'(I7,10F15.8)') iorb,Occ(nOrb-(iorb-1))
   endif
   trace_occ=trace_occ+2d0*Occ(iorb)
  enddo
  write(*,*)
  write(*,'(A33,1X,F16.10,A3)') ' Trace [ 1D ]        = ',trace_occ,'   '
  write(*,*)

  ! Deallocate arrays

  deallocate(J,K,Occ)
  deallocate(Sigma_c,Tmp_mo,Pcorr_mo,wweight2,wcoord2,wcoord2_cpx)

end subroutine

