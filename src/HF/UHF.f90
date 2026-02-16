subroutine UHF(dotest,maxSCF,thresh,max_diis,guess_type,mix,level_shift,writeMOs,nNuc,ZNuc,rNuc,ENuc, & 
               nBas,nO,S,T,V,Hc,ERI,dipole_int,X,EUHF,eHF,c,P,F)

! Perform unrestricted Hartree-Fock calculation

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dotest,writeMOs

  integer,intent(in)            :: maxSCF
  integer,intent(in)            :: max_diis
  integer,intent(in)            :: guess_type
  double precision,intent(in)   :: mix 
  double precision,intent(in)   :: level_shift
  double precision,intent(in)   :: thresh
  integer,intent(in)            :: nBas

  integer,intent(in)            :: nNuc
  double precision,intent(in)   :: ZNuc(nNuc)
  double precision,intent(in)   :: rNuc(nNuc,ncart)
  double precision,intent(in)   :: ENuc

  integer,intent(in)            :: nO(nspin)
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: T(nBas,nBas)
  double precision,intent(in)   :: V(nBas,nBas)
  double precision,intent(in)   :: Hc(nBas,nBas) 
  double precision,intent(in)   :: X(nBas,nBas) 
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)

! Local variables

  logical                       :: file_exists

  integer                       :: iorb
  integer                       :: nSCF
  integer                       :: nBasSq
  integer                       :: n_diis
  double precision              :: Conv
  double precision              :: rcond(nspin)
  double precision              :: ET(nspin)
  double precision              :: EV(nspin)
  double precision              :: EJ(nsp)
  double precision              :: EK(nspin)
  double precision              :: dipole(ncart)

  double precision,allocatable  :: cp(:,:,:)
  double precision,allocatable  :: J(:,:,:)
  double precision,allocatable  :: Fp(:,:,:)
  double precision,allocatable  :: K(:,:,:)
  double precision,allocatable  :: err(:,:,:)
  double precision,allocatable  :: err_diis(:,:,:)
  double precision,allocatable  :: F_diis(:,:,:)
  double precision,allocatable  :: tmp(:,:)
  double precision,external     :: trace_matrix

  integer                       :: ispin

! Output variables

  double precision,intent(out)  :: EUHF
  double precision,intent(out)  :: eHF(nBas,nspin)
  double precision,intent(inout):: c(nBas,nBas,nspin)
  double precision,intent(out)  :: P(nBas,nBas,nspin)
  double precision,intent(out)  :: F(nBas,nBas,nspin)

! Hello world

  write(*,*)
  write(*,*)'*******************************'
  write(*,*)'* Unrestricted HF Calculation *'
  write(*,*)'*******************************'
  write(*,*)

! Useful stuff

  nBasSq = nBas*nBas

! Memory allocation

  allocate(J(nBas,nBas,nspin),K(nBas,nBas,nspin),Fp(nBas,nBas,nspin), &
           err(nBas,nBas,nspin),cp(nBas,nBas,nspin),                  &
           err_diis(nBasSq,max_diis,nspin),F_diis(nBasSq,max_diis,nspin))

! Guess coefficients and demsity matrices

  do ispin=1,nspin
    call mo_guess(nBas,nBas,guess_type,S,Hc,X,c(:,:,ispin))
  end do
  if(guess_type ==6) then
    allocate(tmp(nBas,nBas))
    call read_matin(nBas,nBas,tmp,"real_MOs_alpha.dat")
    c(:,:,1) = tmp
    call read_matin(nBas,nBas,tmp,"real_MOs_beta.dat")
    c(:,:,2) = tmp
    deallocate(tmp)
  end if
  do ispin=1,nspin
    P(:,:,ispin) = matmul(c(:,1:nO(ispin),ispin),transpose(c(:,1:nO(ispin),ispin)))
  end do

! Initialization

  n_diis          = 0
  F_diis(:,:,:)   = 0d0
  err_diis(:,:,:) = 0d0
  rcond(:)        = 0d0

  nSCF = 0
  Conv = 1d0

!------------------------------------------------------------------------
! Main SCF loop
!------------------------------------------------------------------------

  write(*,*)
  write(*,*)'-----------------------------------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A16,1X,A1,1X,A16,1X,A1,1X,A16,1X,A1,1X,A10,1X,A1,1X)') &
            '|','#','|','E(UHF)','|','EJ(UHF)','|','EK(UHF)','|','Conv','|'
  write(*,*)'-----------------------------------------------------------------------------'
  
  do while(Conv > thresh .and. nSCF < maxSCF)

!   Increment 

    nSCF = nSCF + 1

!   Build Hartree repulsion

    do ispin=1,nspin
      call Hartree_matrix_AO_basis(nBas,P(:,:,ispin),ERI(:,:,:,:),J(:,:,ispin))
    end do

!   Compute exchange potential

    do ispin=1,nspin
      call exchange_matrix_AO_basis(nBas,P(:,:,ispin),ERI(:,:,:,:),K(:,:,ispin))
    end do
 
!   Build Fock operator

    do ispin=1,nspin
      F(:,:,ispin) = Hc(:,:) + J(:,:,ispin) + J(:,:,mod(ispin,2)+1) + K(:,:,ispin)
    end do

!   Check convergence 

    do ispin=1,nspin
      err(:,:,ispin) = matmul(F(:,:,ispin),matmul(P(:,:,ispin),S(:,:))) - matmul(matmul(S(:,:),P(:,:,ispin)),F(:,:,ispin))
    end do

    if(nSCF > 2) Conv = maxval(abs(err))
 
!   Kinetic energy

    do ispin=1,nspin
      ET(ispin) = trace_matrix(nBas,matmul(P(:,:,ispin),T))
    end do

!   Potential energy

    do ispin=1,nspin
      EV(ispin) = trace_matrix(nBas,matmul(P(:,:,ispin),V))
    end do

!   Hartree energy

    EJ(1) = 0.5d0*trace_matrix(nBas,matmul(P(:,:,1),J(:,:,1)))
    EJ(2) = trace_matrix(nBas,matmul(P(:,:,1),J(:,:,2)))
    EJ(3) = 0.5d0*trace_matrix(nBas,matmul(P(:,:,2),J(:,:,2)))

!   Exchange energy

    do ispin=1,nspin
      EK(ispin) = 0.5d0*trace_matrix(nBas,matmul(P(:,:,ispin),K(:,:,ispin)))
    end do

!   Total energy

    EUHF = sum(ET) + sum(EV) + sum(EJ) + sum(EK)

!   DIIS extrapolation

    if(max_diis > 1) then

      n_diis = min(n_diis+1,max_diis)
      do ispin=1,nspin
        if(nO(ispin) > 1) call DIIS_extrapolation(rcond(ispin),nBasSq,nBasSq,n_diis,err_diis(:,1:n_diis,ispin), &
                                                F_diis(:,1:n_diis,ispin),err(:,:,ispin),F(:,:,ispin))
      end do

    end if

!   Level-shifting

    if(level_shift > 0d0 .and. Conv > thresh) then

      do ispin=1,nspin
        call level_shifting(level_shift,nBas,nBas,nO(ispin),S,c(:,:,ispin),F(:,:,ispin))
      end do

    end if

!  Transform Fock matrix in orthogonal basis

    do ispin=1,nspin
      Fp(:,:,ispin) = matmul(transpose(X),matmul(F(:,:,ispin),X))
    end do

!  Diagonalize Fock matrix to get eigenvectors and eigenvalues

    cp(:,:,:) = Fp(:,:,:)
    do ispin=1,nspin
      call diagonalize_matrix(nBas,cp(:,:,ispin),eHF(:,ispin))
    end do

!   Back-transform eigenvectors in non-orthogonal basis

    do ispin=1,nspin
      c(:,:,ispin) = matmul(X,cp(:,:,ispin))
    end do

!   Mix guess for UHF solution in singlet states

    if(nSCF == 1 .and. mix > 0d0) call mix_guess(nBas,nO,mix,c)

!   Compute density matrix 

    do ispin=1,nspin
      P(:,:,ispin) = matmul(c(:,1:nO(ispin),ispin),transpose(c(:,1:nO(ispin),ispin)))
    end do
 
!   Dump results

    write(*,'(1X,A1,1X,I3,1X,A1,1X,F16.10,1X,A1,1X,F16.10,1X,A1,1X,F16.10,1X,A1,1X,E10.2,1X,A1,1X)') &
      '|',nSCF,'|',EUHF + ENuc,'|',sum(EJ),'|',sum(EK),'|',Conv,'|'
 
  end do
  write(*,*)'-----------------------------------------------------------------------------'
!------------------------------------------------------------------------
! End of SCF loop
!------------------------------------------------------------------------

! Did it actually converge?

  if(nSCF == maxSCF) then

    write(*,*)
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)'                 Convergence failed                 '
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)

    stop

  end if

! Compute final UHF energy

  call dipole_moment(nBas,P(:,:,1)+P(:,:,2),nNuc,ZNuc,rNuc,dipole_int,dipole)
  call print_UHF(nBas,nO,S,eHF,c,P,ENuc,ET,EV,EJ,EK,EUHF,dipole)

! Write MOs

  if(writeMOs) then
    call write_matout(nBas,nBas,c(:,:,1),'real_MOs_alpha.dat')
    call write_matout(nBas,nBas,c(:,:,2),'real_MOs_beta.dat')
    call write_matout(nBas,nBas,0*c(:,:,1),'imag_MOs_alpha.dat')
    call write_matout(nBas,nBas,0*c(:,:,2),'imag_MOs_beta.dat')
  endif

! Print test values

  if(dotest) then

    call dump_test_value('U','UHF energy',EUHF)
    call dump_test_value('U','UHF HOMOa energy',eHF(nO(1),1))
    call dump_test_value('U','UHF HOMOb energy',eHF(nO(2),2))
    call dump_test_value('U','UHF LUMOa energy',eHF(nO(1)+1,1))
    call dump_test_value('U','UHF LUMOb energy',eHF(nO(2)+1,2))
    call dump_test_value('U','UHF dipole moment',norm2(dipole))

  end if

end subroutine 
