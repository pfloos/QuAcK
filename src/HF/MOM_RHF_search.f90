subroutine MOM_RHF_search(maxSCF,doaordm,thresh,max_diis,guess_type,level_shift,writeMOs,nNuc,ZNuc,rNuc,ENuc,        &
                      nBas,nOrb,nC,nO,nV,nR,nCVS,FC,S,T,V,Hc,ERI_AO,ERI_MO,dipole_int_AO,dipole_int_MO, & 
                      X,ERHF,e,c,P,F,occupations)

! Search for MOM-RHF solutions

  implicit none
  include 'parameters.h'

  integer,intent(in)            :: maxSCF
  integer,intent(in)            :: max_diis
  integer,intent(in)            :: guess_type
  double precision,intent(in)   :: thresh
  double precision,intent(in)   :: level_shift
  logical,intent(in)            :: doaordm
  logical,intent(in)            :: writeMOs

  integer,intent(in)            :: nBas, nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nCVS
  integer,intent(in)            :: FC
  integer,intent(in)            :: nNuc
  double precision,intent(in)   :: ZNuc(nNuc)
  double precision,intent(in)   :: rNuc(nNuc,ncart)
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: T(nBas,nBas)
  double precision,intent(in)   :: V(nBas,nBas)
  double precision,intent(in)   :: Hc(nBas,nBas)
  double precision,intent(in)   :: X(nBas,nOrb)
  double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas)
  double precision,intent(inout):: ERI_MO(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: dipole_int_AO(nBas,nBas,ncart)
  double precision,intent(inout):: dipole_int_MO(nOrb,nOrb,ncart)
  integer,intent(in)            :: occupations(nO)

! Local variables

  double precision              :: start_HF     ,end_HF       ,t_HF
  double precision              :: start_stab   ,end_stab     ,t_stab
  double precision              :: start_AOtoMO ,end_AOtoMO   ,t_AOtoMO

  logical                       :: unstab,found
  integer                       :: guess
  integer                       :: nS,nSCVS,nFC

  integer,parameter             :: maxS = 20
  integer                       :: ia,i,a,mu
  integer                       :: ispin

  double precision,allocatable  :: Aph(:,:)
  double precision,allocatable  :: Bph(:,:)
  double precision,allocatable  :: AB(:,:)
  double precision,allocatable  :: Om(:)
  double precision,allocatable  :: R(:,:)
  double precision,allocatable  :: ExpR(:,:)
  
  integer,allocatable          :: virtuals(:)
  integer,allocatable          :: occupations_fc(:)
  
  integer                       :: ixyz
  integer                       :: eig

! Output variables

  double precision,intent(out)  :: ERHF
  double precision,intent(out)  :: e(nOrb)
  double precision,intent(inout):: c(nBas,nOrb)
  double precision,intent(out)  :: P(nBas,nBas)
  double precision,intent(out)  :: F(nBas,nBas)

! Memory allocation

  write(*,*)
  write(*,*) '****************************'
  write(*,*) '* Search for RHF solutions *'
  write(*,*) '****************************'
  write(*,*)

! CVS
  
  print *, "No exications to the first", nCVS, "orbital(s) are considered."
  if(nC/=0) then
    print *, "Do not use PyDuck frozen core with CVS !"
    stop
  end if

! Frozen Core

  nFC = MERGE(1,0,FC/=0) 
  allocate(occupations_fc(nO-nFC))
  ! remove FC from occupations
  do ispin=1,nspin
    occupations_fc(1:nO-nFC) = occupations(1:nO - nFC) 
    found = .false.
    do i=1,nO-1
      if(.not. found) then
        if(occupations(i)==FC) then
          found = .true.
          occupations_fc(i) = occupations(i+1) 
        else
          occupations_fc(i) = occupations(i)
        endif
      else
        occupations_fc(i) = occupations(i+1) 
      endif 
    enddo
  enddo
  print *, "Not Frozen orbitals:"
  print *,occupations_fc(1:nO-nFC)
  
!-------------------!
! Memory allocation !
!-------------------!

  nSCVS = (nV - nCVS)*(nO - nFC)
  
  allocate(Aph(nSCVS,nSCVS),Bph(nSCVS,nSCVS),AB(nSCVS,nSCVS),Om(nSCVS))
  allocate(R(nOrb,nOrb),ExpR(nOrb,nOrb),virtuals(nV))
  call non_occupied(nO,nBas,occupations,virtuals)

!------------------!
! Search algorithm !
!------------------!

  unstab = .true.
  guess  = 0

  do while(unstab) 

!---------------------!
! Hartree-Fock module !
!---------------------!

    call wall_time(start_HF)
    call MOM_RHF(.false.,doaordm,maxSCF,thresh,max_diis,guess,level_shift,writeMOs,nNuc,ZNuc,rNuc,ENuc,&
             nBas,nOrb,nO,S,T,V,Hc,ERI_AO,dipole_int_AO,X,ERHF,e,c,P,F,occupations)
    call wall_time(end_HF)

    t_HF = end_HF - start_HF
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for RHF = ',t_HF,' seconds'
    write(*,*)

!----------------------------------!
! AO to MO integral transformation !
!----------------------------------!

    call wall_time(start_AOtoMO)
    write(*,*)
    write(*,*) 'AO to MO transformation... Please be patient'
    write(*,*)
    do ixyz = 1,ncart
      call AOtoMO(nBas,nOrb,c,dipole_int_AO(1,1,ixyz),dipole_int_MO(1,1,ixyz))
    end do
    call AOtoMO_ERI_RHF(nBas,nOrb,c,ERI_AO,ERI_MO)
    call wall_time(end_AOtoMO)
 
    t_AOtoMO = end_AOtoMO - start_AOtoMO
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for AO to MO transformation = ',t_AOtoMO,' seconds'
    write(*,*)

!-------------------------------------------------------------!
! Stability analysis: Real RHF -> Real RHF  
!-------------------------------------------------------------!

    ispin = 1
    
    call CVS_phRLR_A(ispin,.false.,nOrb,nC,nO,nV,nR,nSCVS,nCVS,nFC,occupations_fc,virtuals,1d0,e,ERI_MO,Aph)
    call CVS_phRLR_B(ispin,.false.,nOrb,nC,nO,nV,nR,nSCVS,nCVS,nFC,occupations_fc,virtuals,1d0,ERI_MO,Bph)
 
    AB(:,:) = Aph(:,:) + Bph(:,:)
 
    call diagonalize_matrix(nSCVS,AB,Om)
    Om(:) = 2d0*Om(:)

    write(*,*)'-------------------------------------------------------------'
    write(*,*)'|       Stability analysis: Real RHF -> Real RHF            |'
    write(*,*)'-------------------------------------------------------------'
    write(*,'(1X,A1,1X,A5,1X,A1,1X,A23,1X,A1,1X,A23,1X,A1,1X)') &
              '|','State','|',' Excitation energy (au) ','|',' Excitation energy (eV) ','|'
    write(*,*)'-------------------------------------------------------------'
    do ia=1,min(nSCVS,maxS)
      write(*,'(1X,A1,1X,I5,1X,A1,1X,F23.6,1X,A1,1X,F23.6,1X,A1,1X)') &
        '|',ia,'|',Om(ia),'|',Om(ia)*HaToeV,'|'
    end do
    write(*,*)'-------------------------------------------------------------'
 

    write(*,'(1X,A40,1X,F15.10,A3)') 'Smallest eigenvalue:',Om(1),' au'
    write(*,'(1X,A40,1X,F15.10,A3)') 'E(RHF) = ',ENuc + ERHF,' au'
    write(*,*) 
    write(*,'(1X,A40,1X,A10)')       'Which one would you like to follow?','[Exit:0]'
    read(*,*) eig

    if(eig < 0 .or. eig > nSCVS)  then
      write(*,'(1X,A40,1X,A10)')     'Invalid option...','Stop...'
      write(*,*)
      deallocate(Aph,Bph,AB,Om,R,ExpR)
      stop
    end if

    if(eig == 0) return

    R(:,:) = 0d0
    ia = 0
    do i=1,nO-nFC
      do a=1+nCVS,nV
        ia = ia + 1
        R(virtuals(a),occupations_fc(i)) = +AB(ia,eig)
        R(occupations_fc(i),virtuals(a)) = -AB(ia,eig)
      end do
    end do

    call matrix_exponential(nOrb,R,ExpR)
    c = matmul(c,ExpR)
 
    ! Do Frozen core
    do ispin=1,nspin
      occupations_fc(1:nO-nFC) = occupations(1:nO - nFC) 
      found = .false.
      do i=1,nO-1
        if(.not. found) then
          if(occupations(i)==FC) then
            found = .true.
            occupations_fc(i) = occupations(i+1) 
          else
            occupations_fc(i) = occupations(i)
          endif
        else
          occupations_fc(i) = occupations(i+1) 
        endif 
      enddo
    enddo

    write(*,*)'-------------------------------------------------------------'
    write(*,*)

!---------------!
! End of Search !
!---------------!
  end do

  deallocate(Aph,Bph,AB,Om,R,ExpR)

end subroutine
