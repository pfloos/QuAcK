subroutine UHF_search(maxSCF,thresh,max_diis,guess_type,mix,level_shift,nNuc,ZNuc,rNuc,ENuc, &
                      nBas,nC,nO,nV,nR,S,T,V,Hc,ERI_AO,ERI_aaaa,ERI_aabb,ERI_bbbb,           &
                      dipole_int_AO,dipole_int_aa,dipole_int_bb,X,EUHF,e,c,P)

! Search for UHF solutions

  implicit none
  include 'parameters.h'

  integer,intent(in)            :: maxSCF
  integer,intent(in)            :: max_diis
  integer,intent(in)            :: guess_type
  double precision,intent(in)   :: thresh
  double precision,intent(inout):: mix
  double precision,intent(in)   :: level_shift

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC(nspin)
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nV(nspin)
  integer,intent(in)            :: nR(nspin)
  integer,intent(in)            :: nNuc
  double precision,intent(in)   :: ZNuc(nNuc)
  double precision,intent(in)   :: rNuc(nNuc,ncart)
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: T(nBas,nBas)
  double precision,intent(in)   :: V(nBas,nBas)
  double precision,intent(in)   :: Hc(nBas,nBas)
  double precision,intent(in)   :: X(nBas,nBas)
  double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas)
  double precision,intent(inout):: ERI_aaaa(nBas,nBas,nBas,nBas)
  double precision,intent(inout):: ERI_aabb(nBas,nBas,nBas,nBas)
  double precision,intent(inout):: ERI_bbbb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int_AO(nBas,nBas,ncart)
  double precision,intent(inout):: dipole_int_aa(nBas,nBas,ncart)
  double precision,intent(inout):: dipole_int_bb(nBas,nBas,ncart)

! Local variables

  double precision              :: start_HF     ,end_HF       ,t_HF
  double precision              :: start_stab   ,end_stab     ,t_stab
  double precision              :: start_AOtoMO ,end_AOtoMO   ,t_AOtoMO

  logical                       :: unstab
  integer                       :: guess

  integer                       :: nS(nspin)

  integer,parameter             :: maxS = 20
  integer                       :: ia,i,a,mu
  integer                       :: ispin

  integer                       :: nS_aa,nS_bb,nS_sc
  double precision,allocatable  :: Om_sc(:)
  double precision,allocatable  :: A_sc(:,:)
  double precision,allocatable  :: B_sc(:,:)
  double precision,allocatable  :: AB_sc(:,:)
  double precision,allocatable  :: R(:,:)
  double precision,allocatable  :: ExpR(:,:)
  
  integer                       :: ixyz
  integer                       :: eig

! Output variables

  double precision,intent(out)  :: EUHF
  double precision,intent(out)  :: e(nBas,nspin)
  double precision,intent(inout):: c(nBas,nBas,nspin)
  double precision,intent(out)  :: P(nBas,nBas,nspin)

! Memory allocation

  write(*,*)
  write(*,*) '****************************'
  write(*,*) '* Search for UHF solutions *'
  write(*,*) '****************************'
  write(*,*)

!-------------------!
! Memory allocation !
!-------------------!

  nS(:) = (nO(:) - nC(:))*(nV(:) - nR(:))

  nS_aa = nS(1)
  nS_bb = nS(2)
  nS_sc = nS_aa + nS_bb

  allocate(Om_sc(nS_sc),A_sc(nS_sc,nS_sc),B_sc(nS_sc,nS_sc),AB_sc(nS_sc,nS_sc),R(nBas,nBas),ExpR(nBas,nBas))

!------------------!
! Search algorithm !
!------------------!

  unstab = .true.
  guess  = 0
  mix = 0d0

  do while(unstab) 

!---------------------!
! Hartree-Fock module !
!---------------------!

    call wall_time(start_HF)
    call UHF(.false.,maxSCF,thresh,max_diis,guess,mix,level_shift,nNuc,ZNuc,rNuc,ENuc, &
             nBas,nO,S,T,V,Hc,ERI_AO,dipole_int_AO,X,EUHF,e,c,P)
    call wall_time(end_HF)

    t_HF = end_HF - start_HF
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for UHF = ',t_HF,' seconds'
    write(*,*)

!----------------------------------!
! AO to MO integral transformation !
!----------------------------------!

    call wall_time(start_AOtoMO)
    write(*,*)
    write(*,*) 'AO to MO transformation... Please be patient'
    write(*,*)

  ! Transform dipole-related integrals

    do ixyz=1,ncart
      call AOtoMO(nBas,nBas,c(:,:,1),dipole_int_AO(:,:,ixyz),dipole_int_aa(:,:,ixyz))
      call AOtoMO(nBas,nBas,c(:,:,2),dipole_int_AO(:,:,ixyz),dipole_int_bb(:,:,ixyz))
    end do

    ! 4-index transform for (aa|aa) block
    call AOtoMO_ERI_UHF(1,1,nBas,c,ERI_AO,ERI_aaaa)
 
    ! 4-index transform for (aa|bb) block
    call AOtoMO_ERI_UHF(1,2,nBas,c,ERI_AO,ERI_aabb)
 
    ! 4-index transform for (bb|bb) block
    call AOtoMO_ERI_UHF(2,2,nBas,c,ERI_AO,ERI_bbbb)

    call wall_time(end_AOtoMO)
 
    t_AOtoMO = end_AOtoMO - start_AOtoMO
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for AO to MO transformation = ',t_AOtoMO,' seconds'
    write(*,*)

!-------------------------------------------------------------!
! Stability analysis: Real UHF -> Real UHF  
!-------------------------------------------------------------!

    ispin = 1
 
    call phULR_A(ispin,.false.,nBas,nC,nO,nV,nR,nS_aa,nS_bb,nS_sc,1d0,e,ERI_aaaa,ERI_aabb,ERI_bbbb,A_sc)
    call phULR_B(ispin,.false.,nBas,nC,nO,nV,nR,nS_aa,nS_bb,nS_sc,1d0,ERI_aaaa,ERI_aabb,ERI_bbbb,B_sc)
 
    AB_sc(:,:) = A_sc(:,:) + B_sc(:,:)
 
    call diagonalize_matrix(nS_sc,AB_sc,Om_sc)
    Om_sc(:) = 2d0*Om_sc(:)
 
    write(*,*)'-------------------------------------------------------------'
    write(*,*)'|       Stability analysis: Real UHF -> Real UHF            |'
    write(*,*)'-------------------------------------------------------------'
    write(*,'(1X,A1,1X,A5,1X,A1,1X,A23,1X,A1,1X,A23,1X,A1,1X)') &
              '|','State','|',' Excitation energy (au) ','|',' Excitation energy (eV) ','|'
    write(*,*)'-------------------------------------------------------------'
    do ia=1,min(nS_sc,maxS)
      write(*,'(1X,A1,1X,I5,1X,A1,1X,F23.6,1X,A1,1X,F23.6,1X,A1,1X)') &
        '|',ia,'|',Om_sc(ia),'|',Om_sc(ia)*HaToeV,'|'
    end do
    write(*,*)'-------------------------------------------------------------'
 
    if(minval(Om_sc(:)) < 1d-7) then
 
      write(*,'(1X,A40,1X)')           'Too bad, UHF solution is unstable!'
      write(*,'(1X,A40,1X,F15.10,A3)') 'Largest negative eigenvalue:',Om_sc(1),' au'
      write(*,'(1X,A40,1X,F15.10,A3)') 'E(UHF) = ',ENuc + EUHF,' au'
      write(*,*)
      write(*,'(1X,A40,1X,A10)')       'Which one would you like to follow?','[Exit:0]'
      read(*,*) eig

      if(eig < 0 .or. eig > nS_sc)  then
        write(*,'(1X,A40,1X,A10)')     'Invalid option...','Stop...'
        write(*,*)
        stop
      end if

      if(eig == 0) return

      ! Spin-up kick

      R(:,:) = 0d0
      ia = 0
      do i=nC(1)+1,nO(1)
        do a=nO(1)+1,nBas-nR(1)
          ia = ia + 1
          R(a,i) = +AB_sc(ia,eig)
          R(i,a) = -AB_sc(ia,eig)
        end do
      end do

      call matrix_exponential(nBas,R,ExpR)
      c(:,:,1) = matmul(c(:,:,1),ExpR)

      ! Spin-down kick

      R(:,:) = 0d0
      ia = nS_aa
      do i=nC(2)+1,nO(2)
        do a=nO(2)+1,nBas-nR(2)
          ia = ia + 1
          R(a,i) = +AB_sc(ia,eig)
          R(i,a) = -AB_sc(ia,eig)
        end do
      end do

      call matrix_exponential(nBas,R,ExpR)
      c(:,:,2) = matmul(c(:,:,2),ExpR)
 
    else
 
      write(*,'(1X,A40,1X)')           'Well done, UHF solution is stable!'
      write(*,'(1X,A40,1X,F15.10,A3)') 'Smallest eigenvalue: ',Om_sc(1),' au'
      write(*,'(1X,A40,1X,F15.10,A3)') 'E(UHF) = ',ENuc + EUHF,' au'
 
      unstab = .false.
 
    end if
    write(*,*)'-------------------------------------------------------------'
    write(*,*)

!---------------!
! End of Search !
!---------------!
  end do

end subroutine
