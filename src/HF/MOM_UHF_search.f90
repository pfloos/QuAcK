subroutine MOM_UHF_search(maxSCF,thresh,max_diis,guess_type,mix,level_shift,writeMOs,nNuc,ZNuc,rNuc,ENuc, &
                          nBas,nC,nO,nV,nR,nCVS,FC,S,T,V,Hc,ERI_AO,ERI_aaaa,ERI_aabb,ERI_bbbb,           &
                          dipole_int_AO,dipole_int_aa,dipole_int_bb,X,EUHF,e,c,P,F,occupations,working_dir)

! Search for MOM-UHF solutions
  
  implicit none
  include 'parameters.h'

  integer,intent(in)            :: maxSCF
  integer,intent(in)            :: max_diis
  integer,intent(in)            :: guess_type
  double precision,intent(in)   :: thresh
  double precision,intent(inout):: mix
  double precision,intent(in)   :: level_shift
  logical,intent(in)            :: writeMOs
  
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC(nspin)
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nV(nspin)
  integer,intent(in)            :: nR(nspin)
  integer,intent(in)            :: nCVS(nspin)
  integer,intent(in)            :: FC(nspin)
  integer,intent(in)            :: nNuc
  integer,intent(in)            :: occupations(maxval(nO),nspin)
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
  
  character(len=256),intent(in) :: working_dir

! Local variables

  double precision              :: start_HF     ,end_HF       ,t_HF
  double precision              :: start_stab   ,end_stab     ,t_stab
  double precision              :: start_AOtoMO ,end_AOtoMO   ,t_AOtoMO
  double precision              :: kick

  logical                       :: unstab,found
  integer                       :: guess

  integer                       :: nS(nspin),nFC(nspin)

  integer,parameter             :: maxS = 20
  integer                       :: ia,i,a,mu
  integer                       :: ispin

  integer                       :: nSa,nSb,nSt
  double precision,allocatable  :: Om(:)
  double precision,allocatable  :: Aph(:,:)
  double precision,allocatable  :: Bph(:,:)
  double precision,allocatable  :: AB(:,:)
  double precision,allocatable  :: R(:,:)
  double precision,allocatable  :: ExpR(:,:)
  double precision              :: thresh_eig = 0.5d0
  
  integer                       :: ixyz
  integer                       :: eig
  integer,allocatable           :: occupations_fc(:,:)
  integer,allocatable           :: virtuals(:,:)

! Output variables

  double precision,intent(out)  :: EUHF
  double precision,intent(out)  :: e(nBas,nspin)
  double precision,intent(inout):: c(nBas,nBas,nspin)
  double precision,intent(out)  :: P(nBas,nBas,nspin)
  double precision,intent(out)  :: F(nBas,nBas,nspin)





! Memory allocation

  write(*,*)
  write(*,*) '********************************'
  write(*,*) '* Search for MOM-UHF solutions *'
  write(*,*) '********************************'
  write(*,*)

! Memory allocation

  nFC(1) = MERGE(1,0,FC(1)/=0) 
  nFC(2) = MERGE(1,0,FC(2)/=0)
  nSa = (nBas - nO(1) - nCVS(1))*(nO(1) - nFC(1))
  nSb = (nBas - nO(2) - nCVS(2))*(nO(2) - nFC(2))
  nSt = nSa + nSb
  allocate(Om(nSt),Aph(nSt,nSt),Bph(nSt,nSt),AB(nSt,nSt),R(nBas,nBas),ExpR(nBas,nBas))
  allocate(virtuals(nBas - minval(nO),nspin))
  virtuals = 0
  do ispin=1,nspin
    call non_occupied(nO(ispin),nBas,occupations(1:nO(ispin),ispin),virtuals(1:nBas - nO(ispin),ispin))
  end do


  ! CVS

  print *, "No exications to the first", nCVS(1), "alpha orbital(s) are considered."
  print *, "No exications to the first", nCVS(2), "beta orbital(s) are considered."
  if(any(nC/=0)) then
    print *, "Do not use PyDuck frozen core with CVS !"
    stop
  end if

! Frozen Core
  
  allocate(occupations_fc(maxval(nO-nFC),nspin))
  ! remove FC from occupations
  do ispin=1,nspin
    occupations_fc(1:nO(ispin)-nFC(ispin),ispin) = occupations(1:nO(ispin) - nFC(ispin), ispin) 
    found = .false.
    do i=1,nO(ispin)-1
      if(.not. found) then
        if(occupations(i,ispin)==FC(ispin)) then
          found = .true.
          occupations_fc(i,ispin) = occupations(i+1,ispin) 
        else
          occupations_fc(i,ispin) = occupations(i,ispin)
        endif
      else
        occupations_fc(i,ispin) = occupations(i+1,ispin) 
      endif 
    enddo
  enddo
  do ispin=1,nspin
    print *, "Not Frozen orbitals:"
    print *,occupations_fc(1:nO(ispin)-nFC(ispin),ispin)
  end do

!------------------!
! Search algorithm !
!------------------!

  unstab = .true.
  guess  = 0
  mix    = 0d0

  do while(unstab) 

!---------------------!
! Hartree-Fock module !
!---------------------!

    call wall_time(start_HF)
    call MOM_UHF(.false.,maxSCF,thresh,max_diis,guess,mix,level_shift,writeMOs,nNuc,ZNuc,rNuc,ENuc, & 
             nBas,nO,S,T,V,Hc,ERI_AO,dipole_int_AO,X,EUHF,e,c,P,F,occupations,working_dir)
    call wall_time(end_HF)

    t_HF = end_HF - start_HF
    write(*,'(A65,1X,F9.3,A8)') 'Total wall time for UHF = ',t_HF,' seconds'
    write(*,*)

    ! Do Frozen core
    do ispin=1,nspin
      occupations_fc(1:nO(ispin)-nFC(ispin),ispin) = occupations(1:nO(ispin) - nFC(ispin), ispin) 
      found = .false.
      do i=1,nO(ispin)-1
        if(.not. found) then
          if(occupations(i,ispin)==FC(ispin)) then
            found = .true.
            occupations_fc(i,ispin) = occupations(i+1,ispin) 
          else
            occupations_fc(i,ispin) = occupations(i,ispin)
          endif
        else
          occupations_fc(i,ispin) = occupations(i+1,ispin) 
        endif 
      enddo
    enddo

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
    
    call MOM_phULR_A(ispin,.false.,nBas,nC,nO,nV,nR,nSa,nSb,nSt,nCVS,nFC,occupations_fc,virtuals,1d0,e,ERI_aaaa,ERI_aabb,ERI_bbbb,Aph)
    call MOM_phULR_B(ispin,.false.,nBas,nC,nO,nV,nR,nSa,nSb,nSt,nCVS,nFC,occupations_fc,virtuals,1d0,ERI_aaaa,ERI_aabb,ERI_bbbb,Bph)
 
    AB(:,:) = Aph(:,:) + Bph(:,:)
 
    call diagonalize_matrix(nSt,AB,Om)
    Om(:) = 2d0*Om(:)

    write(*,*)'-------------------------------------------------------------'
    write(*,*)'|       Stability analysis: Real UHF -> Real UHF            |'
    write(*,*)'-------------------------------------------------------------'
    write(*,'(1X,A1,1X,A5,1X,A1,1X,A23,1X,A1,1X,A23,1X,A1,1X)') &
              '|','State','|',' Excitation energy (au) ','|',' Excitation energy (eV) ','|'
    write(*,*)'-------------------------------------------------------------'
    do ia=1,min(nSt,maxS)
      write(*,'(1X,A1,1X,I5,1X,A1,1X,F23.6,1X,A1,1X,F23.6,1X,A1,1X)') &
        '|',ia,'|',Om(ia),'|',Om(ia)*HaToeV,'|'
    end do
    write(*,*)'-------------------------------------------------------------'
    do mu = 1,min(nSt,maxS)
      write(*,*) 'Eigenvale no:',mu
      do i=1,nO(1) - nFC(1)
        do a=nCVS(1)+1,nBas - nO(1)
          ia = ia + 1
          if (abs(AB(ia,mu)) > thresh_eig) then
            write(*,*) 'Transition alpha', ':' ,'weight'
            write(*,*) occupations_fc(i,1),'->',virtuals(a,1),':',abs(AB(ia,mu))
          end if
        end do
      end do
      ia = nSa
      do i=1,nO(2) - nFC(2)
        do a=nCVS(2)+1,nBas - nO(2)
          ia = ia + 1
          if (abs(AB(ia,mu)) > thresh_eig) then
            write(*,*) 'Transition beta', ':' ,'weight'
            write(*,*) occupations_fc(i,2),'->',virtuals(a,2),':',abs(AB(ia,mu))
          end if
        end do
      end do

      write(*,*) '--------------------------------'

    end do 
    write(*,'(1X,A40,1X,F15.10,A3)') 'Smallest eigenvalue:',Om(1),' au'
    write(*,'(1X,A40,1X,F15.10,A3)') 'E(UHF) = ',ENuc + EUHF,' au'
    write(*,*)
    write(*,'(1X,A40,1X,A10)')       'Which one would you like to follow?','[Exit:0]'
    read(*,*) eig
    if(eig == 0) return
    write(*,'(1X,A40,1X,A10)')       'How strong do you want to kick?', '[Exit: 0]'
    read(*,*) kick
    if(kick == 0) return
    

    if(eig < 0 .or. eig > nSt)  then
      write(*,'(1X,A40,1X,A10)')     'Invalid option...','Stop...'
      write(*,*)
      deallocate(Aph,Bph,AB,Om,R,ExpR)
      stop
    end if


    ! Spin-up kick

    R(:,:) = 0d0
    ia = 0
    do i=1,nO(1) - nFC(1)
      do a=nCVS(1)+1,nBas - nO(1)
        ia = ia + 1
        R(virtuals(a,1),occupations_fc(i,1)) = +AB(ia,eig)*kick
        R(occupations_fc(i,1),virtuals(a,1)) = -AB(ia,eig)*kick
      end do
    end do

    call matrix_exponential(nBas,R,ExpR)
    c(:,:,1) = matmul(c(:,:,1),ExpR)

    ! Spin-down kick

    R(:,:) = 0d0
    ia = nSa
    do i=1,nO(2) - nFC(2)
      do a=nCVS(2)+1,nBas - nO(2)
        ia = ia + 1
        R(virtuals(a,2),occupations_fc(i,2)) = +AB(ia,eig)*kick
        R(occupations_fc(i,2),virtuals(a,2)) = -AB(ia,eig)*kick
      end do
    end do

    call matrix_exponential(nBas,R,ExpR)
    c(:,:,2) = matmul(c(:,:,2),ExpR)

    write(*,*)'-------------------------------------------------------------'
    write(*,*)
    write(*,*) unstab
!---------------!
! End of Search !
!---------------!
  end do

  deallocate(Aph,Bph,AB,Om,R,ExpR)

end subroutine
