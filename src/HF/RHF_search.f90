subroutine RHF_search(maxSCF,thresh,max_diis,guess_type,level_shift,nNuc,ZNuc,rNuc,ENuc,   &
                      nBas,nC,nO,nV,nR,S,T,V,Hc,ERI_AO,ERI_MO,dipole_int_AO,dipole_int_MO, & 
                      X,ERHF,e,c,P)

! Search for RHF solutions

  implicit none
  include 'parameters.h'

  integer,intent(in)            :: maxSCF
  integer,intent(in)            :: max_diis
  integer,intent(in)            :: guess_type
  double precision,intent(in)   :: thresh
  double precision,intent(in)   :: level_shift

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
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
  double precision,intent(inout):: ERI_MO(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int_AO(nBas,nBas,ncart)
  double precision,intent(inout):: dipole_int_MO(nBas,nBas,ncart)

! Local variables

  double precision              :: start_HF     ,end_HF       ,t_HF
  double precision              :: start_stab   ,end_stab     ,t_stab
  double precision              :: start_AOtoMO ,end_AOtoMO   ,t_AOtoMO

  logical                       :: unstab
  integer                       :: guess
  integer                       :: nS

  integer,parameter             :: maxS = 20
  integer                       :: ia,i,a,mu
  integer                       :: ispin

  double precision,allocatable  :: Aph(:,:)
  double precision,allocatable  :: Bph(:,:)
  double precision,allocatable  :: AB(:,:)
  double precision,allocatable  :: Om(:)
  double precision,allocatable  :: R(:,:)
  double precision,allocatable  :: ExpR(:,:)
  
  integer                       :: ixyz
  integer                       :: eig

! Output variables

  double precision,intent(out)  :: ERHF
  double precision,intent(out)  :: e(nBas)
  double precision,intent(inout):: c(nBas,nBas)
  double precision,intent(out)  :: P(nBas,nBas)

! Memory allocation

  write(*,*)
  write(*,*) '****************************'
  write(*,*) '* Search for RHF solutions *'
  write(*,*) '****************************'
  write(*,*)

!-------------------!
! Memory allocation !
!-------------------!

  nS = (nO - nC)*(nV - nR)
  allocate(Aph(nS,nS),Bph(nS,nS),AB(nS,nS),Om(nS),R(nBas,nBas),ExpR(nBas,nBas))

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
    call RHF(.false.,maxSCF,thresh,max_diis,guess,level_shift,nNuc,ZNuc,rNuc,ENuc, &
             nBas,nO,S,T,V,Hc,ERI_AO,dipole_int_AO,X,ERHF,e,c,P)
    call wall_time(end_HF)

    t_HF = end_HF - start_HF
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for RHF = ',t_HF,' seconds'
    write(*,*)

!----------------------------------!
! AO to MO integral transformation !
!----------------------------------!

    call wall_time(start_AOtoMO)
    write(*,*)
    write(*,*) 'AO to MO transformation... Please be patient'
    write(*,*)
    do ixyz=1,ncart
      call AOtoMO(nBas,c,dipole_int_AO(:,:,ixyz),dipole_int_MO(:,:,ixyz))
    end do
    call AOtoMO_ERI_RHF(nBas,c,ERI_AO,ERI_MO)
    call wall_time(end_AOtoMO)
 
    t_AOtoMO = end_AOtoMO - start_AOtoMO
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for AO to MO transformation = ',t_AOtoMO,' seconds'
    write(*,*)

!-------------------------------------------------------------!
! Stability analysis: Real RHF -> Real RHF  
!-------------------------------------------------------------!

    ispin = 1
 
    call phLR_A(ispin,.false.,nBas,nC,nO,nV,nR,nS,1d0,e,ERI_MO,Aph)
    call phLR_B(ispin,.false.,nBas,nC,nO,nV,nR,nS,1d0,ERI_MO,Bph)
 
    AB(:,:) = Aph(:,:) + Bph(:,:)
 
    call diagonalize_matrix(nS,AB,Om)
    Om(:) = 2d0*Om(:)

    write(*,*)'-------------------------------------------------------------'
    write(*,*)'|       Stability analysis: Real RHF -> Real RHF            |'
    write(*,*)'-------------------------------------------------------------'
    write(*,'(1X,A1,1X,A5,1X,A1,1X,A23,1X,A1,1X,A23,1X,A1,1X)') &
              '|','State','|',' Excitation energy (au) ','|',' Excitation energy (eV) ','|'
    write(*,*)'-------------------------------------------------------------'
    do ia=1,min(nS,maxS)
      write(*,'(1X,A1,1X,I5,1X,A1,1X,F23.6,1X,A1,1X,F23.6,1X,A1,1X)') &
        '|',ia,'|',Om(ia),'|',Om(ia)*HaToeV,'|'
    end do
    write(*,*)'-------------------------------------------------------------'
 
    if(minval(Om(:)) < 1d-7) then
 
      write(*,'(1X,A40,1X)')           'Too bad, RHF solution is unstable!'
      write(*,'(1X,A40,1X,F15.10,A3)') 'Largest negative eigenvalue:',Om(1),' au'
      write(*,'(1X,A40,1X,F15.10,A3)') 'E(RHF) = ',ENuc + ERHF,' au'
      write(*,*) 
      write(*,'(1X,A40,1X,A10)')       'Which one would you like to follow?','[Exit:0]'
      read(*,*) eig

      if(eig < 0 .or. eig > nS)  then
        write(*,'(1X,A40,1X,A10)')     'Invalid option...','Stop...'
        write(*,*)
        stop
      end if

      if(eig == 0) return

      R(:,:) = 0d0
      ia = 0
      do i=nC+1,nO
        do a=nO+1,nBas-nR
          ia = ia + 1
          R(a,i) = +AB(ia,eig)
          R(i,a) = -AB(ia,eig)
        end do
      end do

      call matrix_exponential(nBas,R,ExpR)
      c = matmul(c,ExpR)

    else
 
      write(*,'(1X,A40,1X)')           'Well done, RHF solution is stable!'
      write(*,'(1X,A40,1X,F15.10,A3)') 'Smallest eigenvalue: ',Om(1),' au'
      write(*,'(1X,A40,1X,F15.10,A3)') 'E(RHF) = ',ENuc + ERHF,' au'
 
      unstab = .false.
 
    end if
    write(*,*)'-------------------------------------------------------------'
    write(*,*)

!---------------!
! End of Search !
!---------------!
  end do

end subroutine
