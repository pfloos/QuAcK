subroutine RHF_search(maxSCF,thresh,max_diis,guess_type,level_shift,nNuc,ZNuc,rNuc,ENuc, &
                      nBas,nC,nO,nV,nR,S,T,V,Hc,ERI_AO,dipole_int_AO,X,EHF,e,c,P)

! Restricted branch of QuAcK

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
  double precision,intent(in)   :: dipole_int_AO(nBas,nBas,ncart)

! Local variables

  double precision              :: start_HF     ,end_HF       ,t_HF
  double precision              :: start_stab   ,end_stab     ,t_stab
  double precision              :: start_AOtoMO ,end_AOtoMO   ,t_AOtoMO

  logical                       :: unstab
  integer                       :: guess
  double precision,allocatable  :: ERI_MO(:,:,:,:)
  integer                       :: nS

  integer,parameter             :: maxS = 20
  integer                       :: ia,i,a,mu
  integer                       :: ispin

  double precision,allocatable  :: Aph(:,:)
  double precision,allocatable  :: Bph(:,:)
  double precision,allocatable  :: AB(:,:)
  double precision,allocatable  :: Om(:)
  double precision,allocatable  :: R(:,:)
  
  integer                       :: eig
  double precision              :: kick,step

! Output variables

  double precision,intent(out)  :: EHF
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
  allocate(ERI_MO(nBas,nBas,nBas,nBas),Aph(nS,nS),Bph(nS,nS),AB(nS,nS),Om(nS),R(nBas,nBas))

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
    call RHF(maxSCF,thresh,max_diis,guess,level_shift,nNuc,ZNuc,rNuc,ENuc, &
             nBas,nO,S,T,V,Hc,ERI_AO,dipole_int_AO,X,EHF,e,c,P)
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
    call AOtoMO_integral_transform(1,1,1,1,nBas,c,ERI_AO,ERI_MO)
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
    Om(:) = 0.5d0*Om(:)
 
    write(*,*)'-------------------------------------------------------------'
    write(*,*)'|       Stability analysis: Real RHF -> Real RHF            |'
    write(*,*)'-------------------------------------------------------------'
    write(*,'(1X,A1,1X,A5,1X,A1,1X,A23,1X,A1,1X,A23,1X,A1,1X)') &
              '|','State','|',' Excitation energy (au) ','|',' Excitation energy (eV) ','|'
    write(*,*)'-------------------------------------------------------------'
    do ia=1,min(nS,maxS)
      write(*,'(1X,A1,1X,I5,1X,A1,1X,F23.6,1X,A1,1X,F23.6,1X,A1,1X)') &
        '|',ia,'|',Om(ia),'|',Om(ia)*HaToeV,'|'
    enddo
    write(*,*)'-------------------------------------------------------------'
 
    if(minval(Om(:)) < 0d0) then
 
      write(*,'(1X,A40,1X)')        'Too bad, RHF solution is unstable!'
      write(*,'(1X,A40,1X,F15.10,A3)') 'Largest negative eigenvalue:',Om(1),' au'
      write(*,*) 
      write(*,'(1X,A40,1X)')        'Which one would you like to follow?'
      read(*,*) eig

      step = 1d0

      do mu=1,nBas
        ia = 0
        do i=nC+1,nO
          kick = 0d0
          do a=nO+1,nBas-nR
            ia = ia + 1
            kick = kick + AB(ia,eig)*c(mu,a)
          end do
          c(mu,i) = c(mu,i) + step*kick
        end do
      end do

    else
 
      write(*,'(1X,A40,1X)')        'Well done, RHF solution is stable!'
      write(*,'(1X,A40,1X,F15.10,A3)') 'Smallest eigenvalue: ',Om(1),' au'
 
      unstab = .false.
 
    end if
    write(*,*)'-------------------------------------------------------------'
    write(*,*)

  end do

end subroutine
