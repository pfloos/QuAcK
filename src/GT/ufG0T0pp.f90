subroutine ufG0T0pp(dotest,TDA_W,nBas,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,eHF)

! Upfolded G0T0pp equations

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dotest

  logical,intent(in)            :: TDA_W
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: ERHF
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: eHF(nBas)

! Local variables

  integer                       :: p
  integer                       :: s
  integer                       :: i,j,k,l
  integer                       :: a,b,c,d
  integer                       :: jb,kc,ia,ja
  integer                       :: klc,kcd,ija,ijb,iab,jab

  logical                       :: dRPA
  integer                       :: ispin
  double precision              :: EcRPA
  integer                       :: n2h1p,n2p1h,nH
  double precision,external     :: Kronecker_delta
  double precision,allocatable  :: H(:,:)
  double precision,allocatable  :: eGW(:)
  double precision,allocatable  :: Z(:)
  double precision,allocatable  :: Aph(:,:)
  double precision,allocatable  :: Bph(:,:)
  double precision,allocatable  :: Om(:)
  double precision,allocatable  :: XpY(:,:)
  double precision,allocatable  :: XmY(:,:)
  double precision,allocatable  :: rho(:,:,:)

  logical                       :: verbose = .true.
  double precision,parameter    :: cutoff1 = 0.01d0
  double precision,parameter    :: cutoff2 = 0.01d0
  double precision              :: eF
  double precision,parameter    :: window = 2d0

  double precision              :: start_timing,end_timing,timing

! Output variables

! Hello world

  write(*,*)
  write(*,*)'******************************************'
  write(*,*)'* Restricted Upfolded G0T0pp Calculation *'
  write(*,*)'******************************************'
  write(*,*)

! Dimension of the supermatrix

  n2h1p = nO*nO*nV
  n2p1h = nV*nV*nO
  nH = 1 + n2h1p + n2p1h

! Memory allocation

  allocate(H(nH,nH),eGW(nH),Z(nH))
  
! Initialization

  dRPA = .true.
  EcRPA = 0d0

  eF = 0.5d0*(eHF(nO+1) + eHF(nO))

!-------------------------!
! Main loop over orbitals !
!-------------------------!

  do p=nO-1,nO

    H(:,:) = 0d0

    if (TDA_W) then
 
      ! TDA for W
 
      write(*,*) 'Tamm-Dancoff approximation actived!'
      write(*,*)
 
      !---------------------------!
      !  Compute GW supermatrix   !
      !---------------------------!
      !                           !
      !     |   F   V2h1p V2p1h | ! 
      !     |                   | ! 
      ! H = | V2h1p C2h1p     0 | ! 
      !     |                   | ! 
      !     | V2p1h   0   C2p1h | ! 
      !                           !
      !---------------------------!
 
      call wall_time(start_timing)

      !---------!
      ! Block F !
      !---------!
      
!     H(1,1) = eHF(p)
!
!     !-------------!
!     ! Block V2h1p !
!     !-------------!

!     ija = 0
!     do i=nC+1,nO
!       do j=nC+1,nO
!         do a=nO+1,nBas-nR
!           ija = ija + 1
!              
!           H(1    ,1+ija) = sqrt(2d0)*ERI(p,a,i,j)
!           H(1+ija,1    ) = sqrt(2d0)*ERI(p,a,i,j)
!
!         end do
!       end do
!     end do

!     !-------------!
!     ! Block V2p1h !
!     !-------------!     
!
!     iab = 0
!     do i=nC+1,nO
!       do a=nO+1,nBas-nR
!         do b=nO+1,nBas-nR
!           iab = iab + 1   
!
!           H(1          ,1+n2h1p+iab) = sqrt(2d0)*ERI(p,i,b,a)
!           H(1+n2h1p+iab,1          ) = sqrt(2d0)*ERI(p,i,b,a)
!              
!         end do
!       end do
!     end do
!
!     !-------------!
!     ! Block C2h1p !
!     !-------------!
!
!     ija = 0
!     do i=nC+1,nO
!       do j=nC+1,nO
!         do a=nO+1,nBas-nR
!           ija = ija + 1
!              
!           klc = 0
!           do k=nC+1,nO
!             do l=nC+1,nO
!               do c=nO+1,nBas-nR
!                 klc = klc + 1
!                    
!                 H(1+ija,1+klc) & 
!                      = ((eHF(i) + eHF(j) - eHF(a))*Kronecker_delta(j,l)*Kronecker_delta(a,c) & 
!                      - 2d0*ERI(j,c,a,l))*Kronecker_delta(i,k)
!                    
!               end do
!             end do
!           end do
!
!         end do
!       end do
!     end do
!
!     !-------------!
!     ! Block C2p1h !
!     !-------------!
!     
!     iab = 0
!     do i=nC+1,nO
!       do a=nO+1,nBas-nR
!         do b=nO+1,nBas-nR
!           iab = iab + 1
!              
!           kcd = 0
!           do k=nC+1,nO
!             do c=nO+1,nBas-nR
!               do d=nO+1,nBas-nR
!                 kcd = kcd + 1
!                    
!                 H(1+n2h1p+iab,1+n2h1p+kcd) &
!                      = ((eHF(a) + eHF(b) - eHF(i))*Kronecker_delta(i,k)*Kronecker_delta(a,c) & 
!                      + 2d0*ERI(a,k,i,c))*Kronecker_delta(b,d)
!                    
!               end do
!             end do
!           end do
!       
!         end do
!       end do
!     end do
!      
      call wall_time(end_timing)
 
      timing = end_timing - start_timing
      write(*,*)
      write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for construction of supermatrix = ',timing,' seconds'
      write(*,*)

    else
 
      ! RPA for W
 
      write(*,*) 'Tamm-Dancoff approximation deactivated!'
      write(*,*)
 
      !---------------------------!
      !  Compute GW supermatrix   !
      !---------------------------!
      !                           !
      !     |   F   W2h1p W2p1h | ! 
      !     |                   | ! 
      ! H = | W2h1p D2h1p     0 | ! 
      !     |                   | ! 
      !     | W2p1h   0   D2p1h | ! 
      !                           !
      !---------------------------!
 
      ! Memory allocation

      allocate(Om1s(nVVs),X1s(nVVs,nVVs),Y1s(nOOs,nVVs),    &
               Om2s(nOOs),X2s(nVVs,nOOs),Y2s(nOOs,nOOs),    &
               rho1s(nBas,nBas,nVVs),rho2s(nBas,nBas,nOOs), &
               Om1t(nVVt),X1t(nVVt,nVVt),Y1t(nOOt,nVVt),    &
               Om2t(nOOt),X2t(nVVt,nOOt),Y2t(nOOt,nOOt),    &
               rho1t(nBas,nBas,nVVt),rho2t(nBas,nBas,nOOt))

      ! alpha-beta block
 
      ispin = 1
      iblock = 3

      ! Compute linear response

      allocate(Bpp(nVVs,nOOs),Cpp(nVVs,nVVs),Dpp(nOOs,nOOs))

      if(.not.TDA_T) call ppLR_B(iblock,nBas,nC,nO,nV,nR,nOOs,nVVs,1d0,ERI,Bpp)
                     call ppLR_C(iblock,nBas,nC,nO,nV,nR,nVVs,1d0,eHF,ERI,Cpp)
                     call ppLR_D(iblock,nBas,nC,nO,nV,nR,nOOs,1d0,eHF,ERI,Dpp)

      call ppLR(TDA_T,nOOs,nVVs,Bpp,Cpp,Dpp,Om1s,X1s,Y1s,Om2s,X2s,Y2s,EcRPA(ispin))

      deallocate(Bpp,Cpp,Dpp)

      call print_excitation_energies('ppRPA@RHF','2p (alpha-beta)',nVVs,Om1s(:))
      call print_excitation_energies('ppRPA@RHF','2h (alpha-beta)',nOOs,Om2s(:))

      ! alpha-alpha block

      ispin  = 2
      iblock = 4

      ! Compute linear response

      allocate(Bpp(nVVt,nOOt),Cpp(nVVt,nVVt),Dpp(nOOt,nOOt))
  
      if(.not.TDA_T) call ppLR_B(iblock,nBas,nC,nO,nV,nR,nOOt,nVVt,1d0,ERI,Bpp)
                     call ppLR_C(iblock,nBas,nC,nO,nV,nR,nVVt,1d0,eHF,ERI,Cpp)
                     call ppLR_D(iblock,nBas,nC,nO,nV,nR,nOOt,1d0,eHF,ERI,Dpp)
  
      call ppLR(TDA_T,nOOt,nVVt,Bpp,Cpp,Dpp,Om1t,X1t,Y1t,Om2t,X2t,Y2t,EcRPA(ispin))
  
      deallocate(Bpp,Cpp,Dpp)

      call print_excitation_energies('ppRPA@RHF','2p (alpha-alpha)',nVVt,Om1t)
      call print_excitation_energies('ppRPA@RHF','2h (alpha-beta)',nOOt,Om2t)

      !----------------------------------------------
      ! Compute excitation densities
      !----------------------------------------------

      iblock = 3
      call GTpp_excitation_density(iblock,nBas,nC,nO,nV,nR,nOOs,nVVs,ERI,X1s,Y1s,rho1s,X2s,Y2s,rho2s)

      iblock = 4
      call GTpp_excitation_density(iblock,nBas,nC,nO,nV,nR,nOOt,nVVt,ERI,X1t,Y1t,rho1t,X2t,Y2t,rho2t)

      call wall_time(start_timing)

      !---------!
      ! Block F !
      !---------!
 
      H(1,1) = eHF(p)
 
      !-------------!
      ! Block D2h1p !
      !-------------!
 
      ija = 0
      do i=nC+1,nO
        do ja=1,nS
          ija = ija + 1
 
          H(1+ija,1+ija) = eHF(i) - Om(ja) 
 
        end do
      end do
 
      !-------------!
      ! Block W2h1p !
      !-------------!
 
      ija = 0
      do i=nC+1,nO
        do ja=1,nS
          ija = ija + 1
 
          H(1    ,1+ija) = sqrt(2d0)*rho(p,i,ja)
          H(1+ija,1    ) = sqrt(2d0)*rho(p,i,ja)
 
        end do
      end do
 
      !-------------!
      ! Block D2p1h !
      !-------------!
 
      iab = 0
      do ia=1,nS
        do b=nO+1,nBas-nR
          iab = iab + 1
 
          H(1+n2h1p+iab,1+n2h1p+iab) = eHF(b) + Om(ia)
 
        end do
      end do
 
      !-------------!
      ! Block W2p1h !
      !-------------!
 
      iab = 0
      do ia=1,nS
        do b=nO+1,nBas-nR
          iab = iab + 1
 
          H(1          ,1+n2h1p+iab) = sqrt(2d0)*rho(p,b,ia)
          H(1+n2h1p+iab,1          ) = sqrt(2d0)*rho(p,b,ia)
 
        end do
      end do
       
      ! Memory deallocation

      deallocate(Om,Aph,Bph,XpY,XmY,rho)

      call wall_time(end_timing)
 
      timing = end_timing - start_timing
      write(*,*)
      write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for construction of supermatrix = ',timing,' seconds'
      write(*,*)
       
    end if
    
    !-------------------------!
    ! Diagonalize supermatrix !
    !-------------------------!
 
    call wall_time(start_timing)

    call diagonalize_matrix(nH,H,eGW)
 
    call wall_time(end_timing)

    timing = end_timing - start_timing
    write(*,*)
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for diagonalization of supermatrix = ',timing,' seconds'
    write(*,*)

    !-----------------!
    ! Compute weights !
    !-----------------!
 
      do s=1,nH
        Z(s) = H(1,s)**2
      end do
 
    !--------------!
    ! Dump results !
    !--------------!
 
    write(*,*)'-------------------------------------------'
    write(*,'(1X,A32,I3,A8)')'| G0W0 energies (eV) for orbital',p,'      |'
    write(*,*)'-------------------------------------------'
    write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X)') &
              '|','#','|','e_QP','|','Z','|'
    write(*,*)'-------------------------------------------'
  
    do s=1,nH
      if(eGW(s) < eF .and. eGW(s) > eF - window) then
!     if(Z(s) > cutoff1) then
        write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
        '|',s,'|',eGW(s)*HaToeV,'|',Z(s),'|'
      end if
    end do
  
    write(*,*)'-------------------------------------------'
    write(*,*)
 
    if(verbose) then 
 
      if(TDA_W) then  
 
        ! TDA printing format
 
        do s=1,nH
       
          if(eGW(s) < eF .and. eGW(s) > eF - window) then
       
            write(*,*)'-------------------------------------------------------------'
            write(*,'(1X,A7,1X,I3,A6,I3,A1,1X,A7,F12.6,A13,F6.4,1X)') & 
             'Orbital',p,' and #',s,':','e_QP = ',eGW(s)*HaToeV,' eV and Z = ',Z(s)
            write(*,*)'-------------------------------------------------------------'
            write(*,'(1X,A20,1X,A20,1X,A15,1X)') &
                      ' Configuration ',' Coefficient ',' Weight ' 
            write(*,*)'-------------------------------------------------------------'
           
            if(p <= nO) & 
              write(*,'(1X,A7,I3,A16,1X,F15.6,1X,F15.6)') &
              '      (',p,')               ',H(1,s),H(1,s)**2
            if(p > nO) & 
              write(*,'(1X,A16,I3,A7,1X,F15.6,1X,F15.6)') &
              '               (',p,')      ',H(1,s),H(1,s)**2
    
            ija = 0
            do i=nC+1,nO
              do j=nC+1,nO
                do a=nO+1,nBas-nR
                  ija = ija + 1
  
                  if(abs(H(1+ija,s)) > cutoff2)               &
                  write(*,'(1X,A3,I3,A1,I3,A6,I3,A7,1X,F15.6,1X,F15.6)') &
                  '  (',i,',',j,') -> (',a,')      ',H(1+ija,s),H(1+ija,s)**2
           
                end do
              end do
            end do
           
            iab = 0
            do i=nC+1,nO
              do a=nO+1,nBas-nR
                do b=nO+1,nBas-nR
                  iab = iab + 1
 
                  if(abs(H(1+n2h1p+iab,s)) > cutoff2)           &
                    write(*,'(1X,A7,I3,A6,I3,A1,I3,A3,1X,F15.6,1X,F15.6)') &
                    '      (',i,') -> (',a,',',b,')  ',H(1+n2h1p+iab,s),H(1+n2h1p+iab,s)**2
                  
                end do
              end do
            end do

            write(*,*)'-------------------------------------------------------------'
            write(*,*)

          end if

        end do
 
      else 
  
        ! non-TDA printing format
 
        do s=1,nH
        
          if(eGW(s) < eF .and. eGW(s) > eF - window) then
        
            write(*,*)'-------------------------------------------------------------'
            write(*,'(1X,A7,1X,I3,A6,I3,A1,1X,A7,F12.6,A13,F6.4,1X)') & 
             'Orbital',p,' and #',s,':','e_QP = ',eGW(s)*HaToeV,' eV and Z = ',Z(s)
            write(*,*)'-------------------------------------------------------------'
            write(*,'(1X,A20,1X,A20,1X,A15,1X)') &
                      ' Conf. (p,ia)  ',' Coefficient ',' Weight ' 
            write(*,*)'-------------------------------------------------------------'
           
            if(p <= nO) & 
              write(*,'(1X,A7,I3,A16,1X,F15.6,1X,F15.6)') &
              '      (',p,')               ',H(1,s),H(1,s)**2
            if(p > nO) & 
              write(*,'(1X,A16,I3,A7,1X,F15.6,1X,F15.6)') &
              '               (',p,')      ',H(1,s),H(1,s)**2
    
            ija = 0
            do i=nC+1,nO
              do ja=1,nS
                ija = ija + 1
  
                if(abs(H(1+ija,s)) > cutoff2)                     &
                write(*,'(1X,A7,I3,A1,I3,A12,1X,F15.6,1X,F15.6)') &
                '      (',i,',',ja,')           ',H(1+ija,s),H(1+ija,s)**2
           
              end do
            end do
           
            iab = 0
            do ia=1,nS
              do b=nO+1,nBas-nR
                iab = iab + 1
 
                  if(abs(H(1+n2h1p+iab,s)) > cutoff2)                 &
                    write(*,'(1X,A7,I3,A1,I3,A12,1X,F15.6,1X,F15.6)') &
                    '      (',ia,',',b,')           ',H(1+n2h1p+iab,s),H(1+n2h1p+iab,s)**2
                  
              end do
            end do
 
            write(*,*)'-------------------------------------------------------------'
            write(*,*)

          end if
 
        end do

      end if
 
    end if

  end do

end subroutine 
