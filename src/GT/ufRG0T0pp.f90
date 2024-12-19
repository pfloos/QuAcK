subroutine ufRG0T0pp(dotest,TDA_T,nBas,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,eHF)

! Upfolded G0T0pp equations

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dotest

  logical,intent(in)            :: TDA_T
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
  integer                       :: ij,ab
  integer                       :: klc,kcd,ija,ijb,iab,jab

  logical                       :: print_T = .false.
  logical                       :: dRPA
  integer                       :: ispin
  integer                       :: nOOs,nOOt
  integer                       :: nVVs,nVVt
  double precision              :: EcRPA(nspin)
  integer                       :: n2h1p,n2p1h,nH
  double precision,external     :: Kronecker_delta
  double precision,allocatable  :: H(:,:)
  double precision,allocatable  :: eGT(:)
  double precision,allocatable  :: Z(:)
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

  logical                       :: verbose = .true.
  double precision,parameter    :: cutoff1 = 0.01d0
  double precision,parameter    :: cutoff2 = 0.01d0
  double precision              :: eF
  double precision,parameter    :: window = 1.5d0

  double precision              :: start_timing,end_timing,timing

! Output variables

! Hello world

  write(*,*)
  write(*,*)'******************************************'
  write(*,*)'* Restricted Upfolded G0T0pp Calculation *'
  write(*,*)'******************************************'
  write(*,*)

! Dimensions of the ppRPA linear reponse matrices

  nOOs = nO*(nO + 1)/2
  nVVs = nV*(nV + 1)/2

  nOOt = nO*(nO - 1)/2
  nVVt = nV*(nV - 1)/2

! Dimension of the supermatrix

  n2h1p = (nOOs+nOOt)*nV
  n2p1h = (nVVs+nVVt)*nO
  nH = 1 + n2h1p + n2p1h
  
! Initialization

  dRPA = .true.
  EcRPA = 0d0

  eF = 0.5d0*(eHF(nO+1) + eHF(nO))

!------------------!
! Compute T-matrix !
!------------------!
 
  if(.not. TDA_T) then

    ! Memory allocation

    allocate(Om1s(nVVs),X1s(nVVs,nVVs),Y1s(nOOs,nVVs),    &
             Om2s(nOOs),X2s(nVVs,nOOs),Y2s(nOOs,nOOs),    &
             rho1s(nBas,nBas,nVVs),rho2s(nBas,nBas,nOOs), &
             Om1t(nVVt),X1t(nVVt,nVVt),Y1t(nOOt,nVVt),    &
             Om2t(nOOt),X2t(nVVt,nOOt),Y2t(nOOt,nOOt),    &
             rho1t(nBas,nBas,nVVt),rho2t(nBas,nBas,nOOt))

    ! alpha-beta block
 
    ispin  = 1

    ! Compute linear response

    allocate(Bpp(nVVs,nOOs),Cpp(nVVs,nVVs),Dpp(nOOs,nOOs))

    call ppRLR_B(ispin,nBas,nC,nO,nV,nR,nOOs,nVVs,1d0,ERI,Bpp)
    call ppRLR_C(ispin,nBas,nC,nO,nV,nR,nVVs,1d0,eHF,ERI,Cpp)
    call ppRLR_D(ispin,nBas,nC,nO,nV,nR,nOOs,1d0,eHF,ERI,Dpp)

    call ppRLR(TDA_T,nOOs,nVVs,Bpp,Cpp,Dpp,Om1s,X1s,Y1s,Om2s,X2s,Y2s,EcRPA(ispin))

    if(print_T) call print_excitation_energies('ppRPA@RHF','2p (alpha-beta)',nVVs,Om1s(:))
    if(print_T) call print_excitation_energies('ppRPA@RHF','2h (alpha-beta)',nOOs,Om2s(:))

    ! Compute excitation densities

    call RGTpp_excitation_density(ispin,nBas,nC,nO,nV,nR,nOOs,nVVs,ERI,X1s,Y1s,rho1s,X2s,Y2s,rho2s)

    deallocate(Bpp,Cpp,Dpp,X1s,Y1s,X2s,Y2s)

    ! alpha-alpha block

    ispin  = 2

    ! Compute linear response

    allocate(Bpp(nVVt,nOOt),Cpp(nVVt,nVVt),Dpp(nOOt,nOOt))
  
    call ppRLR_B(ispin,nBas,nC,nO,nV,nR,nOOt,nVVt,1d0,ERI,Bpp)
    call ppRLR_C(ispin,nBas,nC,nO,nV,nR,nVVt,1d0,eHF,ERI,Cpp)
    call ppRLR_D(ispin,nBas,nC,nO,nV,nR,nOOt,1d0,eHF,ERI,Dpp)
  
    call ppRLR(TDA_T,nOOt,nVVt,Bpp,Cpp,Dpp,Om1t,X1t,Y1t,Om2t,X2t,Y2t,EcRPA(ispin))
  
    if(print_T) call print_excitation_energies('ppRPA@RHF','2p (alpha-alpha)',nVVt,Om1t)
    if(print_T) call print_excitation_energies('ppRPA@RHF','2h (alpha-beta)',nOOt,Om2t)

    ! Compute excitation densities

    call RGTpp_excitation_density(ispin,nBas,nC,nO,nV,nR,nOOt,nVVt,ERI,X1t,Y1t,rho1t,X2t,Y2t,rho2t)

    deallocate(Bpp,Cpp,Dpp,X1t,Y1t,X2t,Y2t)

  else

    allocate(rho1s(0,0,0),rho1t(0,0,0),rho2s(0,0,0),rho2t(0,0,0))

  end if

! Memory allocation

  allocate(H(nH,nH),eGT(nH),Z(nH))

!-------------------------!
! Main loop over orbitals !
!-------------------------!

  do p=nO,nO+1

    H(:,:) = 0d0

    if (TDA_T) then
 
      ! TDA for T
 
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
 
      !-------------!
      ! Block V2h1p !
      !-------------!

!     ija = 0
!     do i=nC+1,nO
!       do j=nC+1,nO
!         do a=nO+1,nBas-nR
!           ija = ija + 1
!              
!           H(1    ,1+ija) = ERI(p,a,i,j)
!           H(1+ija,1    ) = ERI(p,a,i,j)
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
 
      ! RPA for T
 
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

      call wall_time(start_timing)

      !---------!
      ! Block F !
      !---------!
 
      H(1,1) = eHF(p)
 
      !-------------!
      ! Block D2h1p !
      !-------------!
 
      ija = 0
      do ij=1,nOOs
        do a=nO+1,nBas-nR
          ija = ija + 1
 
          H(1+ija,1+ija) = - eHF(a) + Om2s(ij) 
 
        end do
      end do

      do ij=1,nOOt
        do a=nO+1,nBas-nR
          ija = ija + 1
 
          H(1+ija,1+ija) = - eHF(a) + Om2t(ij) 
 
        end do
      end do
 
      !-------------!
      ! Block W2h1p !
      !-------------!
 
      ija = 0
      do ij=1,nOOs
        do a=nO+1,nBas-nR
          ija = ija + 1
 
          H(1    ,1+ija) = rho2s(p,a,ij)
          H(1+ija,1    ) = rho2s(p,a,ij)
 
        end do
      end do

      do ij=1,nOOt
        do a=nO+1,nBas-nR
          ija = ija + 1
 
          H(1    ,1+ija) = rho2t(p,a,ij)
          H(1+ija,1    ) = rho2t(p,a,ij)
 
        end do
      end do
 
      !-------------!
      ! Block D2p1h !
      !-------------!
 
      iab = 0
      do ab=1,nVVs
        do i=nC+1,nO
          iab = iab + 1
 
          H(1+n2h1p+iab,1+n2h1p+iab) = - eHF(i) + Om1s(ab)
 
        end do
      end do

      do ab=1,nVVt
        do i=nC+1,nO
          iab = iab + 1
 
          H(1+n2h1p+iab,1+n2h1p+iab) = - eHF(i) + Om1t(ab)
 
        end do
      end do
 
      !-------------!
      ! Block W2p1h !
      !-------------!
 
      iab = 0
      do ab=1,nVVs
        do i=nC+1,nO
          iab = iab + 1
 
          H(1          ,1+n2h1p+iab) = rho1s(p,i,ab)
          H(1+n2h1p+iab,1          ) = rho1s(p,i,ab)
 
        end do
      end do

      do ab=1,nVVt
        do i=nC+1,nO
          iab = iab + 1
 
          H(1          ,1+n2h1p+iab) = rho1t(p,i,ab)
          H(1+n2h1p+iab,1          ) = rho1t(p,i,ab)
 
        end do
      end do

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

    call diagonalize_matrix(nH,H,eGT)
 
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
    write(*,'(1X,A34,I3,A6)')'| G0T0pp energies (eV) for orbital',p,'    |'
    write(*,*)'-------------------------------------------'
    write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X)') &
              '|','#','|','e_QP','|','Z','|'
    write(*,*)'-------------------------------------------'
  
    do s=1,nH
     if(eGT(s) < eF .and. eGT(s) > eF - window) then
      !if(Z(s) > cutoff1) then
        write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
        '|',s,'|',eGT(s)*HaToeV,'|',Z(s),'|'
      end if
    end do
  
    write(*,*)'-------------------------------------------'
    write(*,*)
 
    if(verbose) then 
 
      if(TDA_T) then  
 
        ! TDA printing format
 
!       do s=1,nH
!      
!         if(eGT(s) < eF .and. eGT(s) > eF - window) then
!      
!           write(*,*)'-------------------------------------------------------------'
!           write(*,'(1X,A7,1X,I3,A6,I3,A1,1X,A7,F12.6,A13,F6.4,1X)') & 
!            'Orbital',p,' and #',s,':','e_QP = ',eGT(s)*HaToeV,' eV and Z = ',Z(s)
!           write(*,*)'-------------------------------------------------------------'
!           write(*,'(1X,A20,1X,A20,1X,A15,1X)') &
!                     ' Configuration ',' Coefficient ',' Weight ' 
!           write(*,*)'-------------------------------------------------------------'
!          
!           if(p <= nO) & 
!             write(*,'(1X,A7,I3,A16,1X,F15.6,1X,F15.6)') &
!             '      (',p,')               ',H(1,s),H(1,s)**2
!           if(p > nO) & 
!             write(*,'(1X,A16,I3,A7,1X,F15.6,1X,F15.6)') &
!             '               (',p,')      ',H(1,s),H(1,s)**2
!   
!           ija = 0
!           do i=nC+1,nO
!             do j=nC+1,nO
!               do a=nO+1,nBas-nR
!                 ija = ija + 1
! 
!                 if(abs(H(1+ija,s)) > cutoff2)               &
!                 write(*,'(1X,A3,I3,A1,I3,A6,I3,A7,1X,F15.6,1X,F15.6)') &
!                 '  (',i,',',j,') -> (',a,')      ',H(1+ija,s),H(1+ija,s)**2
!          
!               end do
!             end do
!           end do
!          
!           iab = 0
!           do i=nC+1,nO
!             do a=nO+1,nBas-nR
!               do b=nO+1,nBas-nR
!                 iab = iab + 1
!
!                 if(abs(H(1+n2h1p+iab,s)) > cutoff2)           &
!                   write(*,'(1X,A7,I3,A6,I3,A1,I3,A3,1X,F15.6,1X,F15.6)') &
!                   '      (',i,') -> (',a,',',b,')  ',H(1+n2h1p+iab,s),H(1+n2h1p+iab,s)**2
!                 
!               end do
!             end do
!           end do

!           write(*,*)'-------------------------------------------------------------'
!           write(*,*)

!         end if

!       end do
 
      else 
  
        ! non-TDA printing format
 
        do s=1,nH
        
          if(eGT(s) < eF .and. eGT(s) > eF - window) then
          !if(Z(s) > cutoff2) then
        
            write(*,*)'-------------------------------------------------------------'
            write(*,'(1X,A7,1X,I3,A6,I3,A1,1X,A7,F12.6,A13,F6.4,1X)') & 
             'Orbital',p,' and #',s,':','e_QP = ',eGT(s)*HaToeV,' eV and Z = ',Z(s)
            write(*,*)'-------------------------------------------------------------'
            write(*,'(1X,A24,1X,A20,1X,A15,1X)') &
                      ' Conf. (i,ab) or (a,ij) ',' Coefficient ',' Weight ' 
            write(*,*)'-------------------------------------------------------------'
           
            if(p <= nO) & 
              write(*,'(1X,A7,I3,A16,1X,F15.6,1X,F15.6,1X,F12.6)') &
              '      (',p,')               ',H(1,s),H(1,s)**2,-eHF(p)*HaToeV
            if(p > nO) & 
              write(*,'(1X,A16,I3,A7,1X,F15.6,1X,F15.6,1X,F12.6)') &
              '               (',p,')      ',H(1,s),H(1,s)**2,-eHF(p)*HaToeV
    
            ija = 0
            do ij=1,nOOs+nOOt
              do a=nO+1,nBas-nR
                ija = ija + 1
  
                if(abs(H(1+ija,s)) > cutoff2 .and. ij<nOOs+1)                     &
                write(*,'(1X,A7,I3,A1,I3,A12,1X,F15.6,1X,F15.6,1X,F12.6)') &
                '      (',a,',',ij,')           ',H(1+ija,s),H(1+ija,s)**2,(- eHF(a) + Om2s(ij))*HaToeV
                if(abs(H(1+ija,s)) > cutoff2 .and. ij>nOOs)                     &
                write(*,'(1X,A7,I3,A1,I3,A12,1X,F15.6,1X,F15.6,1X,F12.6)') &
                '      (',a,',',ij,')           ',H(1+ija,s),H(1+ija,s)**2,(- eHF(a) + Om2t(ij-nOOs))*HaToeV
            
              end do
            end do
           
            iab = 0
            do ab=1,nVVs+nVVt
              do i=nC+1,nO
                iab = iab + 1
 
                  if(abs(H(1+n2h1p+iab,s)) > cutoff2 .and. ab<nVVs+1)                 &
                    write(*,'(1X,A7,I3,A1,I3,A12,1X,F15.6,1X,F15.6,1X,F12.6)') &
                    '      (',i,',',ab,')           ',H(1+n2h1p+iab,s),H(1+n2h1p+iab,s)**2,(-eHF(i) + Om1s(ab))*HaToeV
                                    if(abs(H(1+n2h1p+iab,s)) > cutoff2 .and. ab>nVVs)                 &
                    write(*,'(1X,A7,I3,A1,I3,A12,1X,F15.6,1X,F15.6,1X,F12.6)') &
                    '      (',i,',',ab,')           ',H(1+n2h1p+iab,s),H(1+n2h1p+iab,s)**2,(-eHF(i) + Om1t(ab-nVVs))*HaToeV
                  
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
