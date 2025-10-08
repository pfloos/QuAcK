subroutine RGW_ppBSE_qs(TDA_W,TDA,dBSE,dTDA,singlet,triplet,eta,nOrb,nC,nO,nV,nR,nS,ERI,dipole_int,eW,eGW,EcBSE)

! Compute the Bethe-Salpeter excitation energies at the pp level

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: TDA_W
  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA
  logical,intent(in)            :: singlet
  logical,intent(in)            :: triplet

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: eW(nOrb)
  double precision,intent(in)   :: eGW(nOrb)
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: dipole_int(nOrb,nOrb,ncart)

! Local variables

  integer                       :: ispin
  integer                       :: isp_W

  logical                       :: dRPA   = .false.
  logical                       :: dRPA_W = .true.

  integer                       :: nOO
  integer                       :: nVV

  integer                       :: p, q, m
  integer                       :: i_data, supp_data_dbl_size, supp_data_int_size
  integer                       :: n_states, n_states_diag

  double precision              :: tt1, tt2

  double precision,allocatable  :: Aph(:,:)
  double precision,allocatable  :: Bph(:,:)

  double precision              :: EcRPA
  double precision,allocatable  :: OmRPA(:)
  double precision,allocatable  :: XpY_RPA(:,:)
  double precision,allocatable  :: XmY_RPA(:,:)
  double precision,allocatable  :: rho_RPA(:,:,:)

  double precision,allocatable  :: Bpp(:,:)
  double precision,allocatable  :: Cpp(:,:)
  double precision,allocatable  :: Dpp(:,:)

  double precision,allocatable  :: Om1(:)
  double precision,allocatable  :: X1(:,:)
  double precision,allocatable  :: Y1(:,:)

  double precision,allocatable  :: Om2(:)
  double precision,allocatable  :: X2(:,:)
  double precision,allocatable  :: Y2(:,:)

  double precision,allocatable  :: KB_sta(:,:)
  double precision,allocatable  :: KC_sta(:,:)
  double precision,allocatable  :: KD_sta(:,:)
  
  integer,         allocatable  :: supp_data_int(:)
  double precision,allocatable  :: supp_data_dbl(:)
  double precision,allocatable  :: Om(:), R(:,:)

! Output variables

  double precision,intent(out)  :: EcBSE(nspin)

!-----!
! TDA !
!-----!

  if(TDA) then
    write(*,*) 'Tamm-Dancoff approximation activated in ppBSE!'
    write(*,*)
  end if

! Initialization

  EcBSE(:) = 0d0

!---------------------------------
! Compute (singlet) RPA screening 
!---------------------------------

  isp_W = 1

  allocate(OmRPA(nS),XpY_RPA(nS,nS),XmY_RPA(nS,nS),rho_RPA(nOrb,nOrb,nS), &
           Aph(nS,nS),Bph(nS,nS))
 
                 call phRLR_A(isp_W,dRPA_W,nOrb,nC,nO,nV,nR,nS,1d0,eW,ERI,Aph)
  if(.not.TDA_W) call phRLR_B(isp_W,dRPA_W,nOrb,nC,nO,nV,nR,nS,1d0,ERI,Bph)

  call phRLR(TDA_W,nS,Aph,Bph,EcRPA,OmRPA,XpY_RPA,XmY_RPA)
!  call phLR_transition_vectors(.true.,nOrb,nC,nO,nV,nR,nS,dipole_int,OmRPA,XpY_RPA,XmY_RPA)

  call RGW_excitation_density(nOrb,nC,nO,nR,nS,ERI,XpY_RPA,rho_RPA)

  deallocate(XpY_RPA,XmY_RPA,Aph,Bph)

!-------------------
! Singlet manifold
!-------------------

 if(singlet) then

    write(*,*) '****************'
    write(*,*) '*** Singlets ***'
    write(*,*) '****************'
    write(*,*) 

    ispin = 1

    nOO = nO*(nO+1)/2
    nVV = nV*(nV+1)/2

    allocate(Om1(nVV),X1(nVV,nVV),Y1(nOO,nVV))
    allocate(Om2(nOO),X2(nVV,nOO),Y2(nOO,nOO))

    ! Compute BSE excitation energies

    ! ---
    ! LAPACK
    ! ---

    allocate(Bpp(nVV,nOO),Cpp(nVV,nVV),Dpp(nOO,nOO))
    allocate(KB_sta(nVV,nOO),KC_sta(nVV,nVV),KD_sta(nOO,nOO))

    KB_sta(:,:) = 0d0
    KC_sta(:,:) = 0d0
  
    call RGW_ppBSE_static_kernel_C_qs(ispin,eta,nOrb,nC,nO,nV,nR,nS,nVV,1d0,ERI,eGW,OmRPA,rho_RPA,KC_sta)
    call RGW_ppBSE_static_kernel_D_qs(ispin,eta,nOrb,nC,nO,nV,nR,nS,nOO,1d0,ERI,eGW,OmRPA,rho_RPA,KD_sta)
    if(.not.TDA) then
       call RGW_ppBSE_dynamic_kernel_B(ispin,eta,nOrb,nC,nO,nV,nR,nS,nOO,nVV,1d0,eGW,OmRPA,rho_RPA,KB_sta)
    endif
    
                 call ppRLR_C(ispin,nOrb,nC,nO,nV,nR,nVV,1d0,eGW,ERI,Cpp)
                 call ppRLR_D(ispin,nOrb,nC,nO,nV,nR,nOO,1d0,eGW,ERI,Dpp)
    if(.not.TDA) call ppRLR_B(ispin,nOrb,nC,nO,nV,nR,nOO,nVV,1d0,ERI,Bpp)

    Bpp(:,:) = Bpp(:,:) + KB_sta(:,:)
    Cpp(:,:) = Cpp(:,:) + KC_sta(:,:)
    Dpp(:,:) = Dpp(:,:) + KD_sta(:,:)

    call ppRLR(TDA,nOO,nVV,Bpp,Cpp,Dpp,Om1,X1,Y1,Om2,X2,Y2,EcBSE(ispin))

    deallocate(Bpp,Cpp,Dpp,KB_sta,KC_sta,KD_sta)
!
!    print*, 'LAPACK:'
!    print*, Om2
!    print*, Om1

    ! ---



    !! ---
    !! Davidson
    !! ---

    !n_states = nOO + 5
    !n_states_diag = n_states + 4
    !allocate(Om(nOO+nVV), R(nOO+nVV,n_states_diag))

    !supp_data_int_size = 1
    !allocate(supp_data_int(supp_data_int_size))
    !supp_data_int(1) = nS

    !supp_data_dbl_size = nS + nOrb*nOrb*nS + 1
    !allocate(supp_data_dbl(supp_data_dbl_size))
    !! scalars
    !supp_data_dbl(1) = eta
    !i_data = 1
    !! rho_RPA
    !do q = 1, nOrb
    !  do p = 1, nOrb
    !    do m = 1, nS
    !      i_data = i_data + 1
    !      supp_data_dbl(i_data) = rho_RPA(p,q,m)
    !    enddo
    !  enddo
    !enddo
    !! OmRPA
    !do m = 1, nS
    !  i_data = i_data + 1
    !  supp_data_dbl(i_data) = OmRPA(m)
    !enddo

    !call ppLR_davidson(ispin, TDA, nC, nO, nR, nOrb, nOO, nVV, &
    !                   1.d0,                                   & ! lambda
    !                   eGW(1),                                 &
    !                   0.d0,                                   & ! eF
    !                   ERI(1,1,1,1),                           &
    !                   supp_data_int(1), supp_data_int_size,   &
    !                   supp_data_dbl(1), supp_data_dbl_size,   &
    !                   Om(1), R(1,1), n_states, n_states_diag, "GW", 1)

    !deallocate(Om, R)
    !deallocate(supp_data_dbl, supp_data_int)
    !stop

    ! ---

    call ppLR_transition_vectors(.true.,nOrb,nC,nO,nV,nR,nOO,nVV,dipole_int,Om1,X1,Y1,Om2,X2,Y2)

    deallocate(Om1,X1,Y1,Om2,X2,Y2)
    
  end if

!-------------------
! Triplet manifold
!-------------------

 if(triplet) then

    write(*,*) '****************'
    write(*,*) '*** Triplets ***'
    write(*,*) '****************'
    write(*,*) 

    ispin = 2

    nOO = nO*(nO-1)/2
    nVV = nV*(nV-1)/2

    allocate(Om1(nVV),X1(nVV,nVV),Y1(nOO,nVV))
    allocate(Om2(nOO),X2(nVV,nOO),Y2(nOO,nOO))

    ! Compute BSE excitation energies

    ! ---
    ! LAPACK
    ! ---

    allocate(Bpp(nVV,nOO),Cpp(nVV,nVV),Dpp(nOO,nOO))
    allocate(KB_sta(nVV,nOO),KC_sta(nVV,nVV),KD_sta(nOO,nOO))
    
    KB_sta(:,:) = 0d0
    KC_sta(:,:) = 0d0

    call RGW_ppBSE_static_kernel_C_qs(ispin,eta,nOrb,nC,nO,nV,nR,nS,nVV,1d0,ERI,eGW,OmRPA,rho_RPA,KC_sta)
    call RGW_ppBSE_static_kernel_D_qs(ispin,eta,nOrb,nC,nO,nV,nR,nS,nOO,1d0,ERI,eGW,OmRPA,rho_RPA,KD_sta)
    if(.not.TDA) call RGW_ppBSE_dynamic_kernel_B(ispin,eta,nOrb,nC,nO,nV,nR,nS,nOO,nVV,1d0,eGW,OmRPA,rho_RPA,KB_sta)

                 call ppRLR_C(ispin,nOrb,nC,nO,nV,nR,nVV,1d0,eGW,ERI,Cpp)
                 call ppRLR_D(ispin,nOrb,nC,nO,nV,nR,nOO,1d0,eGW,ERI,Dpp)
    if(.not.TDA) call ppRLR_B(ispin,nOrb,nC,nO,nV,nR,nOO,nVV,1d0,ERI,Bpp)

    Bpp(:,:) = Bpp(:,:) + KB_sta(:,:)
    Cpp(:,:) = Cpp(:,:) + KC_sta(:,:)
    Dpp(:,:) = Dpp(:,:) + KD_sta(:,:)

    call ppRLR(TDA,nOO,nVV,Bpp,Cpp,Dpp,Om1,X1,Y1,Om2,X2,Y2,EcBSE(ispin))
    deallocate(Bpp,Cpp,Dpp)

    ! ---
    ! Davidson
    ! ---

    !n_states = nOO + 5
    !n_states_diag = n_states + 4
    !allocate(Om(nOO+nVV), R(nOO+nVV,n_states_diag))

    !supp_data_int_size = 1
    !allocate(supp_data_int(supp_data_int_size))
    !supp_data_int(1) = nS

    !supp_data_dbl_size = nS + nOrb*nOrb*nS + 1
    !allocate(supp_data_dbl(supp_data_dbl_size))
    !! scalars
    !supp_data_dbl(1) = eta
    !i_data = 1
    !! rho_RPA
    !do q = 1, nOrb
    !  do p = 1, nOrb
    !    do m = 1, nS
    !      i_data = i_data + 1
    !      supp_data_dbl(i_data) = rho_RPA(p,q,m)
    !    enddo
    !  enddo
    !enddo
    !! OmRPA
    !do m = 1, nS
    !  i_data = i_data + 1
    !  supp_data_dbl(i_data) = OmRPA(m)
    !enddo

    !call ppLR_davidson(ispin, TDA, nC, nO, nR, nOrb, nOO, nVV, &
    !                   1.d0,                                   & ! lambda
    !                   eGW(1),                                 &
    !                   0.d0,                                   & ! eF
    !                   ERI(1,1,1,1),                           &
    !                   supp_data_int(1), supp_data_int_size,   &
    !                   supp_data_dbl(1), supp_data_dbl_size,   &
    !                   Om(1), R(1,1), n_states, n_states_diag, "GW", 1)

    !deallocate(Om, R)
    !deallocate(supp_data_dbl, supp_data_int)
    !stop

    ! ---

    EcBSE(ispin) = 3d0*EcBSE(ispin)

    call ppLR_transition_vectors(.false.,nOrb,nC,nO,nV,nR,nOO,nVV,dipole_int,Om1,X1,Y1,Om2,X2,Y2)

    deallocate(KB_sta,KC_sta,KD_sta)
    deallocate(Om1,X1,Y1,Om2,X2,Y2)

  end if

end subroutine RGW_ppBSE_qs

subroutine RGW_ppBSE_static_kernel_C_qs(ispin,eta,nOrb,nC,nO,nV,nR,nS,nVV,lambda,ERI,eGW,Om,rho,KC)

! Compute the VVVV block of the static screening W for the pp-BSE

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: ispin
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  integer,intent(in)            :: nVV
  double precision,intent(in)   :: eta
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: eGW(nOrb)
  double precision,intent(in)   :: Om(nS)
  double precision,intent(in)   :: rho(nOrb,nOrb,nS)

! Local variables

  double precision,external     :: Kronecker_delta
  double precision              :: eps,num
  double precision              :: chi,dem
  double precision              :: tmp_ab, lambda4, eta2
  integer                       :: a,b,c,d,ab,cd,m
  integer                       :: a0, aa

  double precision, allocatable :: Om_tmp(:)
  double precision, allocatable :: tmp_m(:,:,:)
  double precision, allocatable :: tmp(:,:,:,:)

! Output variables

  double precision,intent(out)  :: KC(nVV,nVV)

!---------------!
! Singlet block !
!---------------!

  if(ispin == 1) then

!    --- --- ---
!    OpenMP implementation
!    --- --- ---

   a0 = nOrb - nR - nO
   lambda4 = 4.d0 * lambda
   eta2 = eta * eta

   !$OMP PARALLEL DEFAULT(NONE)                              &
   !$OMP          PRIVATE(a, b, aa, ab, c, d, cd, m, tmp_ab) &
   !$OMP          SHARED(eta, nO, nOrb, nR, nS, a0, lambda4, Om, rho, KC, num, dem, eGW)
   !$OMP DO
   do a = nO+1, nOrb-nR
     aa = a0 * (a - nO - 1) - (a - nO - 1) * (a - nO) / 2 - nO
     do b = a, nOrb-nR
       ab = aa + b

       cd = 0
       do c = nO+1, nOrb-nR
         do d = c, nOrb-nR
           cd = cd + 1

           KC(ab,cd) = 0d0
           
           do m = 1, nS

              num = (rho(a,c,m)*rho(b,d,m) + rho(b,c,m)*rho(a,d,m))/2d0
             
              dem = - Om(m) + eGW(d) - eGW(a)
              KC(ab,cd) = KC(ab,cd) + num*dem/(dem**2 + eta**2)
             
              dem = - Om(m) + eGW(a) - eGW(d)
              KC(ab,cd) = KC(ab,cd) + num*dem/(dem**2 + eta**2)
              
              dem = - Om(m) + eGW(d) - eGW(b)
              KC(ab,cd) = KC(ab,cd) + num*dem/(dem**2 + eta**2)
               
              dem = - Om(m) + eGW(b) + eGW(d)
              KC(ab,cd) = KC(ab,cd) + num*dem/(dem**2 + eta**2)
              
              dem = - Om(m) + eGW(c) - eGW(a)
              KC(ab,cd) = KC(ab,cd) + num*dem/(dem**2 + eta**2)
              
              dem = - Om(m) + eGW(a) - eGW(c)
              KC(ab,cd) = KC(ab,cd) + num*dem/(dem**2 + eta**2)
              
              dem = - Om(m) + eGW(c) - eGW(b)
              KC(ab,cd) = KC(ab,cd) + num*dem/(dem**2 + eta**2)
              
              dem = - Om(m) + eGW(b) - eGW(c)
              KC(ab,cd) = KC(ab,cd) + num*dem/(dem**2 + eta**2)
              
           end do

           KC(ab,cd) = KC(ab,cd)/sqrt((1d0 + Kronecker_delta(a,b))*(1d0 + Kronecker_delta(c,d)))
          
         enddo
       enddo
     enddo
   enddo
   !$OMP END DO
   !$OMP END PARALLEL
!    --- --- ---


!    --- --- ---
!    Naive implementation
!    --- --- ---
!
!    ab = 0
!    do a=nO+1,nOrb-nR
!      do b=a,nOrb-nR
!        ab = ab + 1
!        cd = 0
!        do c=nO+1,nOrb-nR
!          do d=c,nOrb-nR
!            cd = cd + 1
!
!              chi = 0d0
!              do m=1,nS
!                eps = Om(m)**2 + eta**2
!                chi = chi - rho(a,c,m)*rho(b,d,m)*Om(m)/eps &
!                          - rho(a,d,m)*rho(b,c,m)*Om(m)/eps
!              end do
!
!              KC(ab,cd) = 4d0*lambda*chi/sqrt((1d0 + Kronecker_delta(a,b))*(1d0 + Kronecker_delta(c,d)))
!
!          end do
!        end do
!      end do
!    end do
!    --- --- ---

  end if

!---------------!
! Triplet block !
!---------------!

  if(ispin == 2) then

!    --- --- ---
!    OpenMP implementation
!    --- --- ---

   a0 = nOrb - nR - nO - 1
   lambda4 = 4.d0 * lambda
   eta2 = eta * eta

   !$OMP PARALLEL DEFAULT(NONE)                      &
   !$OMP          PRIVATE(a, b, aa, ab, c, d, cd, m) &
   !$OMP          SHARED(eta, nO, nOrb, nR, nS, a0, lambda4, Om, rho, KC, num, dem, eGW)
   !$OMP DO
   do a = nO+1, nOrb-nR
     aa = a0 * (a - nO - 1) - (a - nO - 1) * (a - nO) / 2 - nO - 1
     do b = a+1, nOrb-nR
       ab = aa + b

       cd = 0
       do c = nO+1, nOrb-nR
         do d = c+1, nOrb-nR
           cd = cd + 1

           KC(ab,cd) = 0d0
           
           do m = 1, nS

              num = (rho(a,c,m)*rho(b,d,m) - rho(b,c,m)*rho(a,d,m))/2d0
             
              dem = - Om(m) + eGW(d) - eGW(a)
              KC(ab,cd) = KC(ab,cd) + num*dem/(dem**2 + eta**2)
             
              dem = - Om(m) + eGW(a) - eGW(d)
              KC(ab,cd) = KC(ab,cd) + num*dem/(dem**2 + eta**2)
              
              dem = - Om(m) + eGW(d) - eGW(b)
              KC(ab,cd) = KC(ab,cd) + num*dem/(dem**2 + eta**2)
               
              dem = - Om(m) + eGW(b) + eGW(d)
              KC(ab,cd) = KC(ab,cd) + num*dem/(dem**2 + eta**2)
              
              dem = - Om(m) + eGW(c) - eGW(a)
              KC(ab,cd) = KC(ab,cd) + num*dem/(dem**2 + eta**2)
              
              dem = - Om(m) + eGW(a) - eGW(c)
              KC(ab,cd) = KC(ab,cd) + num*dem/(dem**2 + eta**2)
              
              dem = - Om(m) + eGW(c) - eGW(b)
              KC(ab,cd) = KC(ab,cd) + num*dem/(dem**2 + eta**2)
              
              dem = - Om(m) + eGW(b) - eGW(c)
              KC(ab,cd) = KC(ab,cd) + num*dem/(dem**2 + eta**2)
              
           end do
           
         enddo
       enddo
     enddo
   enddo
   !$OMP END DO
   !$OMP END PARALLEL
   
!    --- --- ---


!    --- --- ---
!    Naive implementation
!    --- --- ---
!
!    ab = 0
!    do a=nO+1,nOrb-nR
!      do b=a+1,nOrb-nR
!        ab = ab + 1
!        cd = 0
!        do c=nO+1,nOrb-nR
!          do d=c+1,nOrb-nR
!            cd = cd + 1
!
!            chi = 0d0
!            do m=1,nS
!              eps = Om(m)**2 + eta**2
!              chi = chi - rho(a,c,m)*rho(b,d,m)*Om(m)/eps &
!                        + rho(a,d,m)*rho(b,c,m)*Om(m)/eps
!            end do
!           
!            KC(ab,cd) = 4d0*lambda*chi
!
!          end do
!        end do
!      end do
!    end do
!    --- --- ---

  end if

end subroutine 


subroutine RGW_ppBSE_static_kernel_D_qs(ispin,eta,nOrb,nC,nO,nV,nR,nS,nOO,lambda,ERI,eGW,Om,rho,KD)

! Compute the OOOO block of the static screening W for the pp-BSE

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: ispin
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  integer,intent(in)            :: nOO
  double precision,intent(in)   :: eta
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: eGW(nOrb)
  double precision,intent(in)   :: Om(nS)
  double precision,intent(in)   :: rho(nOrb,nOrb,nS)

! Local variables

  double precision,external     :: Kronecker_delta
  double precision              :: dem
  double precision              :: num
  integer                       :: i,j,k,l,ij,kl,m

! Output variables

  double precision,intent(out)  :: KD(nOO,nOO)

! Initialization

  KD(:,:) = 0d0

!---------------!
! Singlet block !
!---------------!

  if(ispin == 1) then

    ij = 0
    do i=nC+1,nO
      do j=i,nO
        ij = ij + 1
        kl = 0
        do k=nC+1,nO
          do l=k,nO
            kl = kl + 1

            do m=1,nS

               num = (rho(i,k,m)*rho(j,l,m) + rho(j,k,m)*rho(i,l,m))/2d0
             
               dem = - Om(m) + eGW(l) - eGW(i)
               KD(ij,kl) = KD(ij,kl) + num*dem/(dem**2 + eta**2)
             
               dem = - Om(m) + eGW(i) - eGW(l)
               KD(ij,kl) = KD(ij,kl) + num*dem/(dem**2 + eta**2)
               
               dem = - Om(m) + eGW(l) - eGW(j)
               KD(ij,kl) = KD(ij,kl) + num*dem/(dem**2 + eta**2)
               
               dem = - Om(m) + eGW(j) + eGW(l)
               KD(ij,kl) = KD(ij,kl) + num*dem/(dem**2 + eta**2)
               
               dem = - Om(m) + eGW(k) - eGW(i)
               KD(ij,kl) = KD(ij,kl) + num*dem/(dem**2 + eta**2)
               
               dem = - Om(m) + eGW(i) - eGW(k)
               KD(ij,kl) = KD(ij,kl) + num*dem/(dem**2 + eta**2)
               
               dem = - Om(m) + eGW(k) - eGW(j)
               KD(ij,kl) = KD(ij,kl) + num*dem/(dem**2 + eta**2)
               
               dem = - Om(m) + eGW(j) - eGW(k)
               KD(ij,kl) = KD(ij,kl) + num*dem/(dem**2 + eta**2)

            end do
 
            KD(ij,kl) = KD(ij,kl)/sqrt((1d0 + Kronecker_delta(i,j))*(1d0 + Kronecker_delta(k,l)))

          end do
        end do
      end do
    end do

  end if

!---------------!
! Triplet block !
!---------------!

  if(ispin == 2) then

    ij = 0
    do i=nC+1,nO
      do j=i+1,nO
        ij = ij + 1
        kl = 0
        do k=nC+1,nO
          do l=k+1,nO
            kl = kl + 1

            do m=1,nS

               num = (rho(i,k,m)*rho(j,l,m) - rho(j,k,m)*rho(i,l,m))/2d0
             
               dem = - Om(m) + eGW(l) - eGW(i)
               KD(ij,kl) = KD(ij,kl) + num*dem/(dem**2 + eta**2)
             
               dem = - Om(m) + eGW(i) - eGW(l)
               KD(ij,kl) = KD(ij,kl) + num*dem/(dem**2 + eta**2)
               
               dem = - Om(m) + eGW(l) - eGW(j)
               KD(ij,kl) = KD(ij,kl) + num*dem/(dem**2 + eta**2)
               
               dem = - Om(m) + eGW(j) + eGW(l)
               KD(ij,kl) = KD(ij,kl) + num*dem/(dem**2 + eta**2)
               
               dem = - Om(m) + eGW(k) - eGW(i)
               KD(ij,kl) = KD(ij,kl) + num*dem/(dem**2 + eta**2)
               
               dem = - Om(m) + eGW(i) - eGW(k)
               KD(ij,kl) = KD(ij,kl) + num*dem/(dem**2 + eta**2)
               
               dem = - Om(m) + eGW(k) - eGW(j)
               KD(ij,kl) = KD(ij,kl) + num*dem/(dem**2 + eta**2)
               
               dem = - Om(m) + eGW(j) - eGW(k)
               KD(ij,kl) = KD(ij,kl) + num*dem/(dem**2 + eta**2)

            end do

          end do
        end do
      end do
    end do

  end if

end subroutine RGW_ppBSE_static_kernel_D_qs
