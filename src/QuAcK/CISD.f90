subroutine CISD(singlet_manifold,triplet_manifold,nBas,nC,nO,nV,nR,ERI,eHF)

! Perform configuration interaction with singles and doubles

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: singlet_manifold
  logical,intent(in)            :: triplet_manifold
  integer,intent(in)            :: nBas,nC,nO,nV,nR
  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  logical                       :: dump_trans = .false.
  integer                       :: i,j,k,l
  integer                       :: a,b,c,d
  integer                       :: ia,jb,iajb,kcld
  integer                       :: ishift,jshift
  integer                       :: ispin
  integer                       :: nS
  integer                       :: nD
  integer                       :: nSD
  double precision,external     :: Kronecker_delta
  double precision,allocatable  :: H(:,:),Omega(:)

! Hello world

  write(*,*)
  write(*,*)'******************************************************'
  write(*,*)'| Configuration Interaction with Singles and Doubles |'
  write(*,*)'******************************************************'
  write(*,*)

! Compute CIS matrix

  if(singlet_manifold) then

    ispin = 1

    ! Dimensions

    nS  = (nO - nC)*(nV - nR)
    nD  = (nO - nC)*(nO - nC + 1)/2*(nV - nR)*(nV - nR + 1)/2
    nSD = 1 + nS + nD

    print*,'nS  = ',nS 
    print*,'nD  = ',nD 
    print*,'nSD = ',nSD 

    ! Memory allocation

    allocate(H(nSD,nSD),Omega(nSD))
    
    ! 0D block

    ishift = 0
    jshift = 1 + nS

    iajb = 0

    do i=nC+1,nO
      do a=1,nV-nR
        do j=i,nO
          do b=a,nV-nR

            iajb = iajb + 1
            H(ishift+1,jshift+iajb) = ERI(i,j,nO+a,nO+b)

          end do
        end do
      end do
    end do
    
    ! SS block

    ishift = 1
    jshift = 1

    ia = 0
    jb = 0

    do i=nC+1,nO
      do a=1,nV-nR

        ia = ia + 1

        do j=nC+1,nO
          do b=1,nV-nR

            jb = jb + 1

            H(ishift+ia,jshift+jb) & 
              = Kronecker_delta(i,j)*Kronecker_delta(a,b)*(eHF(nO+a) - eHF(i)) & 
              + ERI(nO+a,j,i,nO+b) -  ERI(nO+a,j,nO+b,i)

          end do
        end do

      end do
    end do
    
    ! SD block

    ishift = 1 
    jshift = 1 + nS

    ia = 0
    kcld = 0

    do i=nC+1,nO
      do a=1,nV-nR

        ia = ia + 1

        do k=nC+1,nO
          do c=1,nV-nR
            do l=k,nO
              do d=c,nV-nR

                kcld = kcld + 1

                H(ishift+ia,jshift+kcld)                                               & 
                = Kronecker_delta(i,k)*(ERI(nO+a,l,nO+c,nO+d) - ERI(nO+a,l,nO+d,nO+c)) &
                - Kronecker_delta(i,l)*(ERI(nO+a,k,nO+c,nO+d) - ERI(nO+a,k,nO+d,nO+c)) &
                - Kronecker_delta(a,c)*(ERI(k,l,i,nO+d)       - ERI(k,l,nO+d,i))       &
                + Kronecker_delta(a,d)*(ERI(k,l,i,nO+c)       - ERI(k,l,nO+c,i))

              end do
            end do
          end do
        end do

      end do
    end do

    ! DD block

    ishift = 1 + nS
    jshift = 1 + nS

    iajb = 0
    kcld = 0

    do i=nC+1,nO
      do a=1,nV-nR
        do j=i,nO
          do b=a,nV-nR

            iajb = iajb + 1

            do k=nC+1,nO
              do c=1,nV-nR
                do l=k,nO
                  do d=c,nV-nR
         
                    kcld = kcld + 1
         
!                   H(ishift+iajb,jshift+kcld)                            & 
!                   = Kronecker_delta(i,k)*(ERI(a,l,c,d) -  ERI(a,l,d,c)) &
!                   - Kronecker_delta(i,l)*(ERI(a,k,c,d) -  ERI(a,k,d,c)) &
!                   - Kronecker_delta(a,c)*(ERI(k,l,i,d) -  ERI(k,l,d,i)) &
!                   + Kronecker_delta(a,d)*(ERI(k,l,i,c) -  ERI(k,l,c,i))
         
                  end do
                end do
              end do
            end do

          end do
        end do
      end do
    end do

    call diagonalize_matrix(nSD,H,Omega)
    call print_excitation('CISD  ',ispin,nS,Omega)
 
    if(dump_trans) then
      print*,'Singlet CISD transition vectors'
      call matout(nSD,nSD,H)
      write(*,*)
    endif

  endif

! if(triplet_manifold) then

!   ispin = 2
!
!   call diagonalize_matrix(nSD,H,Omega)
!   call print_excitation('CISD  ',ispin,nSD,Omega)

!   if(dump_trans) then
!     print*,'Triplet CIS transition vectors'
!     call matout(nSD,nSD,H)
!     write(*,*)
!   endif

! endif

end subroutine CISD
