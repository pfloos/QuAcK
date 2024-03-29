subroutine UMP2(dotest,nBas,nC,nO,nV,nR,ERI_aa,ERI_ab,ERI_bb,ENuc,EUHF,eHF,Ec)

! Perform unrestricted second-order Moller-Plesset calculation

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dotest

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC(nspin)
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nV(nspin)
  integer,intent(in)            :: nR(nspin)
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: EUHF
  double precision,intent(in)   :: ERI_aa(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_ab(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_bb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: eHF(nBas,nspin)

! Local variables

  integer                       :: bra,ket
  integer                       :: i,j,a,b
  double precision              :: eps
  double precision              :: Edaa,Exaa,Ecaa
  double precision              :: Edab,Exab,Ecab
  double precision              :: Edbb,Exbb,Ecbb
  double precision              :: Ed,Ex

! Output variables

  double precision,intent(out)  :: Ec(nsp)

! Hello world

  write(*,*)
  write(*,*)'********************************'
  write(*,*)'* Unrestricted MP2 Calculation *'
  write(*,*)'********************************'
  write(*,*)

!---------------------!
! Compute UMP2 energy |
!---------------------!

! aaaa block

  bra = 1
  ket = 1

  Edaa = 0d0
  Exaa = 0d0

  do i=nC(bra)+1,nO(bra)
    do a=nO(bra)+1,nBas-nR(bra)

      do j=nC(ket)+1,nO(ket)
        do b=nO(ket)+1,nBas-nR(ket)

          eps = eHF(i,bra) + eHF(j,ket) - eHF(a,bra) - eHF(b,ket) 
         
          Edaa = Edaa + 0.5d0*ERI_aa(i,j,a,b)*ERI_aa(i,j,a,b)/eps
          Exaa = Exaa - 0.5d0*ERI_aa(i,j,a,b)*ERI_aa(i,j,b,a)/eps


        end do
      end do
    end do
  end do

  Ecaa = Edaa + Exaa
  Ec(1) = Ecaa

! aabb block

  bra = 1
  ket = 2

  Edab = 0d0
  Exab = 0d0

  do i=nC(bra)+1,nO(bra)
    do a=nO(bra)+1,nBas-nR(bra)

      do j=nC(ket)+1,nO(ket)
        do b=nO(ket)+1,nBas-nR(ket)

          eps = eHF(i,bra) + eHF(j,ket) - eHF(a,bra) - eHF(b,ket) 
         
          Edab = Edab + ERI_ab(i,j,a,b)*ERI_ab(i,j,a,b)/eps

        end do
      end do
    end do
  end do

  Ecab = Edab + Exab
  Ec(2) = Ecab

! bbbb block

  bra = 2
  ket = 2

  Edbb = 0d0
  Exbb = 0d0

  do i=nC(bra)+1,nO(bra)
    do a=nO(bra)+1,nBas-nR(bra)

      do j=nC(ket)+1,nO(ket)
        do b=nO(ket)+1,nBas-nR(ket)

          eps = eHF(i,bra) + eHF(j,ket) - eHF(a,bra) - eHF(b,ket) 
         
          Edbb = Edbb + 0.5d0*ERI_bb(i,j,a,b)*ERI_bb(i,j,a,b)/eps
          Exbb = Exbb - 0.5d0*ERI_bb(i,j,a,b)*ERI_bb(i,j,b,a)/eps


        end do
      end do
    end do
  end do

  Ecbb = Edbb + Exbb
  Ec(3) = Ecbb

! Final flush

  Ed = Edaa + Edab + Edbb
  Ex = Exaa + Exab + Exbb

  write(*,*)
  write(*,'(A32)')           '--------------------------'
  write(*,'(A32)')           ' MP2 calculation          '
  write(*,'(A32)')           '--------------------------'
  write(*,'(A32,1X,F16.10)') ' MP2 correlation energy = ',sum(Ec(:))
  write(*,'(A32,1X,F16.10)') '   alpha-alpha          = ',Ecaa
  write(*,'(A32,1X,F16.10)') '   alpha-beta           = ',Ecab
  write(*,'(A32,1X,F16.10)') '    beta-beta           = ',Ecbb
  write(*,*)
  write(*,'(A32,1X,F16.10)') ' Direct part            = ',Ed
  write(*,'(A32,1X,F16.10)') '   alpha-alpha          = ',Edaa
  write(*,'(A32,1X,F16.10)') '   alpha-beta           = ',Edab
  write(*,'(A32,1X,F16.10)') '    beta-beta           = ',Edbb
  write(*,*)
  write(*,'(A32,1X,F16.10)') ' Exchange part          = ',Ex
  write(*,'(A32,1X,F16.10)') '   alpha-alpha          = ',Exaa
  write(*,'(A32,1X,F16.10)') '   alpha-beta           = ',Exab
  write(*,'(A32,1X,F16.10)') '    beta-beta           = ',Exbb
  write(*,'(A32)')           '--------------------------'
  write(*,'(A32,1X,F16.10)') ' MP2 electronic  energy = ',       EUHF + sum(Ec(:))
  write(*,'(A32,1X,F16.10)') ' MP2 total       energy = ',ENuc + EUHF + sum(Ec(:))
  write(*,'(A32)')           '--------------------------'
  write(*,*)

  if(dotest) then

    call dump_test_value('U','MP2 correlation energy',sum(Ec))

  end if

end subroutine 
