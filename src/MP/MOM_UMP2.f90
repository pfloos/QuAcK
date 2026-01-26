subroutine MOM_UMP2(dotest,nBas,nC,nO,nV,nR,nCVS,occupations,ERI_aa,ERI_ab,ERI_bb,ENuc,EUHF,eHF,Ec)

! Perform unrestricted second-order Moller-Plesset calculation on top of a MOM Reference

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dotest

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC(nspin)
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nV(nspin)
  integer,intent(in)            :: nR(nspin)
  integer,intent(in)            :: nCVS(nspin)
  integer,intent(in)            :: occupations(maxval(nO),nspin)
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: EUHF
  double precision,intent(in)   :: ERI_aa(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_ab(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_bb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: eHF(nBas,nspin)

! Local variables

  integer                       :: bra,ket
  integer                       :: i,j,a,b,ispin
  double precision              :: eps
  double precision              :: Edaa,Exaa,Ecaa
  double precision              :: Edab,Exab,Ecab
  double precision              :: Edbb,Exbb,Ecbb
  double precision              :: Ed,Ex
  integer,allocatable           :: virtuals(:,:)

! Output variables

  double precision,intent(out)  :: Ec(nsp)

! Hello world

  write(*,*)
  write(*,*)'********************************'
  write(*,*)'* Unrestricted MP2 Calculation *'
  write(*,*)'********************************'
  write(*,*)


! Check no frozen core
  if(any(nC/=0)) then
    print *, "Do not use PyDuck frozen core with CVS !"
    stop
  end if

! CVS
  if(any(nCVS>0)) then
    print *, "No exications to the first", nCVS(1), "alpha orbital(s) are considered."
    print *, "No exications to the first", nCVS(2), "beta orbital(s) are considered."
  end if

! Get virtuals
  allocate(virtuals(nBas-minval(nO),nspin))
  virtuals = 0
  do ispin=1,nspin
    call non_occupied(nO(ispin),nBas,occupations(1:nO(ispin),ispin),virtuals(1:nBas - nO(ispin),ispin))
  end do

!---------------------!
! Compute UMP2 energy |
!---------------------!

! aaaa block

  bra = 1
  ket = 1

  Edaa = 0d0
  Exaa = 0d0

  do i=1,nO(bra)
    do a=1+nCVS(bra),nBas - nO(bra)

      do j=nCVS(ket)+1,nO(ket)
        do b=1,nBas-nO(ket)

          eps = eHF(occupations(i,bra),bra) + eHF(occupations(j,ket),ket) - eHF(virtuals(a,bra),bra) - eHF(virtuals(b,ket),ket) 
         
          Edaa = Edaa& 
               + 0.5d0*ERI_aa(occupations(i,bra),occupations(j,ket),virtuals(a,bra),virtuals(b,ket))&
               * ERI_aa(occupations(i,bra),occupations(j,ket),virtuals(a,bra),virtuals(b,ket))&
               / eps 
          Exaa = Exaa &
               - 0.5d0*ERI_aa(occupations(i,bra),occupations(j,ket),virtuals(a,bra),virtuals(b,ket))&
               * ERI_aa(occupations(i,bra),occupations(j,ket),virtuals(b,bra),virtuals(a,ket))&
               / eps


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

  do i=1,nO(bra)
    do a=nCVS(bra)+1,nBas-nO(bra)

      do j=1,nO(ket)
        do b=nCVS(ket)+1,nBas-nO(ket)

          eps = eHF(occupations(i,bra),bra) + eHF(occupations(j,ket),ket) - eHF(virtuals(a,bra),bra) - eHF(virtuals(b,ket),ket) 
         
          Edab = Edab                                                                                    &
               + ERI_ab(occupations(i,bra),occupations(j,ket),virtuals(a,bra),virtuals(b,ket))           &
               * ERI_ab(occupations(i,bra),occupations(j,ket),virtuals(a,bra),virtuals(b,ket))           &
               / eps

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

  do i=1,nO(bra)
    do a=nCVS(bra)+1,nBas-nO(bra)

      do j=1,nO(ket)
        do b=nCVS(ket)+1,nBas-nO(ket)

          eps = eHF(occupations(i,bra),bra) + eHF(occupations(j,ket),ket) - eHF(virtuals(a,bra),bra) - eHF(virtuals(b,ket),ket) 
         
          Edbb = Edbb& 
               + 0.5d0*ERI_bb(occupations(i,bra),occupations(j,ket),virtuals(a,bra),virtuals(b,ket))    &
               * ERI_bb(occupations(i,bra),occupations(j,ket),virtuals(a,bra),virtuals(b,ket))          &
               / eps
          Exbb = Exbb&
               - 0.5d0*ERI_bb(occupations(i,bra),occupations(j,ket),virtuals(a,bra),virtuals(b,ket))&
               * ERI_bb(occupations(i,bra),occupations(j,ket),virtuals(b,ket),virtuals(a,bra))&
               / eps

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
