subroutine MP2F12(nBas,nC,nO,nV,nA,ERI,F12,Yuk,EHF,EcMP2,c,cA,EcMP2F12)

! Perform restricted Hartree-Fock calculation

  implicit none

! Input variables

  integer,intent(in)            :: nBas,nC,nO,nV,nA
  double precision,intent(in)   :: EHF,EcMP2
  double precision,intent(in)   :: c(nBas,nBas),cA(nBas,nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas),F12(nBas,nBas,nBas,nBas),Yuk(nBas,nBas,nBas,nBas)

! Local variables

  double precision,allocatable  :: ooCoo(:,:,:,:),ooFoo(:,:,:,:),ooYoo(:,:,:,:)
  double precision,allocatable  :: ooCoa(:,:,:,:),ooFoa(:,:,:,:)
  double precision,allocatable  :: ooCvv(:,:,:,:),ooFvv(:,:,:,:)
  double precision,allocatable  :: cO(:,:),cV(:,:)
  double precision              :: E2a,E2b,E3a,E3b,E4a,E4b,E4c,E4d
  integer                       :: i,j,k,l,a,b,x

! Output variables

  double precision,intent(out)  :: EcMP2F12(3)

! Split MOs into occupied and virtual sets

  allocate(cO(nBas,nO),cV(nBas,nV))
  cO(1:nBas,1:nO) = c(1:nBas,nC+1:nC+nO)
  cV(1:nBas,1:nV) = c(1:nBas,nC+nO+1:nBas)

! Compute the two-electron part of the MP2-F12 energy

  allocate(ooYoo(nO,nO,nO,nO))
  call AOtoMO_oooo(nBas,nO,cO,Yuk,ooYoo)

  E2a = 0d0
  E2b = 0d0
  do i=1,nO
    do j=1,nO
      E2a = E2a + ooYoo(i,j,i,j)
      E2b = E2b + ooYoo(i,j,j,i)
    enddo
  enddo

  deallocate(ooYoo)

! Compute the three-electron part of the MP2-F12 energy

  allocate(ooCoa(nO,nO,nO,nA),ooFoa(nO,nO,nO,nA))
  call AOtoMO_oooa(nBas,nO,nA,cO,cA,ERI,ooCoa)
  call AOtoMO_oooa(nBas,nO,nA,cO,cA,F12,ooFoa)

  E3a = 0d0
  E3b = 0d0
  do i=1,nO
    do j=1,nO
      do k=1,nO
        do x=1,nA
          E3a = E3a + ooCoa(i,j,k,x)*ooFoa(j,i,k,x)
          E3b = E3b + ooCoa(i,j,k,x)*ooFoa(i,j,k,x)
        enddo
      enddo
    enddo
  enddo

  deallocate(ooCoa,ooFoa)

! Compute the four-electron part of the MP2-F12 energy

  allocate(ooCoo(nO,nO,nO,nO),ooFoo(nO,nO,nO,nO))
  call AOtoMO_oooo(nBas,nO,cO,ERI,ooCoo)
  call AOtoMO_oooo(nBas,nO,cO,F12,ooFoo)

  E4a = 0d0
  E4b = 0d0
  do i=1,nO
    do j=1,nO
      do k=1,nO
        do l=1,nO
          E4a = E4a + ooCoo(i,j,k,l)*ooFoo(i,j,k,l)
          E4b = E4b + ooCoo(i,j,k,l)*ooFoo(j,i,k,l)
        enddo
      enddo
    enddo
  enddo

  deallocate(ooCoo,ooFoo)

  allocate(ooCvv(nO,nO,nV,nV),ooFvv(nO,nO,nV,nV))
  call AOtoMO_oovv(nBas,nO,nV,cO,cV,ERI,ooCvv)
  call AOtoMO_oovv(nBas,nO,nV,cO,cV,F12,ooFvv)

  E4c = 0d0
  E4d = 0d0
  do i=1,nO
    do j=1,nO
      do a=1,nV
        do b=1,nV
          E4c = E4c + ooCvv(i,j,a,b)*ooFvv(i,j,a,b)
          E4d = E4d + ooCvv(i,j,a,b)*ooFvv(j,i,a,b)
        enddo
      enddo
    enddo
  enddo

  deallocate(ooCvv,ooFvv)

! Final scaling of the various components

  EcMP2F12(1) = +0.625d0*E2a - 0.125d0*E2b
  EcMP2F12(2) = -1.250d0*E3a + 0.250d0*E3b
  EcMP2F12(3) = +0.625d0*E4a - 0.125d0*E4b - 0.625d0*E4c + 0.125d0*E4d

  write(*,*)
  write(*,'(A32)')           '-----------------------'
  write(*,'(A32)')           ' MP2-F12 calculation   '
  write(*,'(A32)')           '-----------------------'
  write(*,'(A32,1X,F16.10)') ' MP2                   ',EcMP2
  write(*,'(A32,1X,F16.10)') ' MP2-F12 E(2)          ',EcMP2F12(1)
  write(*,'(A32,1X,F16.10)') ' MP2-F12 E(3)          ',EcMP2F12(2)
  write(*,'(A32,1X,F16.10)') ' MP2-F12 E(4)          ',EcMP2F12(3)
  write(*,'(A32)')           '-----------------------'
  write(*,'(A32,1X,F16.10)') ' Total                 ',EcMP2+EcMP2F12(1)+EcMP2F12(2)+EcMP2F12(3)
  write(*,'(A32)')           '-----------------------'
  write(*,*)

  deallocate(cO,cV)

end subroutine MP2F12
