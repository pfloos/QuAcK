subroutine energy_from_rdm(ENuc,N,h,ERI_MO,rdm1,rdm2,Emu)

!input
integer,intent(in)               :: N
double precision                 :: ENuc
double precision,intent(in)      :: h(N,N),ERI_MO(N,N,N,N),rdm1(N,N),rdm2(N,N,N,N)

!local
integer                          ::p,q,r,s
double precision                 :: E1,E2

!output
double precision,intent(out)     :: Emu

E1 = 0d0
E2 = 0d0
Emu = 0d0

do p = 1, N
  do q = 1, N
    E1 = E1 + h(p,q) * rdm1(p,q)
  end do
end do

E2 = 0d0
do p = 1, N
  do q = 1, N
    do r = 1, N
      do s = 1, N
        E2 = E2 + eri_mo(p,q,r,s) * rdm2(p,q,r,s)
      end do
    end do
  end do
end do

Emu = E1 + 0.25*E2
write(*,'(A25,F16.10)') ' One-electron energy = ',E1
write(*,'(A25,F16.10)') ' Two-electron energy = ',E2
write(*,'(A25,F16.10)') ' Electronic   energy = ',Emu
write(*,*)
end subroutine
