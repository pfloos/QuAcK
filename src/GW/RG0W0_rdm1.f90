subroutine RG0W0_rdm1(nOrb,nS,lampl,rampl,rdm1)

! Compute 1-Reduced-Density-Matrix based in RG0W0

! Input
integer,intent(in)               :: nOrb,nS
double precision, intent(in)     :: lampl(nS,nOrb),rampl(nS,nOrb)

! Output
double precision,intent(out)     :: rdm1(nOrb,nOrb)

write(*,*) "1RDM Quack Quack 1RDM"
rdm1(:,:) = 0d0

end subroutine
