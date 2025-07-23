subroutine RG0W0_rdm2(nOrb,rdm2)

! Compute 2-Reduced-Density-Matrix based in RG0W0

! Input
integer,intent(in)               :: nOrb

! Output
double precision,intent(out)     :: rdm2(nOrb,nOrb,nOrb,nOrb)

write(*,*) "2RDM Quack Quack 2RDM"
rdm2(:,:,:,:) = 0d0 

end subroutine
