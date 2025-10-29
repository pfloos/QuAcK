#!/bin/bash

name=`git rev-parse HEAD`

echo "subroutine gitversion(sha)" > gitver.f90
echo "!Arguments ------------------------------------ " >> gitver.f90
echo "!scalars " >> gitver.f90
echo "!arrays " >> gitver.f90
echo "character(100),intent(inout)::sha" >> gitver.f90  >> gitver.f90
echo "!Local variables ------------------------------ " >> gitver.f90
echo "!scalars " >> gitver.f90
echo "!arrays " >> gitver.f90
echo "!************************************************************************ " >> gitver.f90
echo " " >> gitver.f90
echo "  write(sha,'(a)')'$name'" >> gitver.f90
echo " " >> gitver.f90
echo "end subroutine gitversion" >> gitver.f90
