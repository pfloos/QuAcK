#!/bin/bash 

if [ $# != 2 ]
then
  echo "Two arguments required [name of function] [name of result]" 
fi
if [ $# = 2 ]
then

  NAME=$1

echo "function ${NAME}() result(${RES})

! Description of the function 

  implicit none

  include 'parameters.h'

! Input variables

  integer,intent(in)            :: 
  double precision,intent(in)   :: 

! Local variables

  integer                       :: 
  double precision              :: 

! Output variables

  integer,intent(out)           :: 
  double precision,intent(out)  :: 

! Initalization

end function " > ${NAME}.f90

fi



