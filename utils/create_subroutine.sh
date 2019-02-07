#!/bin/bash 

if [ $# != 1 ]
then
  echo "One argument required [name of subroutine]" 
fi
if [ $# = 1 ]
then

  NAME=$1

echo "subroutine ${NAME}()

! Description of the subroutine

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

end subroutine ${NAME}" > ${NAME}.f90

fi



