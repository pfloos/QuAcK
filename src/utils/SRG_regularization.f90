double precision function SRG_reg(dem,flow)

  implicit none

  double precision,intent(in)  :: dem,flow
  double precision             :: reg,eps

  eps = 1d-12

  if(abs(dem) < eps) then
    reg = 2d0*flow*dem
  else
    reg = (1d0 - exp(-2d0*flow*dem*dem))/dem
  end if

  SRG_reg = reg

end function 

double precision function SRG_reg2(dem1,dem2,flow)

  implicit none

  double precision,intent(in)  :: dem1,dem2,flow
  double precision             :: reg1,reg2,eps

  eps = 1d-12

  if(abs(dem1) < eps) then
    reg1 = 2d0*flow*dem1
  else
    reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))/dem1
  end if

  if(abs(dem2) < eps) then
    reg2 = 2d0*flow*dem2
  else
    reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))/dem2
  end if

  SRG_reg2 = reg1*reg2

end function

double precision function SRG_reg3(dem1,dem2,dem3,flow)

  implicit none

  double precision,intent(in)  :: dem1,dem2,dem3,flow
  double precision             :: reg1,reg2,reg3,eps

  eps = 1d-12

  if(abs(dem1) < eps) then
    reg1 = 2d0*flow*dem1
  else
    reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))/dem1
  end if

  if(abs(dem2) < eps) then
    reg2 = 2d0*flow*dem2
  else
    reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))/dem2
  end if

  if(abs(dem3) < eps) then
    reg3 = 2d0*flow*dem3
  else
    reg3 = (1d0 - exp(-2d0*flow*dem3*dem3))/dem3
  end if

  SRG_reg3 = reg1*reg2*reg3

end function 
