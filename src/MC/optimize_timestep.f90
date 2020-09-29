subroutine optimize_timestep(nWalk,iMC,Acc,dt)

! Optimize dt to get 50% of accepted moves

  implicit none

! Input variables

  integer,intent(in)            :: nWalk,iMC
  double precision,intent(inout):: Acc,dt

! Local variables

  double precision              :: TotAcc,Current_Acc,Target_Acc,delta

  TotAcc = Acc/dble(nWalk)
  Current_Acc = 100d0*TotAcc/dble(iMC)

  Target_Acc  = 50.0d0

  delta = dt*abs(Target_Acc - Current_Acc)/100.d0
  if(Current_Acc > Target_Acc + 0.5d0)then
    dt = dt + delta
  elseif(Current_Acc < Target_Acc - 0.5d0)then
    dt = dt - delta
  endif

end subroutine optimize_timestep
