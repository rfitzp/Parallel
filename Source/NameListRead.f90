! NameListRead.f90

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Function to read PARALLEL namelist
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine NameListRead (MUE, LAMBDAD, SIGMA,&
     XMAX, NX,&
     GMAX, NG,&
     TMAX, NT,&
     GAMMA, KSMAX, FLG)&
     bind (c, name = 'NameListRead')

  use, intrinsic :: iso_c_binding, only: c_int, c_double
  implicit none

  real    (kind = c_double), intent (inout) :: MUE
  real    (kind = c_double), intent (inout) :: LAMBDAD
  real    (kind = c_double), intent (inout) :: SIGMA

  real    (kind = c_double), intent (inout) :: XMAX
  integer (kind = c_int),    intent (inout) :: NX

  real    (kind = c_double), intent (inout) :: GMAX
  integer (kind = c_int),    intent (inout) :: NG

  real    (kind = c_double), intent (inout) :: TMAX
  integer (kind = c_int),    intent (inout) :: NT

  real    (kind = c_double), intent (inout) :: GAMMA
  real    (kind = c_double), intent (inout) :: KSMAX

  integer (kind = c_int),    intent (inout) :: FLG
 
  namelist /PARALLEL_CONTROL/ MUE, LAMBDAD, SIGMA, XMAX, NX, GMAX, NG, TMAX, NT, GAMMA, KSMAX, FLG
       
  open  (unit = 100, file = 'Inputs/Namelist.nml', status = 'old')
  read  (unit = 100, nml  = PARALLEL_CONTROL)
  close (unit = 100)

endsubroutine namelistRead
