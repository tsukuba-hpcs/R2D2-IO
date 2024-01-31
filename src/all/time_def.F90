!========================================================================================|
!
!  R2D2 radiation MHD code
!  Copyright(C) Hideyuki Hotta <hotta.hieyuki@gmail.com>
!  
!========================================================================================|
module time_def

  character(LEN=8), save :: cnd
  character(LEN=8), save :: cnd_tau
  character(LEN=8), save :: cno
  character(LEN=5), save :: cnou
  character(LEN=8), save ::  cj
  integer, save :: mcont,nd,nd_tau,nd0,ns
  real(KIND(0.d0)), save :: dt,time,timep

  ! initial time
  real(KIND(0.e0)), save :: time00

end module time_def
