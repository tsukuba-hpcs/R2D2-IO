!========================================================================================|
!
!  R2D2 radiation MHD code
!  Copyright(C) Hideyuki Hotta <hotta.hieyuki@gmail.com>
!  
!========================================================================================|
module yinyang_def
  use geometry_def, only: nyg,nzg,yinyang_t,npe0
#ifdef YinYang
  
  real(KIND(0.d0)), dimension(nyg,nzg), save :: yo,zo
  real(KIND(0.d0)), dimension(nyg,nzg), save :: coss,sins

  ! for MPI split
  integer, save :: mpi_comm_world_yy,myrank_yy,npe_yy

  ! for Yin-Yang MPI
  integer, parameter :: four = 4
  integer, parameter :: two = 2
  integer, parameter :: yinyang_parameter = 8
  integer, dimension(0:npe0-1), save :: yinyangcount_rcv,yinyangcount_snd
  real(KIND(0.d0)), dimension(four,yinyang_t,0:npe0-1), save :: yinyang_yz_snd,yinyang_yz_rcv
  integer, dimension(two,yinyang_t,0:npe0-1), save :: yinyang_jk_snd,yinyang_jk_rcv
  real(KIND(0.d0)), dimension(yinyang_parameter,yinyang_t,0:npe0-1), save :: yinyang_dyz

#endif  
real(KIND(0.d0)), dimension(nyg,nzg), save :: yzmask = 1.d0
real(KIND(0.d0)), dimension(nyg,nzg), save :: yzmask_act = 1.d0

end module yinyang_def
