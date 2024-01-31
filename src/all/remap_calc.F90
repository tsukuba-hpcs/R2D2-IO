!========================================================================================|
!
!  R2D2 radiation MHD code
!  Copyright(C) Hideyuki Hotta <hotta.hieyuki@gmail.com>
!
!========================================================================================|
!----------------------------------------------------------------------------------|
! Japanese explanation is found at
! https://hottahd.github.io/R2D2-manual/notation.html#self-vc-dictionary
!----------------------------------------------------------------------------------|
subroutine remap_calc(qq, ix, jx, kx, ix0, jx0)
!----------------------------------------------------------------------------------|
!==================================================================================|
  use geometry_def, only: margin, xmin, xmax, ymin, ymax, zmin, zmax, mtype &
      & , x00, y00, z00 &
      & , ix00, jx00, kx00
  use time_def, only: cnd, cno,cnou

  use const_def, only: pi, pii8, pii4
#ifdef ideal
  use const_def, only: gm, cv, rr
#endif
  use, intrinsic :: ieee_arithmetic
  use remap_def, only: m2d_xy, m2d_xz, m2d_flux, m2d_spex, kkx0 &
       & , is0, js0, ir, jr, iixl, js,je &
       & , disp, fha, global_xy, global_xz, global_flux, global_spex &
       & , mpi_comm_world_remap &
       & , mpi_comm_world_remap_xz &
       & , cor, pdf2d, mtype_remap &
       & , mpi_comm_world_ir &
       & , mpi_comm_world_spex_io,jc,kc

#ifndef spherical
  use remap_def, only: jxr
#endif

#ifdef YinYang
  use remap_def, only: yt, zt
#endif

  
  use comm_def, only: myrank
  use background_def, only: &
       &  pr00, te00, ro00, se00, en00 &
       & , dprdro00, dprdse00, dtedro00, dtedse00, dendro00, dendse00 &
       & , gx00, dsedr00, dtedr00, cp00 &
       & , fa00

  use implicit_def

!$ use omp_lib

  implicit none
  integer :: ix, jx, kx
  real(KIND=KIND(0.e0)), dimension(ix, jx, kx, mtype_remap) :: qq
  real(KIND=KIND(0.e0)), dimension(ix, jx, kx) :: pr, te, op
  integer :: ix0, jx0
  integer :: ijx, ijx0, jkx, jkx0, ist0, iet0, jst0, jet0, kst0, ket0

  integer :: i1 = 1, i2 = 2
  integer :: j1 = 1, j2 = 2
  integer :: k1 = 1, k2 = 2

  logical exist
  include "mpif.h"
  integer, dimension(mpi_status_size) :: mstatus
  integer :: mmx

  ijx = ix*jx
  ijx0 = ix0*jx0

  jkx = jx*kx
  jkx0 = jx0*(kx - margin*2)
  
  ist0 = 1 + margin
  iet0 = ix - margin

  jst0 = 1 + margin
  jet0 = jx - margin

  kst0 = 1 + margin
  ket0 = kx - margin

  open (idf, file="data/remap/qq/"//cnou//"/"//cno//"/qq.dac."//cnd//"."//cno, form="unformatted", access="stream")
  write (idf) qq(ist0:iet0, jst0:jet0, kst0:ket0, 1:mtype) &
       &            , pr(ist0:iet0, jst0:jet0, kst0:ket0) &
       &            , te(ist0:iet0, jst0:jet0, kst0:ket0) &
       &            , op(ist0:iet0, jst0:jet0, kst0:ket0)
  close (idf)

return
end subroutine remap_calc
