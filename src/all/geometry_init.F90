!========================================================================================|
!
!  R2D2 radiation MHD code
!  Copyright(C) Hideyuki Hotta <hotta.hieyuki@gmail.com>
!
!========================================================================================|
subroutine geometry_init
  use const_def, only: pi
  use geometry_def, only: ix00, jx00, kx00, ix0, i1, j1, k1 &
       & , xdcheck, ydcheck, zdcheck, marginx, marginy, marginz &
       & , xmax, xmin, ymax, ymin, zmax, zmin &
       & , nx, ny, nz &
       & , nxg, nyg, nzg &
       & , x00, y00, z00 &
       & , dx, dy, dz, dxi, dyi, dzi &
       & , dx0, dy0, dz0 &
       & , dxi2, dyi2, dzi2 &
       & , rkx1, rkx2, rkx3, rkx4 &
       & , rky1, rky2, rky3, rky4 &
       & , rkz1, rkz2, rkz3, rkz4 &
       & , rkx, rky, rkz &
       & , x, y, z, cflfacx, cflfacy, cflfacz &
       & , ix_ununi, dx00, ununiform_flag &
       & , dxyz_min &
       & , mpi_comm_world_top, mpi_comm_world_bot

#ifdef spherical
  use geometry_def, only: x2, x3, x4, xo, xo2, siny, siny2, siny3, sinyo, coty &
   & , rkx_x2xo2, rkx_xxo, rky_sinyosiny
   use mhd_def, only: ipa
#endif

  use comm_def, only: ib, jb, kb

  use implicit_def
  implicit none
  include "mpif.h"

   integer :: ip,jp,kp
  real(KIND(0.d0)) :: dx_unif, dy_unif, dz_unif
  real(KIND(0.d0)) :: fdx
  integer :: nxx
  real(KIND(0.d0)) :: dx11
  real(KIND(0.d0)) :: xrange, xrange0, xrange1

  ! for inhomogeneous grid
  !
  real(KIND(0.d0)) :: a, b, c, d, e, f, g, h
  real(KIND(0.d0)) :: dx_min, dy_min, dz_min
  real(KIND(0.d0)) :: dxyz_ming

  ! temporary variable for MPI_split
  integer :: color,key
!----------------------------------------------------------------------------------------|
  !
  ! define global geometry (Samples are for uniform)
  dx_unif = (xmax - xmin)/real(ix00 - 2*marginx)
  x00(1) = xmin + (0.5d0 - dble(marginx))*dx_unif

  if (xdcheck == 2) then
    do i = 1 + i1, ix00
      x00(i) = x00(i - i1) + dx_unif
    enddo
  endif

  dy_unif = (ymax - ymin)/real(jx00 - 2*marginy)
  y00(1) = ymin + (0.5d0 - dble(marginy))*dy_unif
  if (ydcheck == 2) then
    do j = 1 + j1, jx00
      y00(j) = y00(j - j1) + dy_unif
    enddo
  endif

  dz_unif = (zmax - zmin)/real(kx00 - 2*marginz)
  z00(1) = zmin + (0.5d0 - dble(marginz))*dz_unif
  if (zdcheck == 2) then
    do k = 1 + k1, kx00
      z00(k) = z00(k - k1) + dz_unif
    enddo
  endif

  if (ununiform_flag) then
!-----
    ! inhomogeneous spacing
    !
    ! homogeneous at top mpi thread
    ! This inhomogeneous setting width of the domain and number of grid
    ! Then the grid spacing increases downward linearly
    ! see an explanation in Powerpoint file provided by H. Hotta

    !nm = 1                     ! number of MPI thread for fine grid around top
    xrange = xmax - xmin         ! x domain size
!     xrange0 = dx00*dble(nx*nm) ! x domain size in top MPI thread, in which grid is uniform
    xrange0 = dx00*dble(ix_ununi) ! x domain size in top MPI thread, in which grid is uniform
    xrange1 = xrange - xrange0 ! x domain size in the other MPI thread than top
!     nxx = (ix0 - nm)*nx        ! number of grid for inhomogeneous grid
    nxx = ix0*nx - ix_ununi        ! number of grid for inhomogeneous grid

    fdx = 2.d0*(xrange1 - dx00*dble(nxx))/dble((nxx - 4)*nxx) ! increase rate of the grid

    !-------------------------------------
    ! grid for uniform area
    ! set uniform area
    ! set the geometry at the top margin

    ! set the geometry around the top boundary
    ! the geometry is defined downward from the top
    x00(ix00 - marginx) = xmax - 0.5d0*dx00

    do i = ix00 - marginx + 1, ix00
      x00(i) = x00(i - 1) + dx00
    enddo

    do i = ix00 - marginx - 1, ix00 - marginx - ix_ununi + 1, -1
      x00(i) = x00(i + 1) - dx00
    enddo

    ! uniform area finish
    !-------------------------------------

    ! first 2 uniform grids in ununiform area
    do i = ix00 - marginx - ix_ununi, ix00 - marginx - ix_ununi - 1, -1
      x00(i) = x00(i + 1) - dx00
    enddo

    ! ordinary ununiform gridding
    dx11 = dx00
    do i = ix00 - marginx - ix_ununi - 2, 5, -1
      x00(i) = x00(i + 1) - dx11
      dx11 = dx11 + fdx
    enddo

    ! final 4 uniform grid including margin
    do i = 4, 1, -1
      x00(i) = x00(i + 1) - dx11
    enddo
  endif

!=======================================================================================|
! transformation from global to local

  ! Radius or x
  do i = 1, nxg
    x(i) = x00(ib*nx + i)
  enddo

  ! Theta or y
  do j = 1, nyg
    y(j) = y00(jb*ny + j)
  enddo

  ! Phi or z
  do k = 1, nzg
    z(k) = z00(kb*nz + k)
  enddo

!----------------------------------------------------------------------------------------|
  ! grid spacing

  !============
  ! x direction
  if (xdcheck == 2) then
    do i = 1 + marginx, nxg - marginx
      a = x(i + 2) - x(i)
      b = x(i + 1) - x(i)
      c = x(i) - x(i - 1)
      d = x(i) - x(i - 2)
      rkx1(i) = (c*d - b*(c + d))/(a - b)/(a + c)/(a + d)
      rkx2(i) = (a*(c + d) - c*d)/(a - b)/(b + c)/(b + d)
      rkx3(i) = (b*d + a*(d - b))/(a + c)/(b + c)/(c - d)
      rkx4(i) = (a*(b - c) - b*c)/(c - d)/(a + d)/(b + d)

      rkx(i, +2) = rkx1(i)
      rkx(i, +1) = rkx2(i)
      rkx(i, -1) = rkx3(i)
      rkx(i, -2) = rkx4(i)
    enddo

    do i = 1+i1,nxg
      dx0 = x(i) - x(i - 1)
      dx(i) = dx0
      dxi(i) = 1.d0/(x(i) - x(i - 1))
    enddo

    do i = 1 + i1, nxg - i1
      dxi2(i) = 1.d0/(x(i + i1) - x(i - i1))
    enddo
  else

    dx0 = 1.d20
    dx(1) = 1.d20
    dxi(1) = 0.d0
    dxi2(1) = 0.d0
    rkx1(1) = 0.d0
    rkx2(1) = 0.d0
    rkx3(1) = 0.d0
    rkx4(1) = 0.d0

    rkx = 0.d0
  endif

  !============
  ! y direction
  if (ydcheck == 2) then
    do j = 1 + marginy, nyg - marginy
      a = y(j + 2) - y(j)
      b = y(j + 1) - y(j)
      c = y(j) - y(j - 1)
      d = y(j) - y(j - 2)
      rky1(j) = (c*d - b*(c + d))/(a - b)/(a + c)/(a + d)
      rky2(j) = (a*(c + d) - c*d)/(a - b)/(b + c)/(b + d)
      rky3(j) = (b*d + a*(d - b))/(a + c)/(b + c)/(c - d)
      rky4(j) = (a*(b - c) - b*c)/(c - d)/(a + d)/(b + d)

      rky(j, +2) = rky1(j)
      rky(j, +1) = rky2(j)
      rky(j, -1) = rky3(j)
      rky(j, -2) = rky4(j)

      dy0 = y(j) - y(j - 1)
      dy(j) = dy0
      dyi(j) = 1.d0/dy(j)

      dy_min = min(dy_min, dy0)
    enddo

    do j = 1+j1,nyg
      dy0 = y(j) - y(j - 1)
      dy(j) = dy0
      dyi(j) = 1.d0/dy(j)
    enddo

    do j = 1 + j1, nyg - j1
      dyi2(j) = 1.d0/(y(j + j1) - y(j - j1))
    enddo
  else

    dyi(1) = 0.d0
    dy(1) = 1.d20
    dyi(1) = 0.d0
    dyi2(1) = 0.d0
    rky1(1) = 0.d0
    rky2(1) = 0.d0
    rky3(1) = 0.d0
    rky4(1) = 0.d0

    rky = 0.d0
  endif

  !============
  ! z direction
  if (zdcheck == 2) then
    do k = 1 + marginz, nzg - marginz
      a = z(k + 2) - z(k)
      b = z(k + 1) - z(k)
      c = z(k) - z(k - 1)
      d = z(k) - z(k - 2)
      rkz1(k) = (c*d - b*(c + d))/(a - b)/(a + c)/(a + d)
      rkz2(k) = (a*(c + d) - c*d)/(a - b)/(b + c)/(b + d)
      rkz3(k) = (b*d + a*(d - b))/(a + c)/(b + c)/(c - d)
      rkz4(k) = (a*(b - c) - b*c)/(c - d)/(a + d)/(b + d)

      rkz(k, +2) = rkz1(k)
      rkz(k, +1) = rkz2(k)
      rkz(k, -1) = rkz3(k)
      rkz(k, -2) = rkz4(k)

      dz0 = z(k) - z(k - 1)
      dz(k) = dz0
      dzi(k) = 1.d0/dz(k)

      dz_min = min(dz_min, dz0)
    enddo

    do k = 1+k1,nzg
      dz0 = z(k) - z(k - 1)
      dz(k) = dz0
      dzi(k) = 1.d0/dz(k)
    enddo

    do k = 1 + k1, nzg - k1
      dzi2(k) = 1.d0/(z(k + k1) - z(k - k1))
    enddo
  else

    dz(1) = 1.d20
    dzi(1) = 0.d0
    dzi2(1) = 0.d0
    rkz1(1) = 0.d0
    rkz2(1) = 0.d0
    rkz3(1) = 0.d0
    rkz4(1) = 0.d0

    rkz = 0.d0
  endif

  !----------------------------------------------------------------------------------------|
  ! for CFL condition
#ifdef spherical
  do i = 1, nxg
    xo(i) = 1.d0/x(i)
    x2(i) = x(i)**2
    x3(i) = x(i)**3
    x4(i) = x(i)**4
    xo2(i) = 1.d0/x2(i)
    do l = 1,2*marginx
    enddo
  enddo

  do l = 1,2*marginx
    ip = ipa(l)
    do i = 1+marginx,nxg-marginx
      rkx_x2xo2(i,ip) = rkx(i,ip)*x2(i+ip)*xo2(i)
      rkx_xxo(i,ip) = rkx(i,ip)*x(i+ip)*xo(i)
    enddo
  enddo 

  do j = 1, nyg
    siny(j) = sin(y(j))
    siny2(j) = sin(y(j))**2
    siny3(j) = sin(y(j))**3
    sinyo(j) = 1.d0/siny(j)
    coty(j) = 1.d0/tan(y(j))
  enddo

  do l = 1,2*marginy
    jp = ipa(l)
    do j = 1+marginy,nyg-marginy
      rky_sinyosiny(j,jp) = rky(j,jp)*sinyo(j)*siny(j+jp)
    enddo
  enddo

#endif

  !----------------------------------------------------------------------------------------|
  ! for CFL condition

  cflfacx = 1.d0
  cflfacy = 1.d0
  cflfacz = 1.d0

  cflfacx(1:marginx) = 1.d20
  cflfacx(nxg - marginx + 1:nxg) = 1.d20

  cflfacy(1:marginy) = 1.d20
  cflfacy(nyg - marginy + 1:nyg) = 1.d20

  cflfacz(1:marginz) = 1.d20
  cflfacz(nzg - marginz + 1:nzg) = 1.d20

  !----------------------------------------------------------------------------------------|
  ! for divb cleaning
  ! smallest grid is evaluate

  dx_min = 1.d30
  do i = 1 + marginx, nxg - marginx
    dx_min = min(dx_min, dx(i))
  enddo

#ifdef spherical
  dy_min = 1.d30
  do j = 1 + marginy, nyg - marginy
  do i = 1 + marginx, nxg - marginx
    dy_min = min(dy_min, x(i)*dy(j))
  enddo
  enddo

  dz_min = 1.d30
  do k = 1 + marginz, nzg - marginz
  do j = 1 + marginy, nyg - marginy
  do i = 1 + marginx, nxg - marginx
    dz_min = min(dz_min, x(i)*sin(y(j))*dz(k))
  enddo
  enddo
  enddo
#else
  dy_min = 1.d30
  do j = 1 + marginy, nyg - marginy
    dy_min = min(dy_min, dy(j))
  enddo

  dz_min = 1.d30
  do k = 1 + marginz, nzg - marginz
    dz_min = min(dz_min, dz(k))
  enddo
#endif

  dxyz_min = min(dx_min, dy_min, dz_min)
  call mpi_allreduce(dxyz_min, dxyz_ming, 1, mpi_real8, mpi_min, mpi_comm_world, merr)
  dxyz_min = dxyz_ming

!=================================================================================
! Create MPI split

  !===================
  ! MPI world for only top thread
  if (ib /= ix0 - 1) then
    color = 1
  else
    color = 0
  endif
! split MPI_COMM_WORLD
  call mpi_comm_split(mpi_comm_world, color, key, mpi_comm_world_top, merr)

  !===================
  if (ib /= 0) then
    color = 1
  else
    color = 0
  endif
! split MPI_COMM_WORLD
  call mpi_comm_split(mpi_comm_world, color, key, mpi_comm_world_bot, merr)

  return
end subroutine geometry_init
