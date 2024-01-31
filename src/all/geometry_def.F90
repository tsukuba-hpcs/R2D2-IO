!========================================================================================|
!
!  R2D2 radiation MHD code
!  Copyright(C) Hideyuki Hotta <hotta.hieyuki@gmail.com>
!
!========================================================================================|
! Define parameters for 3D full radiation transfer
module geometry_def
  use star_def, only: rstar, lstar
  use const_def, only: pi

  integer, parameter :: margin = 2 ! number of margin
  real(KIND(0.d0)), parameter :: n_int = dble(margin + 1) ! interpolation for Yin-Yang

  ! number of grid in a MPI thread
  ! both ny0/kx0 and nz0/jx0 must be integer for bc_potential.F90
  !integer, parameter, private :: nx0 = 64, ny0 = 96, nz0 = 144
  integer, parameter, private :: nx0 = 256, ny0 = 128, nz0 = 128

  ! number of MPI thread for each directions
  integer, parameter :: ix0 = 8, jx0 = 16, kx0 = 16
  integer, parameter, private :: one = 1

  ! check for dimension
  integer, parameter :: xdcheck = 2
  integer, parameter :: ydcheck = 2
  integer, parameter :: zdcheck = 2
  integer, parameter, private :: idd = xdcheck - 1
  integer, parameter, private :: jdd = ydcheck - 1
  integer, parameter, private :: kdd = zdcheck - 1

  ! The number of grid in a block (i.e., an MPI thread)
  ! dummy value for nx, ny, nz
 
  ! for changing of dimension
  integer, parameter :: nx = nx0*idd + abs(idd - 1)
  integer, parameter :: ny = ny0*jdd + abs(jdd - 1)
  integer, parameter :: nz = nz0*kdd + abs(kdd - 1)
  ! When the value is 1, the program becomes 1D or 2D code
  integer, parameter :: i1 = idd, i2 = idd*2
  integer, parameter :: j1 = jdd, j2 = jdd*2
  integer, parameter :: k1 = kdd, k2 = kdd*2
  integer, parameter :: marginx = margin*idd
  integer, parameter :: marginy = margin*jdd
  integer, parameter :: marginz = margin*kdd
  integer, parameter :: nxg = nx + marginx*2
  integer, parameter :: nyg = ny + marginy*2
  integer, parameter :: nzg = nz + marginz*2

  integer, parameter :: yinyang_t = (nyg + nzg)*margin*2

  ! for better OpenMP parallalization
  integer, parameter :: nyz = ny*nz
  integer, parameter :: nyzg = nyg*nzg

#ifdef YinYang
  integer, parameter :: yinyang = 2
#else
  integer, parameter :: yinyang = 1
#endif

! calculation domain
#ifdef ideal
!==================================
! setting for idealized calculation

  real(KIND(0.d0)), parameter :: xmax = 0.d0
  real(KIND(0.d0)), parameter :: xmin = -1.d0
#ifdef spherical
  ! Spherical geometry (ideal setting)
  ! ideal calculation
  real(KIND(0.d0)), parameter :: ymin =  30.d0/180.d0*pi
  real(KIND(0.d0)), parameter :: ymax = 150.d0/180.d0*pi
  real(KIND(0.d0)), parameter :: zmin =   0.d0/180.d0*pi
  real(KIND(0.d0)), parameter :: zmax = 120.d0/180.d0*pi
#else
  !! Cartesian geometry (ideal setting)
  real(KIND(0.d0)), parameter :: ymax = 4.d0
  real(KIND(0.d0)), parameter :: ymin = 0.d0
  real(KIND(0.d0)), parameter :: zmax = 4.d0
  real(KIND(0.d0)), parameter :: zmin = 0.d0
#endif

#else

!==================================
!==================================
!==================================
! setting for non ideal (solar and stellar) calculation

#ifdef deep
  real(KIND(0.d0)), parameter :: xmax = 0.96d0*rstar
  real(KIND(0.d0)), parameter :: xmin = 0.71d0*rstar
#else
  real(KIND(0.d0)), parameter :: xmax = rstar + 94.1d5
  !real(KIND(0.d0)), parameter :: xmin = rstar - 17.732d8
  real(KIND(0.d0)), parameter :: xmin = rstar - 494.3d5
  !real(KIND(0.d0)), parameter :: xmin = rstar - 11.588d8
  !real(KIND(0.d0)), parameter :: xmin = 0.71d0*rstar
#endif

!=====================
! setting for spherical geometry start
#ifdef spherical
#ifdef YinYang
  ! calculate the Yin-Yang grid
  real(KIND(0.d0)), parameter, private :: ymax_tmp = &
       & 2.d0*(dble(jx0*ny) - n_int + 0.5)/dble(jx0*ny)
  real(KIND(0.d0)), parameter :: ymax = (3.d0*pi/2.d0 &
       & - (n_int - 0.5d0)/dble(jx0*ny)*pi)/ymax_tmp
  real(KIND(0.d0)), parameter :: ymin = pi - ymax
  real(KIND(0.d0)), parameter, private :: dy_tmp = (ymax - ymin)/dble(jx0*ny)

  real(KIND(0.d0)), parameter, private :: zmax_tmp0 = 1.d0 + 1.d0/dble(kx0*nz)
  real(KIND(0.d0)), parameter, private :: zmax_tmp1 = 1.d0 - 2.d0*n_int/dble(kx0*nz)
  real(KIND(0.d0)), parameter :: zmax0 = (3.d0*pi/2.d0 - ymax + n_int*dy_tmp)/zmax_tmp0
  real(KIND(0.d0)), parameter :: zmax1 = (3.d0*pi/2.d0 - ymax - 0.5*dy_tmp)/zmax_tmp1

  real(KIND(0.d0)), parameter :: zmax = max(zmax0, zmax1)
  real(KIND(0.d0)), parameter :: zmin = -zmax
#else
  ! large domain
  real(KIND(0.d0)), parameter :: ymin = (90.d0 - 60.d0)/180.d0*pi
  real(KIND(0.d0)), parameter :: ymax = (90.d0 + 60.d0)/180.d0*pi
  real(KIND(0.d0)), parameter :: zmin = 0.d0
  real(KIND(0.d0)), parameter :: zmax = 120.d0/180.d0*pi

  ! small domain  
!  real(KIND(0.d0)), parameter :: ymin = (90.d0 - 16.d0)/180.d0*pi
!  real(KIND(0.d0)), parameter :: ymax = (90.d0 + 16.d0)/180.d0*pi
!  real(KIND(0.d0)), parameter :: zmin = 0.d0
!  real(KIND(0.d0)), parameter :: zmax = 32.d0/180.d0*pi
  
#endif
!=====================
! setting for cartesian geometry start
#else

#ifdef deep
  ! for deep CZ calculation
  real(KIND(0.d0)), parameter :: ymin = 0.d0
  real(KIND(0.d0)), parameter :: ymax = 6.144d8*8.d0
  real(KIND(0.d0)), parameter :: zmin = 0.d0
  real(KIND(0.d0)), parameter :: zmax = 6.144d8*8.d0
#else
  ! for surface calculation
  real(KIND(0.d0)), parameter :: ymin = 0.d0
  real(KIND(0.d0)), parameter :: ymax = 770.3d5
  real(KIND(0.d0)), parameter :: zmin = 0.d0
  real(KIND(0.d0)), parameter :: zmax = 770.3d5
#endif

#endif

!==================================
!==================================
!==================================
#endif
  ! for ununiform grid
#ifdef deep
  logical, parameter :: ununiform_flag = .false.
  integer, parameter :: ix_ununi = 0
  ! uniform grid for ununiformgird_setting
  real(KIND(0.d0)) :: dx00 = 2.d7
#else
  logical, parameter :: ununiform_flag = .false.
  integer, parameter :: ix_ununi = 96
  ! uniform grid for ununiformgird_setting
  real(KIND(0.d0)) :: dx00 = 48.d5
#endif

  integer, parameter :: mtype = 9 ! the number of the kinds of variable

  ! The number of blocks each parent block have
#ifdef spherical

#ifdef YinYang
  integer, parameter :: iper = 0, jper = 0, kper = 0
#else
  integer, parameter :: iper = 0, jper = 0, kper = 1
#endif

#else
  integer, parameter :: iper = 0, jper = 1, kper = 1
#endif

  ! flag for periodic boundary condition 0:without 1:with

  integer, parameter :: ix00 = ix0*nx + 2*marginx
  integer, parameter :: kx00 = kx0*nz + 2*marginz
  integer, parameter :: jx00 = jx0*ny + 2*marginy

  integer, parameter :: npe0 = ix0*jx0*kx0*yinyang

  ! data in global region
  real(KIND(0.d0)), dimension(ix00), save :: x00
  real(KIND(0.d0)), dimension(jx00), save :: y00
  real(KIND(0.d0)), dimension(kx00), save :: z00

  ! data in local region(each MPI thread)
  real(KIND(0.d0)), dimension(nxg), save :: x, dx, dxi, dxi2
  real(KIND(0.d0)), dimension(nyg), save :: y, dy, dyi, dyi2
  real(KIND(0.d0)), dimension(nzg), save :: z, dz, dzi, dzi2
#ifdef spherical
  real(KIND(0.d0)), dimension(nxg), save :: x2, x3, x4, xo, xo2
  real(KIND(0.d0)), dimension(nyg), save :: siny, coty, sinyo, siny2, siny3
#endif

  real(KIND(0.d0)), save :: dx0, dy0, dz0
  ! inhomogeneous Runge-Kutta Factor
  real(KIND(0.d0)), dimension(nxg), save :: rkx1, rkx2, rkx3, rkx4
  real(KIND(0.d0)), dimension(nyg), save :: rky1, rky2, rky3, rky4
  real(KIND(0.d0)), dimension(nzg), save :: rkz1, rkz2, rkz3, rkz4

  real(KIND(0.d0)), dimension(nxg, -margin:margin), save :: rkx
  real(KIND(0.d0)), dimension(nyg, -margin:margin), save :: rky
  real(KIND(0.d0)), dimension(nzg, -margin:margin), save :: rkz

#ifdef spherical
  real(KIND(0.d0)), dimension(nxg, -margin:margin), save :: rkx_x2xo2, rkx_xxo
  real(KIND(0.d0)), dimension(nyg, -margin:margin), save :: rky_sinyosiny
#endif

  ! for CFL condition
  real(KIND(0.d0)), dimension(nxg), save :: cflfacx
  real(KIND(0.d0)), dimension(nyg), save :: cflfacy
  real(KIND(0.d0)), dimension(nzg), save :: cflfacz

  ! for divb cleaning
  real(KIND(0.d0)), save :: dxyz_min

  integer, save :: mpi_comm_world_top    ! defined at remap_init.F90 (for top MPI thread)
  integer, save :: mpi_comm_world_bot    ! defined at remap_init.F90 (for bottom MPI thread)

end module geometry_def
