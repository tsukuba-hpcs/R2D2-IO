!========================================================================================|
!
!  R2D2 radiation MHD code
!  Copyright(C) 2019 Hideyuki Hotta <hotta@chiba-u.jp>
!  
!========================================================================================|
subroutine io_init
  use geometry_def, only: mtype,nxg,nyg,nzg,marginx,marginy,marginz,nx,ny,nz,ix0,jx0,kx0
  use comm_def, only: ib,jb,kb
  use io_def, only: gsize_io,ssize_io,start_io,global_io,mmx_io &
       & ,gsize_io_memory,ssize_io_memory,start_io_memory,global_io_memory,mmx_io_memory &
       & ,mpi_comm_world_jk  
  use implicit_def

  implicit none
  integer :: color,key
  
  include "mpif.h"
  

  gsize_io(1) = nx*ix0 + 2*marginx
  gsize_io(2) = ny*jx0 + 2*marginy
  gsize_io(3) = nz*kx0 + 2*marginz
  gsize_io(4) = mtype

  ssize_io(1) = nxg
  ssize_io(2) = nyg
  ssize_io(3) = nzg
  ssize_io(4) = mtype

  start_io(1) = nx*ib
  start_io(2) = ny*jb
  start_io(3) = nz*kb
  start_io(4) = 0
 
  !mmx_io = mtype*nx*ny*nz
  mmx_io = mtype*nxg*nyg*nzg

  call mpi_type_create_subarray(4,gsize_io,ssize_io,start_io, &
       & mpi_order_fortran,mpi_real8,global_io,merr)  
  call mpi_type_commit(global_io,merr)

  !======================================================
  !======================================================
  !======================================================
  
  gsize_io_memory(1) = nx*ix0 + 2*marginx
  gsize_io_memory(2) = ny*jx0 + 2*marginy
  gsize_io_memory(3) = nz*kx0 + 2*marginz

  ssize_io_memory(1) = nxg
  ssize_io_memory(2) = nyg
  ssize_io_memory(3) = nzg

  start_io_memory(1) = nx*ib
  start_io_memory(2) = ny*jb
  start_io_memory(3) = nz*kb
 
  mmx_io_memory = nxg*nyg*nzg

  call mpi_type_create_subarray(3,gsize_io_memory,ssize_io_memory,start_io_memory, &
       & mpi_order_fortran,mpi_real8,global_io_memory,merr)  
  call mpi_type_commit(global_io_memory,merr)

  !======================================================
  !======================================================
  !======================================================

  color = jb + jx0*kb
  key = ib

  ! split MPI_COMM_WORLD
  ! MPI_COMM_WORLD_jk is the world for the same horizontal position
  ! used for the RTE output
  call mpi_comm_split(mpi_comm_world,color,key,mpi_comm_world_jk,merr)
  call mpi_comm_rank(mpi_comm_world_jk, i, merr)

  return
end subroutine io_init
