!========================================================================================|
!
!  R2D2 radiation MHD code
!  Copyright(C) Hideyuki Hotta <hotta.hieyuki@gmail.com>
!
!========================================================================================|
program main
!========================================================================================|

  use comm_def, only: myrank, npe
  use remap_def, only: jc, kc
  use geometry_def, only: npe0, mtype &
       & , nxg, nyg, nzg &
       & , nx, ny, nz  &
       & , ix0, jx0, kx0 &
       & , xdcheck, ydcheck, zdcheck, yinyang
  use mhd_def, only: qq
#ifndef non_IO
  use remap_def, only: iix0, jjx0, kkx0
#endif
  use time_def, only: cno, cnou,nd, nd_tau, ns, time, timep, mcont &
       & , dt
  use implicit_def
  implicit none
  include "mpif.h"

  ! directory category number
  integer :: ictg
  
  ! directory information
  logical exist

  ! for profiling
  real(KIND(0.d0)) :: time0, time1

  ! time control parameters
  real(KIND(0.d0)) :: dtout, dtout_tau, tend, ifac
  integer :: nstop

!----------------------------------------------------------------------------------------|
! initial MPI setting

  call mpi_init(merr)
  call mpi_comm_size(mpi_comm_world, npe, merr)
  call mpi_comm_rank(mpi_comm_world, myrank, merr)
  write (cno, '(i8.8)') myrank

!----------------------------------------------------------------------------------------|
! define parameter for reading

  mcont = 1 !read or not
  ! 0: without reading even if there is data
  ! 1: with reading if there is no data, no reading

!----------------------------------------------------------------------------------------|
! make data directory

#ifndef K
  if (myrank == 0) then ! K does not require this line.
#endif
    inquire (file='data/param/nd.dac', exist=exist)
    if (.not. exist) then
      call system("mkdir data")
      call system("mkdir data/time")
      call system("mkdir data/time/mhd")
      call system("mkdir data/time/tau")
      call system("mkdir data/qq")
      call system("mkdir data/param")
      call system("mkdir data/remap")
      call system("mkdir data/remap/qq")
      call system("mkdir data/remap/vl")
      call system("mkdir data/tau")
      call system("mkdir data/slice")
#ifdef particle
      call system("mkdir data/pos")
#endif
      
      mcont = 0
    end if
#ifndef K
  end if ! K does not require this line
  call mpi_bcast(mcont, 1, mpi_integer, 0, mpi_comm_world, merr) ! FX
#endif
  !  FX requires this line


  ictg = 1000
  write(cnou, '(i5.5)') myrank/ictg
#ifdef posixio_write
  !n_exist = access('data/qq/'//cno,'r')
  inquire (file='data/qq/'//cnou//'/'//cno//'/qq.dac.e.'//cno, exist=exist)
  
  if(.not. exist) then        
     call system("mkdir -p data/qq/"//cnou//'/'//cno)
  endif
#endif

  ! remap always use POSIX I/O
  inquire (file='data/remap/qq/'//cnou//'/'//cno//'/qq.dac.00000000.'//cno, exist=exist)
  if(.not. exist) then        
     call system("mkdir -p data/remap/qq/"//cnou//'/'//cno)
  endif  
  
!----------------------------------------------------------------------------------------|
! define block information

  call comm_init

  qq = 0.d0

!----------------------------------------------------------------------------------------|
! initialize remap

  call remap_init

!----------------------------------------------------------------------------------------|
  ! initialize slice

!----------------------------------------------------------------------------------------|

  call io_init

!----------------------------------------------------------------------------------------|
! read data (modification has not completed now)

  call io(dtout, dtout_tau, ifac, tend)

!========================================================================================|
!     epilogue
!========================================================================================|
  if (myrank == 0) write (*, *) "### NORMAL END ###"
! time develope end <---
!----------------------------------------------------------------------------------------|
  call mpi_finalize(merr)
  stop

end program main
