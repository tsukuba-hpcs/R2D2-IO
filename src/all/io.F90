!----------------------------------------------------------------------------------------|
subroutine io(dtout, dtout_tau, ifac, tend)
!----------------------------------------------------------------------------------------|
! This subrotine is for data I/O
!----------------------------------------------------------------------------------------|
#ifndef non_IO
#ifdef ifort
  use iflport, only: hostnm
#endif

#ifdef fujitsu
  use service_routines, only: hostnm
#endif

  use const_def, only: sb, pi
#ifdef ideal
  use const_def, only: gm, theta, m_ad, k_cond, sigma
#endif
  use geometry_def, only: &
       &  xdcheck, ydcheck, zdcheck &
       & , marginx, marginy, marginz &
       & , margin, nx, ny, nz, nxg, nyg, nzg &
       & , ix0, jx0, kx0, mtype &
       & , xmax, ymax, zmax, xmin, ymin, zmin &
       & , i1 &
       & , ununiform_flag, ix_ununi, dx00

  use implicit_def

  use io_def, only: io_params_write, fh_io, global_io, mmx_io &
       & , fh_io_memory, global_io_memory, mmx_io_memory &
       & , memory_saving
  use remap_def, only: jc, kc
  use remap_def, only: ixr, jxr, m_in, m_tu, disp

  use comm_def, only: npe, myrank, ib, mpi_yinyang
  use time_def, only: cj, cnd, cnd_tau, nd, nd_tau, nd0, ns, time, timep, mcont, time00,cno,cnou

  use mhd_def, only: qq

#ifdef YinYang  
  use yinyang_def, only: mpi_comm_world_yy
#endif  

  implicit none

  include "mpif.h"

  real(KIND(0.d0)), dimension(nxg, nyg, nzg) :: qqs
  real(KIND(0.d0)) :: dtout, dtout_tau, tend, ifac

  ! for server name detect
  integer server_status
  character(len=20) :: server_name

  ! if true we read nd data, else 'e' or 'o' data.
  logical :: read_nd_flag = .false.
  character(len=8) :: cstep ! cnd or cj
  character(len=2) :: cm

  ! UUID for input_data
  character(len=36) :: uuid_stratification, uuid_stratification_cont
  character(len=36) :: uuid_eos

  logical :: mpi_stop_flag = .false.
  
  ! directory information
  logical :: exist

  integer :: swap
  integer :: deep_flag

  logical :: io_rte_fr_flag, io_rte_in_flag
  real(KIND=KIND(0.d0)), dimension(nxg, nyg, nzg) :: qrad
  integer, dimension(mpi_status_size) :: mstatus
  real(KIND=KIND(0.d0)), dimension(nxg,nyg,nzg) :: tu,fr_io

!----------------------------------------------------------------------------------------|
!----------------------------------------------------------------------------------------|

  write(cnd, '(i8.8)') 0
  cj = 'e'
  if (mcont == 1) then
!----------------------------------------------------------------------------------------|
#ifdef posixio_read

    ! 1. output variables
  open(idf,file='data/qq/'//cnou//'/'//cno//'/qq.dac.'//trim(cstep)//'.'//cno,form='unformatted',access="stream")
  read(idf) qq
  close(idf)
  
#else
! posixio_read

    if (memory_saving .and. nd == 0) then
      do m = 1, mtype
        write (cm, '(i2.2)') m - 1
        !call mpi_file_open(mpi_comm_world,'data/qq/qq.dac.'//cnd, &
#ifdef YinYang
        if (mpi_yinyang == 0) then
          call mpi_file_open(mpi_comm_world_yy, 'data/qq/qq_yin'//cm//'.dac.'//trim(cstep), &
               & mpi_mode_rdonly, mpi_info_null, fh_io_memory, merr)
        else
          call mpi_file_open(mpi_comm_world_yy, 'data/qq/qq_yan'//cm//'.dac.'//trim(cstep), &
               & mpi_mode_rdonly, mpi_info_null, fh_io_memory, merr)
        endif
#else
        call mpi_file_open(mpi_comm_world, 'data/qq/qq'//cm//'.dac.'//trim(cstep), &
             & mpi_mode_rdonly, mpi_info_null, fh_io_memory, merr)
#endif
        call mpi_file_set_view(fh_io_memory, disp, mpi_real8, global_io_memory, "native", mpi_info_null, merr)
#ifdef no_collective_mpiio
        call mpi_file_read(fh_io_memory, qqs, mmx_io_memory, mpi_real8, mstatus, merr)
#else
        call mpi_file_read_all(fh_io_memory, qqs, mmx_io_memory, mpi_real8, mstatus, merr)
#endif
        call mpi_file_close(fh_io_memory, merr)

        do k = 1, nzg
        do j = 1, nyg
        do i = 1, nxg
          qq(i, j, k, m) = qqs(i, j, k)
        enddo
        enddo
        enddo
      enddo
    else
#ifdef YinYang
      if (mpi_yinyang == 0) then
        call mpi_file_open(mpi_comm_world_yy, 'data/qq/qq_yin.dac.'//trim(cstep), &
             & mpi_mode_rdonly, mpi_info_null, fh_io, merr)
      else
        call mpi_file_open(mpi_comm_world_yy, 'data/qq/qq_yan.dac.'//trim(cstep), &
             & mpi_mode_rdonly, mpi_info_null, fh_io, merr)
      endif
#else
      call mpi_file_open(mpi_comm_world, 'data/qq/qq.dac.'//trim(cstep), &
           & mpi_mode_rdonly, mpi_info_null, fh_io, merr)
      !call mpi_file_open(mpi_comm_world,'data/qq/qq.dac.'//cnd, &
      !     & mpi_mode_rdonly,mpi_info_null,fh_io,merr)
#endif

      call mpi_file_set_view(fh_io, disp, mpi_real8, global_io, "native", mpi_info_null, merr)
#ifdef no_collective_mpiio
      call mpi_file_read(fh_io, qq, mmx_io, mpi_real8, mstatus, merr)
#else
      call mpi_file_read_all(fh_io, qq, mmx_io, mpi_real8, mstatus, merr)
#endif
      call mpi_file_close(fh_io, merr)
    endif

#endif
! posixio_read

!----------------------------------------------------------------------------------------|
!----------------------------------------------------------------------------------------|
!----------------------------------------------------------------------------------------|

  endif ! without reading (not continuous calculation)
!----------------------------------------------------------------------------------------|
!----------------------------------------------------------------------------------------|
!----------------------------------------------------------------------------------------|
  

  if ((int(timep/dtout) < int(time/dtout)) .or. ns == 0) then

    !----------------------------------------------------------------------------------------|

#ifdef posixio_write
  !-----
  ! 1. output variables
  open(idf,file='data/qq/'//cnou//'/'//cno//'/qq.dac.'//trim(cj)//'.'//cno,form='unformatted',access="stream")
  write(idf) qq
  close(idf)

  if(mod(int(nd),10) == 0) then
    open(idf,file='data/qq/'//cnou//'/'//cno//'/qq.dac.'//cnd//'.'//cno,form='unformatted',access="stream")
    write(idf) qq
    close(idf)
 endif

#else
!posixio_write
  
#ifdef YinYang
    if (mpi_yinyang == 0) then
      call mpi_file_open(mpi_comm_world_yy, 'data/qq/qq_yin.dac.'//trim(cj), &
           & ior(mpi_mode_create, mpi_mode_wronly), mpi_info_null, fh_io, merr)
    else
      call mpi_file_open(mpi_comm_world_yy, 'data/qq/qq_yan.dac.'//trim(cj), &
           & ior(mpi_mode_create, mpi_mode_wronly), mpi_info_null, fh_io, merr)
    endif
#else
    call mpi_file_open(mpi_comm_world, 'data/qq/qq.dac.'//trim(cj), &
         & ior(mpi_mode_create, mpi_mode_wronly), mpi_info_null, fh_io, merr)
#endif
    call mpi_file_set_view(fh_io, disp, mpi_real8, global_io, "native", mpi_info_null, merr)
#ifdef no_collective_mpiio
    call mpi_file_write(fh_io, qq, mmx_io, mpi_real8, mstatus, merr)
#else
    call mpi_file_write_all(fh_io, qq, mmx_io, mpi_real8, mstatus, merr)
#endif
    call mpi_file_close(fh_io, merr)

    if (mod(int(nd), 10) == 0) then
#ifdef YinYang
      if (mpi_yinyang == 0) then
        call mpi_file_open(mpi_comm_world_yy, 'data/qq/qq_yin.dac.'//cnd, &
             & ior(mpi_mode_create, mpi_mode_wronly), mpi_info_null, fh_io, merr)
      else
        call mpi_file_open(mpi_comm_world_yy, 'data/qq/qq_yan.dac.'//cnd, &
             & ior(mpi_mode_create, mpi_mode_wronly), mpi_info_null, fh_io, merr)
      endif
#else
      call mpi_file_open(mpi_comm_world, 'data/qq/qq.dac.'//cnd, &
           & ior(mpi_mode_create, mpi_mode_wronly), mpi_info_null, fh_io, merr)
#endif
      call mpi_file_set_view(fh_io, disp, mpi_real8, global_io, "native", mpi_info_null, merr)
#ifdef no_collective_mpiio
      call mpi_file_write(fh_io, qq, mmx_io, mpi_real8, mstatus, merr)
#else
      call mpi_file_write_all(fh_io, qq, mmx_io, mpi_real8, mstatus, merr)
#endif
      call mpi_file_close(fh_io, merr)
    endif

#endif
!posixio_write
    !----------------------------------------------------------------------------------------|
    call remap(tu, fr_io)
 endif

1000 continue
  !----------------------------------------------------------------------------------------|

  return
913 format(1x, '### write    ', 'step=', i8, ' time=', e10.3, ' nd =', i4, ' ###')
! non_IO
#endif
end subroutine io

