!========================================================================================|
!
!  R2D2 radiation MHD code
!  Copyright(C) Hideyuki Hotta <hotta.hideyuki@gmail.com>
!  
!========================================================================================|
subroutine remap_init
!======================================================================
#ifndef non_IO  
  use const_def, only: pi
  use geometry_def, only: ix0,jx0,kx0,nx,ny,nz,nxg,nyg,nzg,margin &
       & ,xmax,xmin,ymax,ymin,zmax,zmin,ix00,jx00,kx00,mtype &
       & ,zdcheck
  
  use implicit_def
  use comm_def, only: myrank,npe,ib,jb,kb
  use remap_def, only: &
       &  iix,jjx,kkx,iix0,jjx0,kkx0 &
       & ,i2ir,j2jr &
       & ,ir,jr,iixl,jjxl,iixlg,jjxlg &
       & ,np_ijr &
       & ,is,ie,js,je &
       & ,is0,ie0,js0,je0 &
       & ,remapcount_rcv,remapcount_snd &
       & ,remapijk_rcv,remapijk_snd &
       & ,qre &
       & ,ixr,jxr &
       & ,ssize ,gsize ,start &
       & ,ssize_in,gsize_in,start_in &
       & ,ssize_2d,gsize_2d,start_2d &
       & ,global_xy,global_xz,global_flux,global_top,global_2d &
       & ,m2d_xy,m2d_xz,m2d_flux,m_in,m_tu &
       & ,mpi_comm_world_remap &
       & ,mpi_comm_world_remap_xz &
       & ,mpi_comm_world_ir &
       & ,mtype_remap &
       & ,jc
  implicit none
  include "mpif.h"

  integer :: color,key
  integer :: np,np0
  integer :: il,jl
  integer :: ii0,ii1,iii0,iii1
  integer :: jj0,jj1,jjj0,jjj1
  integer :: kk0,kk1,kkk0,kkk1
  integer :: marginx0,marginx1
  integer :: marginy0,marginy1
  integer :: marginz0,marginz1
  integer :: iixlg_tmp,jjxlg_tmp

!======================================================================

  if(zdcheck == 2) then
  ixr = min(int(iix0/24),int(npe))   ! number of thread in x direction (remap)
  jxr = max(int(npe/ixr),1)  ! number of thread in y direction (remap)
  
#ifndef remap_2d_assign
  ixr = min(int(iix0),int(npe))
  jxr = 1
#endif
  
  allocate(np_ijr(ixr,jxr))

  do i = 1,margin
     i2ir(i) = 1
  enddo

  do i = iix-margin+1,iix
     i2ir(i) = ixr
  enddo

  do j = 1,margin
     j2jr(j) = 1
  enddo

  do j = jjx-margin+1,jjx
     j2jr(j) = jxr
  enddo
  
  i = margin
  j = margin
  do np = 0,npe-1
     jr(np) = np/ixr + 1 ! my x location (remap)
#ifndef remap_2d_assign
     jr(np) = 1
#endif     
     ir(np) = np - (jr(np)-1)*ixr + 1 ! my y location (remap)

     iixl(np) = iix0/ixr
     jjxl(np) = jjx0/jxr

     if(mod(iix0,ixr) /= 0) then
        if(ir(np) <= mod(iix0,ixr)) iixl(np) = iixl(np) + 1
     endif
     if(ir(np) >  ixr) then
        iixl(np) = 0
        jjxl(np) = 0
     endif

     !
     if(mod(jjx0,jxr) /= 0) then
        if(jr(np) <= mod(jjx0,jxr)) jjxl(np) = jjxl(np) + 1
     endif

     if(jr(np) >  jxr) then
        iixl(np) = 0
        jjxl(np) = 0
     endif

     if(jr(np) == 1 .and. iixl(np) /= 0) then
        do il = 1,iixl(np)
           i = i + 1
           i2ir(i) = ir(np)
        enddo
     endif
     
     if(ir(np) == 1 .and. jjxl(np) /= 0) then
        do jl = 1,jjxl(np)
           j = j + 1
           j2jr(j) = jr(np)
        enddo
     endif

     if(iixl(np) > 0 .and. jjxl(np) > 0) then
        np_ijr(ir(np),jr(np)) = np
     endif
  enddo

  is = 0
  is0 = 0
  js = 0
  js0 = 0
  do np = 0,npe-1
     if(iixl(np) /= 0 .and. jjxl(np) /= 0) then
        if(ir(np) == 1) then
           is(np)  = 1
           is0(np) = 1
        else
           np0 = np_ijr(ir(np)-1,jr(np))
           is(np)  = ie(np0) + 1             ! without margin
           is0(np) = ie0(np0) - 2*margin + 1 ! with margin
        endif
        ie(np)  = is(np)  + iixl(np) - 1
        ie0(np) = is0(np) + iixl(np) + 2*margin - 1

        if(jr(np) == 1) then
           js(np)  = 1
           js0(np) = 1
        else
           np0 = np_ijr(ir(np),jr(np)-1)
           js(np)  = je(np0) + 1              ! without margin
           js0(np) = je0(np0) - 2*margin + 1  ! with margin
        endif
        je(np)  = js(np)  + jjxl(np) - 1     
        je0(np) = js0(np) + jjxl(np) + 2*margin - 1
     endif
  enddo

  do np = 0,npe-1
     iixlg(np) = iixl(np) + 2*margin
     jjxlg(np) = jjxl(np) + 2*margin
  enddo

  if(iixl(myrank) /= 0 .and. jjxl(myrank) /= 0) then
     iixlg_tmp = iixlg(myrank)
     jjxlg_tmp = jjxlg(myrank)
     allocate(qre(iixlg_tmp,jjxlg_tmp,kkx,mtype_remap))
  endif
  
  remapcount_snd = 0 ! numbger of grids to send
  remapcount_rcv = 0

  !!!----------------------------!
  !!! x-range
  ii0 = ib*nx + 1  + margin
  ii1 = ib*nx + nx + margin
  marginx0 = margin
  marginx1 = margin
  if(ib == 0) then
     ii0 = ib*nx + 1
     marginx0 = 0
  endif

  if(ib == ix0 - 1) then
     ii1 = ib*nx + nx + 2*margin
     marginx1 = 0
  endif

  !!!----------------------------!
  !!! y-range
  jj0 = jb*ny + 1  + margin
  jj1 = jb*ny + ny + margin
  marginy0 = margin
  marginy1 = margin
  if(jb == 0) then
     jj0 = jb*ny + 1
     marginy0 = 0
  endif

  if(jb == jx0 - 1) then
     jj1 = jb*ny + ny + 2*margin
     marginy1 = 0
  endif

  !!!----------------------------!
  !!! z-range
  kk0 = kb*nz + 1  + margin
  kk1 = kb*nz + nz + margin
  marginz0 = margin
  marginz1 = margin
  if(kb == 0) then
     kk0 = kb*nz + 1
     marginz0 = 0
  endif

  if(kb == kx0 -1) then
     kk1 = kb*nz + nz + 2*margin
     marginz1 = 0
  endif

  ! count
!  do j = jj0,jj1
!  do i = ii0,ii1
!     np = np_ijr(i2ir(i),j2jr(j))
!     remapcount_snd(np) = remapcount_snd(np) + 1
!  enddo
!  enddo

  do np = 0,npe-1
!     if(remapcount_snd(np) /= 0) then
        iii0 = max(is0(np),ii0)
        iii1 = min(ie0(np),ii1)
        jjj0 = max(js0(np),jj0)
        jjj1 = min(je0(np),jj1)

        if(iii0 <= iii1 .and. jjj0 <= jjj1) then
           remapcount_snd(np) = 1
        endif

!        iii0 = max(is0(np) + marginx0,ii0)
!        iii1 = min(ie0(np) - marginx1,ii1)
!        jjj0 = max(js0(np) + marginy0,jj0)
!        jjj1 = min(je0(np) - marginy1,jj1)

        kkk0 = kk0
        kkk1 = kk1

        remapijk_snd(1,np) = iii0
        remapijk_snd(2,np) = jjj0
        remapijk_snd(3,np) = kkk0

        remapijk_snd(4,np) = iii1
        remapijk_snd(5,np) = jjj1
        remapijk_snd(6,np) = kkk1
!     endif
  enddo
  
  call mpi_alltoall(remapcount_snd(0:npe-1),1,mpi_integer &
       &           ,remapcount_rcv(0:npe-1),1,mpi_integer &
       &           ,mpi_comm_world,merr)

  call mpi_alltoall(remapijk_snd,6,mpi_integer &
       &           ,remapijk_rcv,6,mpi_integer &
       &           ,mpi_comm_world,merr)

  if(myrank == 0) then
     open(idf,file="data/remap/remap_info.dac",form="unformatted",access="stream")
     write(idf) is,ie,js,je,iixl,jjxl,np_ijr,ir,jr,i2ir,j2jr
     close(idf)
  endif

!======================================================================
! preparaing mpi io for remap

  if(iixl(myrank) /= 0 .and. jjxl(myrank) /= 0) then
     !--------------
     gsize(1) = iix0
     gsize(2) = jjx0
     gsize(3) = m2d_xy
     ssize(1) = iixl(myrank)
     ssize(2) = jjxl(myrank)
     ssize(3) = m2d_xy
     start(1) = is(myrank) - 1
     start(2) = js(myrank) - 1
     start(3) = 0     

     call mpi_type_create_subarray(3,gsize,ssize,start, &
          & mpi_order_fortran,mpi_real,global_xy,merr)  
     call mpi_type_commit(global_xy,merr)

     !-----------------
     gsize(1) = iix0
     gsize(2) = kkx0
     gsize(3) = m2d_xz
     ssize(1) = iixl(myrank)
     ssize(2) = kkx0
     ssize(3) = m2d_xz
     start(1) = is(myrank) - 1
     start(2) = 0
     start(3) = 0     
     
     call mpi_type_create_subarray(3,gsize,ssize,start, &
          & mpi_order_fortran,mpi_real,global_xz,merr)  
     call mpi_type_commit(global_xz,merr)
     
     !--------------
     gsize(1) = iix0 + 1
     gsize(2) = jjx0
     gsize(3) = m2d_flux
     ssize(1) = iixl(myrank) + 1
     ssize(2) = jjxl(myrank)
     ssize(3) = m2d_flux

     start(1) = is(myrank) - 1
     start(2) = js(myrank) - 1
     start(3) = 0

     call mpi_type_create_subarray(3,gsize,ssize,start, &
          & mpi_order_fortran,mpi_real,global_flux,merr)  
     call mpi_type_commit(global_flux,merr)
  endif
    
  !!!========================================
  !!! For intensity output
  !!! 
  if(ib == ix0 - 1) then
!!! define global data type !!!
     gsize_in(1) = m_tu
     gsize_in(2) = m_in
     gsize_in(3) = jjx0
     gsize_in(4) = kkx0
     ssize_in(1) = m_tu
     ssize_in(2) = m_in
     ssize_in(3) = ny
     ssize_in(4) = nz
     start_in(1) = 0
     start_in(2) = 0
     start_in(3) = jb*ny
     start_in(4) = kb*nz

     call mpi_type_create_subarray(4,gsize_in,ssize_in,start_in, &
          & mpi_order_fortran,mpi_real,global_top,merr)  
     call mpi_type_commit(global_top,merr)
  endif

  key = 0

  !===================
  color = 0
  if(iixl(myrank) == 0 .or. jjxl(myrank) == 0) color = 1
  ! split MPI_COMM_WORLD
  call mpi_comm_split(mpi_comm_world,color,key,mpi_comm_world_remap,merr)

  !===================
  color = 0
  if(iixl(myrank) == 0 .or. jjxl(myrank) == 0) color = 1
  if(js(myrank) > jc .or. je(myrank) < jc) color = 1
  ! split MPI_COMM_WORLD
  call mpi_comm_split(mpi_comm_world,color,key,mpi_comm_world_remap_xz,merr)

  endif ! zdcheck == 2

  !===================
  if(iixl(myrank) == 0 .or. jjxl(myrank) == 0) then
     color = 0
  else
     color = ir(myrank)
  endif
  call mpi_comm_split(mpi_comm_world,color,key,mpi_comm_world_ir,merr)

  if(zdcheck == 1) then
!!! define global data type for 2D data !!!
     gsize_2d(1) = mtype + 5
     gsize_2d(2) = iix0
     gsize_2d(3) = jjx0
     ssize_2d(1) = mtype + 5
     ssize_2d(2) = nx
     ssize_2d(3) = ny
     start_2d(1) = 0
     start_2d(2) = ib*nx
     start_2d(3) = jb*ny
  
     call mpi_type_create_subarray(3,gsize_2d,ssize_2d,start_2d, &
          & mpi_order_fortran,mpi_real,global_2d,merr)  
     call mpi_type_commit(global_2d,merr)
  endif
  
#endif
  return
end subroutine remap_init
