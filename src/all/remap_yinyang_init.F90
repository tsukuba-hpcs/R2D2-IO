!======================================================================
subroutine remap_yinyang_init
!======================================================================
#ifndef non_IO
#ifdef YinYang
  use const_def, only: pi
  use geometry_def, only: ix0, jx0, kx0, nx, ny, nz, nxg, nyg, nzg, margin &
       & , dx, dy, dz &
       & , xmax, xmin, ymax, ymin, zmax, zmin, ix00, jx00, kx00, mtype, y, z
  use yinyang_def, only: yo, zo, yzmask
  use implicit_def
  use remap_def, only: mtype_remap &
       & , iix, jjx, kkx, iix0, jjx0, kkx0 &
       & , xt, yt, zt, i2ir, j2jr &
       & , ir, jr &
       & , iixl, jjxl, iixlg, jjxlg &
       & , iixl_myrank, jjxl_myrank, iixlg_myrank, jjxlg_myrank &
       & , jkc, np_ijr &
       & , is, ie, js, je &
       & , is0, ie0, js0, je0 &
       & , remapcount_rcv, remapcount_snd &
       & , qre &
       & , icount, ixr, jxr &
       & , ssize_in, gsize_in, start_in &
       & , ssize, gsize, start &
       & , global_xy, global_xz, global_flux, global_spex, global_top &
       & , m_tu, m_in&
       & , m2d_xy, m2d_xz, m2d_spex, m2d_flux &
       & , mpi_comm_world_remap &
       & , mpi_comm_world_remap_xz &
       & , mpi_comm_world_spex, myrank_spex &
       & , mpi_comm_world_spex_io &
       & , mpi_comm_world_ir &
       & , mpi_comm_world_top_yy &
       & , jc

  use comm_def, only: npe, ib, jb, kb, myrank, mpi_yinyang

  implicit none
  include "mpif.h"

  integer :: color, key
  integer :: jj, kk
  real(KIND(0.d0)) :: xmax0 = xmax, xmin0 = xmin
  real(KIND(0.d0)) :: ymax0 = pi, ymin0 = 0.d0
  real(KIND(0.d0)) :: zmax0 = pi, zmin0 = -pi
  real(KIND(0.d0)) :: dxt, dyt, dzt
  integer :: ic
  integer :: jmin, jmax, kmin, kmax
  real(KIND(0.d0)) :: yo1, yo2, yo3, yo4
  real(KIND(0.d0)) :: zo1, zo2, zo3, zo4
  real(KIND(0.d0)) :: ymin00, ymax00, zmin00, zmax00
  real(KIND(0.d0)) :: yo0, zo0, sct, sco

  integer :: np, np0
  integer :: il, jl
  integer :: ii0, ii1

  integer :: ix, jxr_tmp, jx
  integer :: one = 1, two = 2, three = 3, four = 4
!======================================================================
  ! setting geometry in ordinary spherical geoemetry

  dxt = (xmax0 - xmin0)/dble(iix - 2*margin)
  dyt = (ymax0 - ymin0)/dble(jjx - 2*margin)
  dzt = (zmax0 - zmin0)/dble(kkx - 2*margin)

  xt(1) = 0.5d0*dxt + xmin0 - dble(margin)*dxt
  yt(1) = 0.5d0*dyt + ymin0 - dble(margin)*dyt
  zt(1) = 0.5d0*dzt + zmin0 - dble(margin)*dzt

  do i = 2, iix
    xt(i) = xt(i - 1) + dxt ! transformed geometry (ordinary spherical)
  enddo

  do j = 2, jjx
    yt(j) = yt(j - 1) + dyt ! transformed geometry (ordinary spherical)
  enddo

  do k = 2, kkx
    zt(k) = zt(k - 1) + dzt ! transformed geometry (ordinary spherical)
  enddo

  !! counting the number for sending
  if (mpi_yinyang == 0) then !!! YIN grid
    icount = 0
    !do kk = 1+margin,nzg-margin-1
    !do jj = 1+margin,nyg-margin-1
    do kk = 1 + margin - 1, nzg - margin
    do jj = 1 + margin - 1, nyg - margin
      jmin = int((y(jj) - yt(1))/dyt) + 1
      jmax = int((y(jj + 1) - yt(1))/dyt) + 2
      kmin = int((z(kk) - zt(1))/dzt) + 1
      kmax = int((z(kk + 1) - zt(1))/dzt) + 2

      do k = kmin, kmax
      do j = jmin, jmax
        icount = icount + 1
      enddo
      enddo
    enddo
    enddo

    allocate (jkc(four, icount))

    !! again counting the number for sending
    !! and storing
    ic = 0
    !do kk = 1+margin,nzg-margin-1
    !do jj = 1+margin,nyg-margin-1
    do kk = 1 + margin - 1, nzg - margin
    do jj = 1 + margin - 1, nyg - margin
      jmin = int((y(jj) - yt(1))/dyt) + 1
      jmax = int((y(jj + 1) - yt(1))/dyt) + 2
      kmin = int((z(kk) - zt(1))/dzt) + 1
      kmax = int((z(kk + 1) - zt(1))/dzt) + 2

      do k = kmin, kmax
      do j = jmin, jmax
        ic = ic + 1
        jkc(one, ic) = j
        jkc(two, ic) = k
        jkc(three, ic) = jj
        jkc(four, ic) = kk
      enddo
      enddo
    enddo
    enddo
  else !!!! For Yang grid
    icount = 0
    !do kk = 1+margin,nzg-margin-1
    !do jj = 1+margin,nyg-margin-1
    do kk = 1 + margin - 1, nzg - margin
    do jj = 1 + margin - 1, nyg - margin
      if (yzmask(jj, kk) == 1.d0) then
        yo1 = yo(jj, kk)
        yo2 = yo(jj, kk + 1)
        yo3 = yo(jj + 1, kk + 1)
        yo4 = yo(jj + 1, kk)

        zo1 = zo(jj, kk)
        zo2 = zo(jj, kk + 1)
        zo3 = zo(jj + 1, kk + 1)
        zo4 = zo(jj + 1, kk)

        ymin00 = min(yo1, yo2, yo3, yo4)
        ymax00 = max(yo1, yo2, yo3, yo4)

        zmin00 = min(zo1, zo2, zo3, zo4)
        zmax00 = max(zo1, zo2, zo3, zo4)

        jmin = int((ymin00 - yt(1))/dyt) + 1
        jmax = int((ymax00 - yt(1))/dyt) + 2
        kmin = int((zmin00 - zt(1))/dzt) + 1
        kmax = int((zmax00 - zt(1))/dzt) + 2

        if ((y(jj) - 0.5d0*pi)*(y(jj + 1) - 0.5d0*pi) < 0.d0) then
          ! Yang grid over the pole
          if ((z(kk) - 0.5d0*pi)*(z(kk + 1) - 0.5d0*pi) < 0.d0) then
            jmax = max(jmin, jmax)
            jmin = 1
          endif
          if ((z(kk) + 0.5d0*pi)*(z(kk + 1) + 0.5d0*pi) < 0.d0) then
            jmin = min(jmin, jmax)
            jmax = jjx
          endif

          do k = 1, kkx
          do j = jmin, jmax

            yo0 = acos(sin(yt(j))*sin(zt(k))) ! grid in Yang
            zo0 = asin(cos(yt(j))/sin(yo0))
            sct = sin(yt(j))*cos(zt(k))
            sco = -sin(yo0)*cos(zo0)
            if (sign(1.d0, sct) /= sign(1.d0, sco)) then
              zo0 = sign(1.d0, zo0)*pi - zo0 ! grid in Yang
            endif

            if (((y(jj) - yo0)*(y(jj + 1) - yo0) < 0.d0) .and. &
                 &   ((z(kk) - zo0)*(z(kk + 1) - zo0) < 0.d0) &
                 & ) then
              icount = icount + 1
            endif
          enddo
          enddo
        else
          do k = kmin, kmax
          do j = jmin, jmax

            yo0 = acos(sin(yt(j))*sin(zt(k))) ! grid in Yang
            zo0 = asin(cos(yt(j))/sin(yo0))
            sct = sin(yt(j))*cos(zt(k))
            sco = -sin(yo0)*cos(zo0)
            if (sign(1.d0, sct) /= sign(1.d0, sco)) then
              zo0 = sign(1.d0, zo0)*pi - zo0 ! grid in Yang
            endif

            if (((y(jj) - yo0)*(y(jj + 1) - yo0) < 0.d0) .and. &
                 &   ((z(kk) - zo0)*(z(kk + 1) - zo0) < 0.d0) &
                 & ) then
              icount = icount + 1
            endif
          enddo
          enddo
        endif
      endif ! yzmask
    enddo
    enddo

    if (allocated(jkc)) deallocate (jkc)
    allocate (jkc(four, icount))

    ic = 0
    !do kk = 1+margin,nzg-margin-1
    !do jj = 1+margin,nyg-margin-1
    do kk = 1 + margin - 1, nzg - margin
    do jj = 1 + margin - 1, nyg - margin
      if (yzmask(jj, kk) == 1.d0) then
        yo1 = yo(jj, kk)
        yo2 = yo(jj, kk + 1)
        yo3 = yo(jj + 1, kk + 1)
        yo4 = yo(jj + 1, kk)

        zo1 = zo(jj, kk)
        zo2 = zo(jj, kk + 1)
        zo3 = zo(jj + 1, kk + 1)
        zo4 = zo(jj + 1, kk)

        ymin00 = min(yo1, yo2, yo3, yo4)
        ymax00 = max(yo1, yo2, yo3, yo4)

        zmin00 = min(zo1, zo2, zo3, zo4)
        zmax00 = max(zo1, zo2, zo3, zo4)

        jmin = int((ymin00 - yt(1))/dyt) + 1
        jmax = int((ymax00 - yt(1))/dyt) + 2
        kmin = int((zmin00 - zt(1))/dzt) + 1
        kmax = int((zmax00 - zt(1))/dzt) + 2

        if ((y(jj) - 0.5d0*pi)*(y(jj + 1) - 0.5d0*pi) < 0.d0) then
          ! Yang grid over the pole
          if ((z(kk) - 0.5d0*pi)*(z(kk + 1) - 0.5d0*pi) < 0.d0) then
            jmax = max(jmin, jmax)
            jmin = 1
          endif
          if ((z(kk) + 0.5d0*pi)*(z(kk + 1) + 0.5d0*pi) < 0.d0) then
            jmin = min(jmin, jmax)
            jmax = jjx
          endif

          do k = 1, kkx
          do j = jmin, jmax
            yo0 = acos(sin(yt(j))*sin(zt(k))) ! grid in Yang
            zo0 = asin(cos(yt(j))/sin(yo0))
            sct = sin(yt(j))*cos(zt(k))
            sco = -sin(yo0)*cos(zo0)
            if (sign(1.d0, sct) /= sign(1.d0, sco)) then
              zo0 = sign(1.d0, zo0)*pi - zo0 ! grid in Yang
            endif

            if (((y(jj) - yo0)*(y(jj + 1) - yo0) < 0.d0) .and. &
                 &   ((z(kk) - zo0)*(z(kk + 1) - zo0) < 0.d0) &
                 & ) then
              ic = ic + 1
              jkc(one, ic) = j
              jkc(two, ic) = k
              jkc(three, ic) = jj
              jkc(four, ic) = kk
            endif
          enddo
          enddo
        else

          do k = kmin, kmax
          do j = jmin, jmax
            yo0 = acos(sin(yt(j))*sin(zt(k))) ! grid in Yang
            zo0 = asin(cos(yt(j))/sin(yo0))
            sct = sin(yt(j))*cos(zt(k))
            sco = -sin(yo0)*cos(zo0)
            if (sign(1.d0, sct) /= sign(1.d0, sco)) then
              zo0 = sign(1.d0, zo0)*pi - zo0 ! grid in Yang
            endif

            if (((y(jj) - yo0)*(y(jj + 1) - yo0) < 0.d0) .and. &
                 &   ((z(kk) - zo0)*(z(kk + 1) - zo0) < 0.d0) &
                 & ) then
              ic = ic + 1
              jkc(one, ic) = j
              jkc(two, ic) = k
              jkc(three, ic) = jj
              jkc(four, ic) = kk
            endif
          enddo
          enddo
        endif
      endif ! yzmask
    enddo
    enddo
  endif

  do ix = iix0/2, 1, -1
    if ((npe/ix)*ix == npe .and. (iix0/ix)*ix == iix0) then
      ixr = ix
      goto 1001
    endif
  enddo
1001 continue

  !ixr = min(int(iix0/8),int(npe))   ! number of thread in x direction (remap)
  jxr_tmp = max(int(npe/ixr), 1)  ! number of thread in y direction (remap)

  ! jjx must be divided by jxr for spex.f90 in values.f90
  !

  ! serch jxr which can divid jjx
  do jx = jxr_tmp, 1, -1
    if ((jjx0/jx)*jx == jjx0) then
      jxr = jx
      goto 1000
    endif
  enddo
1000 continue

#ifndef remap_2d_assign
  jxr = 1
#endif

  allocate (np_ijr(ixr, jxr))

  do i = 1, margin
    i2ir(i) = 1
  enddo

  do i = iix - margin + 1, iix
    i2ir(i) = ixr
  enddo

  do j = 1, margin
    j2jr(j) = 1
  enddo

  do j = jjx - margin + 1, jjx
    j2jr(j) = jxr
  enddo

  i = margin
  j = margin
  do np = 0, npe - 1
    jr(np) = np/ixr + 1 ! my x location (remap)
#ifndef remap_2d_assign
    jr(np) = 1
#endif
    ir(np) = np - (jr(np) - 1)*ixr + 1 ! my y location (remap)

    iixl(np) = iix0/ixr
    jjxl(np) = jjx0/jxr

    if (mod(iix0, ixr) /= 0) then
      if (ir(np) <= mod(iix0, ixr)) iixl(np) = iixl(np) + 1
    endif
    if (ir(np) > ixr) then
      iixl(np) = 0
      jjxl(np) = 0
    endif

    !
    if (mod(jjx0, jxr) /= 0) then
      if (jr(np) <= mod(jjx0, jxr)) jjxl(np) = jjxl(np) + 1
    endif
    if (jr(np) > jxr) then
      iixl(np) = 0
      jjxl(np) = 0
    endif

    if (jr(np) == 1 .and. iixl(np) /= 0) then
      do il = 1, iixl(np)
        i = i + 1
        i2ir(i) = ir(np)
      enddo
    endif

    if (ir(np) == 1 .and. jjxl(np) /= 0) then
      do jl = 1, jjxl(np)
        j = j + 1
        j2jr(j) = jr(np)
      enddo
    endif

    if (iixl(np) > 0 .and. jjxl(np) > 0) then
      np_ijr(ir(np), jr(np)) = np
    endif
  enddo

  is = 0
  is0 = 0
  js = 0
  js0 = 0
  do np = 0, npe - 1
    if (iixl(np) /= 0 .and. jjxl(np) /= 0) then
      if (ir(np) == 1) then
        is(np) = 1
        is0(np) = 1
      else
        np0 = np_ijr(ir(np) - 1, jr(np))
        is(np) = ie(np0) + 1
        is0(np) = ie0(np0) - 2*margin + 1
      endif
      ie(np) = is(np) + iixl(np) - 1
      ie0(np) = is0(np) + iixl(np) + 2*margin - 1

      if (jr(np) == 1) then
        js(np) = 1
        js0(np) = 1
      else
        np0 = np_ijr(ir(np), jr(np) - 1)
        js(np) = je(np0) + 1 ! without margin
        js0(np) = je0(np0) - 2*margin + 1 ! with margin
      endif
      je(np) = js(np) + jjxl(np) - 1
      je0(np) = js0(np) + jjxl(np) + 2*margin - 1
    endif
  enddo

  do np = 0, npe - 1
    iixlg(np) = iixl(np) + 2*margin
    jjxlg(np) = jjxl(np) + 2*margin
  enddo

  iixl_myrank = iixl(myrank)
  jjxl_myrank = jjxl(myrank)
  iixlg_myrank = iixlg(myrank)
  jjxlg_myrank = jjxlg(myrank)

  if (iixl(myrank) /= 0 .and. jjxl(myrank) /= 0) then
    allocate (qre(iixlg_myrank, jjxlg_myrank, kkx, mtype_remap))
  endif

  remapcount_snd = 0 ! numbger of grids to send
  remapcount_rcv = 0

  ii0 = ib*nx + 1 + margin
  ii1 = ib*nx + nx + margin
  if (ib == 0) then
    ii0 = ib*nx + 1
  endif

  if (ib == ix0 - 1) then
    ii1 = ib*nx + nx + 2*margin
  endif

  do ic = 1, icount
    j = jkc(one, ic)
    k = jkc(two, ic)
!     do ii = 1,nxg
!        i = ib*nx + ii
    do i = ii0, ii1
      np = np_ijr(i2ir(i), j2jr(j))
      remapcount_snd(np) = remapcount_snd(np) + 1
    enddo
  enddo

  call mpi_alltoall(remapcount_snd(0:npe - 1), 1, mpi_integer &
       &           , remapcount_rcv(0:npe - 1), 1, mpi_integer &
       &           , mpi_comm_world, merr)

  if (myrank == 0) then
    open (idf, file="data/remap/remap_info.dac", form="unformatted", access="stream")
    write (idf) is, ie, js, je, iixl, jjxl, np_ijr, ir, jr, i2ir, j2jr
    close (idf)
  endif

!======================================================================
! preparaing mpi io

  if (iixl(myrank) /= 0 .and. jjxl(myrank) /= 0) then
    !-----------------
    gsize(1) = iix0
    gsize(2) = jjx0
    gsize(3) = m2d_xy
    ssize(1) = iixl(myrank)
    ssize(2) = jjxl(myrank)
    ssize(3) = m2d_xy
    start(1) = is(myrank) - 1
    start(2) = js(myrank) - 1
    start(3) = 0

    call mpi_type_create_subarray(3, gsize, ssize, start, &
         & mpi_order_fortran, mpi_real, global_xy, merr)
    call mpi_type_commit(global_xy, merr)

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

    call mpi_type_create_subarray(3, gsize, ssize, start, &
         & mpi_order_fortran, mpi_real, global_xz, merr)
    call mpi_type_commit(global_xz, merr)

    !-----------------
    gsize(1) = iix0 + 1
    gsize(2) = jjx0
    gsize(3) = m2d_flux
    ssize(1) = iixl(myrank) + 1
    ssize(2) = jjxl(myrank)
    ssize(3) = m2d_flux

    start(1) = is(myrank) - 1
    start(2) = js(myrank) - 1
    start(3) = 0

    call mpi_type_create_subarray(3, gsize, ssize, start, &
         & mpi_order_fortran, mpi_real, global_flux, merr)
    call mpi_type_commit(global_flux, merr)

    !!!========================================
    !!! for spherical harmonic transformation
!!! define global data type !!!
    gsize(1) = iix0
    gsize(2) = kkx0/4
    gsize(3) = m2d_spex
    ssize(1) = iixl(myrank)
    ssize(2) = kkx0/4
    ssize(3) = m2d_spex
    start(1) = is(myrank) - 1
    start(2) = 0
    start(3) = 0

    call mpi_type_create_subarray(3, gsize, ssize, start, &
         & mpi_order_fortran, mpi_real, global_spex, merr)
    call mpi_type_commit(global_spex, merr)

  endif

  !!!========================================
  !!! For intensity output
  !!!
  if (ib == ix0 - 1) then
!!! define global data type !!!
    gsize_in(1) = m_tu
    gsize_in(2) = m_in
    gsize_in(3) = jx0*ny
    gsize_in(4) = kx0*nz
    ssize_in(1) = m_tu
    ssize_in(2) = m_in
    ssize_in(3) = ny
    ssize_in(4) = nz
    start_in(1) = 0
    start_in(2) = 0
    start_in(3) = jb*ny
    start_in(4) = kb*nz

    call mpi_type_create_subarray(4, gsize_in, ssize_in, start_in, &
         & mpi_order_fortran, mpi_real, global_top, merr)
    call mpi_type_commit(global_top, merr)
  endif

  key = 0

  !===================
  color = 0
  if (iixl(myrank) == 0 .or. jjxl(myrank) == 0) color = 1
  ! split MPI_COMM_WORLD
  call mpi_comm_split(mpi_comm_world, color, key, mpi_comm_world_remap, merr)


  !===================
  color = 0
  if(iixl(myrank) == 0 .or. jjxl(myrank) == 0) color = 1
  if(js(myrank) > jc .or. je(myrank) < jc) color = 1
  ! split MPI_COMM_WORLD
  call mpi_comm_split(mpi_comm_world,color,key,mpi_comm_world_remap_xz,merr)

  !===================
  color = ir(myrank)
  if (iixl(myrank) == 0 .or. jjxl(myrank) == 0) color = 0
  call mpi_comm_split(mpi_comm_world, color, key, mpi_comm_world_spex, merr)
  call mpi_comm_rank(mpi_comm_world_spex, myrank_spex, merr)

  !===================
  color = ir(myrank)
  if (iixl(myrank) == 0 .or. jjxl(myrank) == 0) color = 0
  call mpi_comm_split(mpi_comm_world, color, key, mpi_comm_world_ir, merr)

  !===================
  color = 0
  if (jr(myrank) /= 1 .or. jjxl(myrank) == 0) color = 1
  ! split MPI_COMM_WORLD
  call mpi_comm_split(mpi_comm_world, color, key, mpi_comm_world_spex_io, merr)

  !===================
  if (ib == ix0 - 1) then
    if (mpi_yinyang == 0) then
      color = 1
    else
      color = 2
    endif
  else
    color = 0
  endif
  ! split MPI_COMM_WORLD
  call mpi_comm_split(mpi_comm_world, color, key, mpi_comm_world_top_yy, merr)

  ! Yin-Yang
#endif
  ! non_IO
#endif
  return
end subroutine remap_yinyang_init
