!======================================================================
subroutine remap_yinyang(tu, fr)
!======================================================================
#ifndef non_IO
#ifdef YinYang
  use const_def, only: pi
  use geometry_def, only: ix0, jx0, kx0, nx, ny, nz, nxg, nyg, nzg, margin &
       & , xmax, xmin, ymax, ymin, zmax, zmin, ix00, jx00, kx00 &
       & , dx, dy, dz &
       & , y, z, mtype
  use time_def, only: cnd, cno
  use implicit_def
  use remap_def, only: &
       & iix, jjx, kkx, iix0, jjx0, kkx0 &
       & , yt, zt, i2ir, j2jr &
       & , ir, jr &
       & , jkc, np_ijr &
       & , iixl_myrank, jjxl_myrank &
       & , iixlg_myrank, jjxlg_myrank &
       & , is0, ie0, js0, je0 &
       & , remapcount_rcv, remapcount_snd &
       & , qre &
       & , icount, ixr, jxr, mtype_remap

  use comm_def, only: ib, myrank, npe, mpi_yinyang
  use mhd_def, only: qq

!$ use omp_lib
  implicit none
  include "mpif.h"
  integer, dimension(mpi_status_size) :: mstatus

  real(KIND(0.d0)) :: yo0, zo0
  integer :: ii, jj, kk
  integer :: kp, np
  integer :: ic

  real(KIND(0.d0)), dimension(nxg, nyg, nzg, mtype_remap) :: qq_remap
  real(KIND(0.d0)), dimension(nxg, nyg, nzg) :: tu, fr
  real(KIND(0.e0)) :: sins, coss
  real(KIND(0.d0)) :: sco, sct

  real(KIND(0.e0)), dimension(6) :: qzzz
  real(KIND(0.e0)), dimension(mtype_remap) :: qqq

  real(KIND(0.e0)), dimension(:, :), allocatable :: bufsnd_qq, bufrcv_qq
  integer, dimension(:, :), allocatable :: bufsnd_ijk, bufrcv_ijk
  real(KIND(0.e0)), dimension(mtype_remap, margin, jjxlg_myrank, kkx) :: bufsndxd, bufsndxu
  real(KIND(0.e0)), dimension(mtype_remap, iixlg_myrank, margin, kkx) :: bufsndyd, bufsndyu
  real(KIND(0.e0)), dimension(mtype_remap, margin, jjxlg_myrank, kkx) :: bufrcvxd, bufrcvxu
  real(KIND(0.e0)), dimension(mtype_remap, iixlg_myrank, margin, kkx) :: bufrcvyd, bufrcvyu
  integer, dimension(:, :), allocatable :: iiic
  integer :: mmx
  integer :: ireq_snd_qq, ireq_rcv_qq
  integer :: ireq_snd_ijk, ireq_rcv_ijk
  integer :: ireq_sndxd, ireq_sndxu
  integer :: ireq_sndyd, ireq_sndyu
  integer :: ireq_rcvxd, ireq_rcvxu
  integer :: ireq_rcvyd, ireq_rcvyu

  integer :: remapcount_snd0, remapcount_rcv0
  integer :: jst, jet, kst, ket
  integer :: il, jl
  integer :: ii0, ii1
  integer :: jj0, jj1
  integer :: iii0, iii1
  integer :: margin0, margin1

  real(KIND(0.e0)) :: ym1, yn0, yp1, yp2
  real(KIND(0.e0)) :: zm1, zn0, zp1, zp2

  integer :: remapcount_snd_np, remapcount_rcv_np
  integer :: one = 1, two = 2, three = 3, four = 4
!======================================================================

  do k = 1, nzg
  do j = 1, nyg
  do i = 1, nxg
    qq_remap(i, j, k, 1:mtype) = qq(i, j, k, 1:mtype)
    qq_remap(i, j, k, mtype + 1) = fr(i, j, k)
    qq_remap(i, j, k, mtype + 2) = tu(i, j, k)
  enddo
  enddo
  enddo
  
  ii0 = ib*nx + 1 + margin
  ii1 = ib*nx + nx + margin
  margin0 = margin
  margin1 = margin
  if (ib == 0) then
    ii0 = ib*nx + 1
    margin0 = 0
  endif

  if (ib == ix0 - 1) then
    ii1 = ib*nx + nx + 2*margin
    margin1 = 0
  endif

  jj0 = minval(jkc(one, :))
  jj1 = maxval(jkc(one, :))

  do np = npe - 1, 0, -1
    ! counting and constructind snd data
    remapcount_snd_np = remapcount_snd(np)
    remapcount_rcv_np = remapcount_rcv(np)
    if (remapcount_snd_np /= 0) then
      !! allocate the sending data
      ! location in sent geometry
      ! 1: i, 2: j, 3: k
      if (allocated(bufsnd_ijk)) deallocate (bufsnd_ijk)
      allocate (bufsnd_ijk(three, remapcount_snd_np))
      ! transformed data
      if (allocated(bufsnd_qq)) deallocate (bufsnd_qq)
      allocate (bufsnd_qq(remapcount_snd_np, mtype_remap))
      !! 1: ii, 2: ic
      if (allocated(iiic)) deallocate (iiic)
      allocate (iiic(two, remapcount_snd_np))

      if ((ie0(np) + margin0 >= ii0) .and. (is0(np) - margin1 <= ii1) .and. &
           & (je0(np) + margin >= jj0) .and. (js0(np) - margin <= jj1)) then
        iii0 = max(is0(np) + margin0, ii0)
        iii1 = min(ie0(np) - margin1, ii1)

        ! serching sending directory
        remapcount_snd0 = 0
        do ic = 1, icount
          j = jkc(one, ic)
          k = jkc(two, ic)
          do i = iii0, iii1
            if (np == np_ijr(i2ir(i), j2jr(j))) then
              il = i - is0(np) + 1
              jl = j - js0(np) + 1
              remapcount_snd0 = remapcount_snd0 + 1
              bufsnd_ijk(one, remapcount_snd0) = il
              bufsnd_ijk(two, remapcount_snd0) = jl
              bufsnd_ijk(three, remapcount_snd0) = k

              iiic(one, remapcount_snd0) = i - ib*nx
              iiic(two, remapcount_snd0) = ic
            endif
          enddo ! ii = 1,nxg
        enddo ! ic = 1,icount
      endif ! is ...

      mmx = 3*remapcount_snd(np)
      call mpi_isend(bufsnd_ijk, mmx, mpi_integer, np, 10 &
           & , mpi_comm_world, ireq_snd_ijk, merr)

      if (mpi_yinyang == 0) then
        !$omp parallel do private(ic,j,k,ii,jj,kk,remapcount_snd0,m,n &
        !$omp & ,kp &
        !$omp & ,qzzz,ym1,yn0,yp1,yp2,zm1,zn0,zp1,zp2)
        do remapcount_snd0 = 1, remapcount_snd_np
          ic = iiic(two, remapcount_snd0)
          j = jkc(one, ic)
          k = jkc(two, ic)

          ii = iiic(one, remapcount_snd0)
          jj = jkc(three, ic)
          kk = jkc(four, ic)

          ! Lagrange interpolation 3rd order
          ! for j - 1
          ym1 = real(&
               &  (yt(j) - y(jj))*(yt(j) - y(jj + 1))*(yt(j) - y(jj + 2))/ &
               & ((y(jj - 1) - y(jj))*(y(jj - 1) - y(jj + 1))*(y(jj - 1) - y(jj + 2))))

          ! for j
          yn0 = real(&
               &  (yt(j) - y(jj - 1))*(yt(j) - y(jj + 1))*(yt(j) - y(jj + 2))/ &
               & ((y(jj) - y(jj - 1))*(y(jj) - y(jj + 1))*(y(jj) - y(jj + 2))))

          ! for j + 1
          yp1 = real(&
               &  (yt(j) - y(jj - 1))*(yt(j) - y(jj))*(yt(j) - y(jj + 2))/ &
               & ((y(jj + 1) - y(jj - 1))*(y(jj + 1) - y(jj))*(y(jj + 1) - y(jj + 2))))

          ! for j + 2
          yp2 = real(&
               &  (yt(j) - y(jj - 1))*(yt(j) - y(jj))*(yt(j) - y(jj + 1))/ &
               & ((y(jj + 2) - y(jj - 1))*(y(jj + 2) - y(jj))*(y(jj + 2) - y(jj + 1))))

          ! for k - 1
          zm1 = real(&
               &  (zt(k) - z(kk))*(zt(k) - z(kk + 1))*(zt(k) - z(kk + 2))/ &
               & ((z(kk - 1) - z(kk))*(z(kk - 1) - z(kk + 1))*(z(kk - 1) - z(kk + 2))))

          ! for k
          zn0 = real(&
               &  (zt(k) - z(kk - 1))*(zt(k) - z(kk + 1))*(zt(k) - z(kk + 2))/ &
               & ((z(kk) - z(kk - 1))*(z(kk) - z(kk + 1))*(z(kk) - z(kk + 2))))

          ! for k + 1
          zp1 = real(&
               &  (zt(k) - z(kk - 1))*(zt(k) - z(kk))*(zt(k) - z(kk + 2))/ &
               & ((z(kk + 1) - z(kk - 1))*(z(kk + 1) - z(kk))*(z(kk + 1) - z(kk + 2))))

          ! for k + 2
          zp2 = real(&
               &  (zt(k) - z(kk - 1))*(zt(k) - z(kk))*(zt(k) - z(kk + 1))/ &
               & ((z(kk + 2) - z(kk - 1))*(z(kk + 2) - z(kk))*(z(kk + 2) - z(kk + 1))))

          do m = 1, mtype_remap
            ! 3rd order Lagrange interpolation
            do n = 1, 4
              kp = n - 2
              qzzz(n) = real( &
                   & qq_remap(ii, jj - 1, kk + kp, m)*ym1 + &
                   & qq_remap(ii, jj, kk + kp, m)*yn0 + &
                   & qq_remap(ii, jj + 1, kk + kp, m)*yp1 + &
                   & qq_remap(ii, jj + 2, kk + kp, m)*yp2)
            enddo

            bufsnd_qq(remapcount_snd0, m) = &
                 & qzzz(1)*zm1 + &
                 & qzzz(2)*zn0 + &
                 & qzzz(3)*zp1 + &
                 & qzzz(4)*zp2

          enddo
        enddo
      else
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! FOR YANG GRID
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !$omp parallel do private(remapcount_snd0,ic,j,k &
        !$omp & ,yo0,zo0,sct,sco,coss,sins &
        !$omp & ,ii,jj,kk &
        !$omp & ,ym1,yn0,yp1,yp2,zm1,zn0,zp1,zp2 &
        !$omp & ,kp,m,n &
        !$omp & ,qzzz,qqq)
        do remapcount_snd0 = 1, remapcount_snd_np
          ic = iiic(two, remapcount_snd0)
          j = jkc(one, ic)
          k = jkc(two, ic)
          yo0 = acos(sin(yt(j))*sin(zt(k))) ! grid in Yang
          zo0 = asin(cos(yt(j))/sin(yo0))
          sct = sin(yt(j))*cos(zt(k))
          sco = -sin(yo0)*cos(zo0)
          if (sign(1.d0, sct) /= sign(1.d0, sco)) then
            zo0 = sign(1.d0, zo0)*pi - zo0 ! grid in Yang
          endif
          coss = real(-sin(zo0)*sin(zt(k)))
          sins = real(cos(zo0)/sin(yt(j)))

          ii = iiic(one, remapcount_snd0)
          jj = jkc(three, ic)
          kk = jkc(four, ic)

          ! Lagrange interpolation 3rd order
          ! for j - 1
          ym1 = real( &
               &  (yo0 - y(jj))*(yo0 - y(jj + 1))*(yo0 - y(jj + 2))/ &
               & ((y(jj - 1) - y(jj))*(y(jj - 1) - y(jj + 1))*(y(jj - 1) - y(jj + 2))))

          ! for j
          yn0 = real(&
               &  (yo0 - y(jj - 1))*(yo0 - y(jj + 1))*(yo0 - y(jj + 2))/ &
               & ((y(jj) - y(jj - 1))*(y(jj) - y(jj + 1))*(y(jj) - y(jj + 2))))

          ! for j + 1
          yp1 = real(&
               &  (yo0 - y(jj - 1))*(yo0 - y(jj))*(yo0 - y(jj + 2))/ &
               & ((y(jj + 1) - y(jj - 1))*(y(jj + 1) - y(jj))*(y(jj + 1) - y(jj + 2))))

          ! for j + 2
          yp2 = real(&
               &  (yo0 - y(jj - 1))*(yo0 - y(jj))*(yo0 - y(jj + 1))/ &
               & ((y(jj + 2) - y(jj - 1))*(y(jj + 2) - y(jj))*(y(jj + 2) - y(jj + 1))))

          ! for k - 1
          zm1 = real(&
               &  (zo0 - z(kk))*(zo0 - z(kk + 1))*(zo0 - z(kk + 2))/ &
               & ((z(kk - 1) - z(kk))*(z(kk - 1) - z(kk + 1))*(z(kk - 1) - z(kk + 2))))

          ! for k
          zn0 = real(&
               &  (zo0 - z(kk - 1))*(zo0 - z(kk + 1))*(zo0 - z(kk + 2))/ &
               & ((z(kk) - z(kk - 1))*(z(kk) - z(kk + 1))*(z(kk) - z(kk + 2))))

          ! for k + 1
          zp1 = real(&
               &  (zo0 - z(kk - 1))*(zo0 - z(kk))*(zo0 - z(kk + 2))/ &
               & ((z(kk + 1) - z(kk - 1))*(z(kk + 1) - z(kk))*(z(kk + 1) - z(kk + 2))))

          ! for k + 2
          zp2 = real(&
               &  (zo0 - z(kk - 1))*(zo0 - z(kk))*(zo0 - z(kk + 1))/ &
               & ((z(kk + 2) - z(kk - 1))*(z(kk + 2) - z(kk))*(z(kk + 2) - z(kk + 1))))

          do m = 1, mtype_remap

            do n = 1, 4
              kp = n - 2
              qzzz(n) = real(&
                   & qq_remap(ii, jj - 1, kk + kp, m)*ym1 + &
                   & qq_remap(ii, jj, kk + kp, m)*yn0 + &
                   & qq_remap(ii, jj + 1, kk + kp, m)*yp1 + &
                   & qq_remap(ii, jj + 2, kk + kp, m)*yp2)
            enddo
            qqq(m) = &
                 & qzzz(1)*zm1 + &
                 & qzzz(2)*zn0 + &
                 & qzzz(3)*zp1 + &
                 & qzzz(4)*zp2
          enddo

          bufsnd_qq(remapcount_snd0, 1) = qqq(1)
          bufsnd_qq(remapcount_snd0, 2) = qqq(2)
          bufsnd_qq(remapcount_snd0, 5) = qqq(5)
          bufsnd_qq(remapcount_snd0, 8) = qqq(8)
          bufsnd_qq(remapcount_snd0, 9) = qqq(9)

          bufsnd_qq(remapcount_snd0, 3) = qqq(3)*coss - qqq(4)*sins
          bufsnd_qq(remapcount_snd0, 4) = qqq(3)*sins + qqq(4)*coss
          bufsnd_qq(remapcount_snd0, 6) = qqq(6)*coss - qqq(7)*sins
          bufsnd_qq(remapcount_snd0, 7) = qqq(6)*sins + qqq(7)*coss

          ! deal with other variable
          ! especially flux
          do m = mtype + 1, mtype_remap
            bufsnd_qq(remapcount_snd0, m) = qqq(m)
          enddo
        enddo
      endif

      mmx = mtype_remap*remapcount_snd_np
      call mpi_isend(bufsnd_qq, mmx, mpi_real, np, 20 &
           & , mpi_comm_world, ireq_snd_qq, merr)
    endif

    ! recieving
    if (remapcount_rcv_np /= 0) then
      if (allocated(bufrcv_ijk)) deallocate (bufrcv_ijk)
      allocate (bufrcv_ijk(three, remapcount_rcv_np))
      if (allocated(bufrcv_qq)) deallocate (bufrcv_qq)
      allocate (bufrcv_qq(remapcount_rcv_np, mtype_remap))

      mmx = three*remapcount_rcv_np
      call mpi_irecv(bufrcv_ijk, mmx, mpi_integer, np, 10 &
           & , mpi_comm_world, ireq_rcv_ijk, merr)

      mmx = mtype_remap*remapcount_rcv_np
      call mpi_irecv(bufrcv_qq, mmx, mpi_real, np, 20 &
           & , mpi_comm_world, ireq_rcv_qq, merr)
    endif

    if (remapcount_snd_np /= 0) then
      call mpi_wait(ireq_snd_ijk, mstatus, merr)
      call mpi_wait(ireq_snd_qq, mstatus, merr)
    endif

    if (remapcount_rcv_np /= 0) then
      call mpi_wait(ireq_rcv_ijk, mstatus, merr)
      call mpi_wait(ireq_rcv_qq, mstatus, merr)

      do remapcount_rcv0 = 1, remapcount_rcv_np
        i = bufrcv_ijk(one, remapcount_rcv0)
        j = bufrcv_ijk(two, remapcount_rcv0)
        k = bufrcv_ijk(three, remapcount_rcv0)
        do m = 1, mtype_remap
          qre(i, j, k, m) = bufrcv_qq(remapcount_rcv0, m)
        enddo
      enddo
    endif

  enddo

  !!====================================
  !!====================================
  !!====================================
  ! boundary condition

  if (iixl_myrank /= 0 .and. jjxl_myrank /= 0) then
    !$omp parallel do private(m,i,j,k)
    do k = 1, kkx
    do j = 1, jjxlg_myrank
    do i = 1, margin
    do m = 1, mtype_remap
      bufsndxd(m, i, j, k) = qre(i + margin, j, k, m)
      bufsndxu(m, i, j, k) = qre(iixlg_myrank - 2*margin + i, j, k, m)
    enddo
    enddo
    enddo
    enddo

    !$omp parallel do private(m,i,j,k)
    do k = 1, kkx
    do j = 1, margin
    do i = 1, iixlg_myrank
    do m = 1, mtype_remap
      bufsndyd(m, i, j, k) = qre(i, j + margin, k, m)
      bufsndyu(m, i, j, k) = qre(i, jjxlg_myrank - 2*margin + j, k, m)
    enddo
    enddo
    enddo
    enddo

    if (ir(myrank) /= 1) then
      np = np_ijr(ir(myrank) - 1, jr(myrank))
      mmx = mtype_remap*margin*jjxlg_myrank*kkx
      call mpi_isend(bufsndxd, mmx, mpi_real, np, 10 &
           & , mpi_comm_world, ireq_sndxd, merr)
      call mpi_irecv(bufrcvxd, mmx, mpi_real, np, 20 &
           & , mpi_comm_world, ireq_rcvxd, merr)
    endif

    if (ir(myrank) /= ixr) then
      np = np_ijr(ir(myrank) + 1, jr(myrank))
      mmx = mtype_remap*margin*jjxlg_myrank*kkx
      call mpi_isend(bufsndxu, mmx, mpi_real, np, 20 &
           & , mpi_comm_world, ireq_sndxu, merr)
      call mpi_irecv(bufrcvxu, mmx, mpi_real, np, 10 &
           & , mpi_comm_world, ireq_rcvxu, merr)
    endif

    if (jr(myrank) /= 1) then
      np = np_ijr(ir(myrank), jr(myrank) - 1)
      mmx = mtype_remap*iixlg_myrank*margin*kkx
      call mpi_isend(bufsndyd, mmx, mpi_real, np, 30 &
           & , mpi_comm_world, ireq_sndyd, merr)
      call mpi_irecv(bufrcvyd, mmx, mpi_real, np, 40 &
           & , mpi_comm_world, ireq_rcvyd, merr)
    else
      !$omp parallel do private(m,i,j,k)
      do k = 1, kkx
      do j = 1, margin
      do i = 1, iixlg_myrank
      do m = 1, mtype_remap
        qre(i, j, k, m) = qre(i, 2*margin - j + 1, k, m)
      enddo
      enddo
      enddo
      enddo
    endif

    if (jr(myrank) /= jxr) then
      np = np_ijr(ir(myrank), jr(myrank) + 1)
      mmx = mtype_remap*iixlg_myrank*margin*kkx
      call mpi_isend(bufsndyu, mmx, mpi_real, np, 40 &
           & , mpi_comm_world, ireq_sndyu, merr)
      call mpi_irecv(bufrcvyu, mmx, mpi_real, np, 30 &
           & , mpi_comm_world, ireq_rcvyu, merr)
    else
      !$omp parallel do private(m,i,j,k)
      do k = 1, kkx
      do j = 1, margin
      do i = 1, iixlg_myrank
      do m = 1, mtype_remap
        qre(i, jjxlg_myrank - margin + j, k, m) = qre(i, jjxlg_myrank - margin - j + 1, k, m)
      enddo
      enddo
      enddo
      enddo
    endif

    if (ir(myrank) /= 1) then
      call mpi_wait(ireq_sndxd, mstatus, merr)
      call mpi_wait(ireq_rcvxd, mstatus, merr)
      !$omp parallel do private(m,i,j,k)
      do k = 1, kkx
      do j = 1, jjxlg_myrank
      do i = 1, margin
      do m = 1, mtype_remap
        qre(i, j, k, m) = bufrcvxd(m, i, j, k)
      enddo
      enddo
      enddo
      enddo
    endif

    if (ir(myrank) /= ixr) then
      call mpi_wait(ireq_sndxu, mstatus, merr)
      call mpi_wait(ireq_rcvxu, mstatus, merr)
      !$omp parallel do private(m,i,j,k)
      do k = 1, kkx
      do j = 1, jjxlg_myrank
      do i = 1, margin
      do m = 1, mtype_remap
        qre(iixlg_myrank - margin + i, j, k, m) = bufrcvxu(m, i, j, k)
      enddo
      enddo
      enddo
      enddo
    endif

    if (jr(myrank) /= 1) then
      call mpi_wait(ireq_sndyd, mstatus, merr)
      call mpi_wait(ireq_rcvyd, mstatus, merr)
      !$omp parallel do private(m,i,j,k)
      do k = 1, kkx
      do j = 1, margin
      do i = 1, iixlg_myrank
      do m = 1, mtype_remap
        qre(i, j, k, m) = bufrcvyd(m, i, j, k)
      enddo
      enddo
      enddo
      enddo
    endif

    if (jr(myrank) /= jxr) then
      call mpi_wait(ireq_sndyu, mstatus, merr)
      call mpi_wait(ireq_rcvyu, mstatus, merr)
      !$omp parallel do private(m,i,j,k)
      do k = 1, kkx
      do j = 1, margin
      do i = 1, iixlg_myrank
      do m = 1, mtype_remap
        qre(i, jjxlg_myrank - margin + j, k, m) = bufrcvyu(m, i, j, k)
      enddo
      enddo
      enddo
      enddo
    endif

    do k = 1, margin
    do j = 1, jjxlg_myrank
    do i = 1, iixlg_myrank
    do m = 1, mtype_remap
      qre(i, j, kkx - margin + k, m) = qre(i, j, margin + k, m)
      qre(i, j, k, m) = qre(i, j, kkx - 2*margin + k - 1, m)
    enddo
    enddo
    enddo
    enddo
  endif

  if (iixl_myrank /= 0 .and. jjxl_myrank /= 0) then
    ist = 1 + margin
    iet = iixlg_myrank - margin
    jst = 1 + margin
    jet = jjxlg_myrank - margin
    kst = 1 + margin
    ket = kkx - margin

    call remap_calc(qre, iixlg_myrank, jjxlg_myrank, kkx, iixl_myrank, jjxl_myrank)
  endif
#endif
#endif

  return
end subroutine remap_yinyang
