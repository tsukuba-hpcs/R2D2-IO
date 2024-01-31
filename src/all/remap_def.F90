!========================================================================================|
!
!  R2D2 radiation MHD code
!  Copyright(C) Hideyuki Hotta <hotta.hieyuki@gmail.com>
!  
!========================================================================================|
module remap_def
#ifndef non_IO
  use geometry_def, only: ix0,jx0,kx0,nx,ny,nz,marginx,marginy,marginz,npe0,mtype
  include "mpif.h"

  !================================================================  
  ! for remaping
  integer, parameter :: iix0 = ix0*nx
#ifdef YinYang
  integer, parameter :: jjx0 = jx0*ny*2
  integer, parameter :: kkx0 = jx0*ny*4  
#else  
  integer, parameter :: jjx0 = jx0*ny
  integer, parameter :: kkx0 = kx0*nz
#endif  
  integer, parameter :: iix = iix0 + 2*marginx
  integer, parameter :: jjx = jjx0 + 2*marginy
  integer, parameter :: kkx = kkx0 + 2*marginz
  real(KIND=KIND(0.d0)), dimension(iix),save :: xt ! transformed ordinal geometry
  real(KIND=KIND(0.d0)), dimension(jjx),save :: yt ! transformed ordinal geometry
  real(KIND=KIND(0.d0)), dimension(kkx),save :: zt ! transformed ordinal geometry
  integer, dimension(iix),save :: i2ir
  integer, dimension(jjx),save :: j2jr
  
  integer, dimension(0:npe0-1), save :: ir,jr,iixl,jjxl,iixlg,jjxlg
  integer, dimension(0:npe0-1), save :: is ,ie ,js ,je  ! start and end without margin
  integer, dimension(0:npe0-1), save :: is0,ie0,js0,je0 ! start and end with margin

  integer, dimension(:,:), allocatable, save :: np_ijr
  integer, dimension(0:npe0-1), save :: remapcount_snd,remapcount_rcv
  integer, dimension(6,0:npe0-1), save :: remapijk_snd,remapijk_rcv
  real(KIND=KIND(0.e0)), dimension(:,:,:,:), allocatable, save :: qre  
  integer, save :: icount
  integer, save :: ixr,jxr
  
  integer, parameter :: mtype_remap = mtype + 2

  integer, save :: iixl_myrank,jjxl_myrank,iixlg_myrank,jjxlg_myrank
  integer, dimension(:,:), allocatable, save :: jkc
  integer, parameter :: m2d_spex = 14
  integer, public, save :: mpi_comm_world_spex, mpi_comm_world_spex_io
  integer, public, save :: myrank_spex

  ! MPI IO for values.F95
  integer, dimension(3), public, save :: gsize,ssize,start
  integer, dimension(4), public, save :: gsize_in,ssize_in,start_in
  integer, dimension(3), public, save :: gsize_2d,ssize_2d,start_2d
  integer, dimension(2), public, save :: gsize0,ssize0,start0
  integer, public, save :: fha,fhtop,fh2d
  integer(Kind=MPI_OFFSET_KIND), public, parameter :: disp = 0
  ! ***_a for ordinal average or correlation
  ! ***_b for spherical harmonic expansion
  integer, public, save :: global_xy,global_xz,global_flux,global_spex,global_top,global_2d
  integer, save :: mpi_comm_world_remap  ! defined at remap_init.F90 (for remap output)
  integer, save :: mpi_comm_world_remap_xz ! defined ad remap_init.F90 (for remap xz slice)
  integer, save :: mpi_comm_world_ir     ! defined at remap_init.F90 (for ir thred)
  integer, save :: mpi_comm_world_top_yy ! defined at remap_yinyang_init.F90 (for top output)

  ! slice position
  integer, parameter :: jc = ny*jx0/2 ! for xz slice, sometimes used in model_*** i.e., initial condition
  integer, parameter :: kc = nz*kx0/2 ! for xy slice, sometimes used in model_*** i.e., initial condition

  integer, parameter :: m_in = 13 ! output for tau = 1
  integer, parameter :: m_tu = 3  ! number of optical depth for output
  integer, parameter :: m2d_xy = 100
  integer, parameter :: m2d_xz = 13
  integer, parameter :: m2d_flux = 5
contains

!========================================================================================|
!========================================================================================|
  subroutine cor(qq1,qq1m,qq2,qq2m,qq1qq2co,ix,jx,kx,margin)
!========================================================================================|
    implicit none
    
    integer :: i,j,k,ix,jx,kx,margin
    real(KIND(0.e0)), dimension(ix,jx,kx) :: qq1,qq2
    real(KIND(0.e0)), dimension(ix,jx) :: qq1m,qq2m,qq1qq2co
!----------------------------------------------------------------------------------------|

    do i = 1,ix
    do j = 1,jx
       qq1qq2co(i,j) = 0.d0
       do k = 1+margin,kx-margin
          qq1qq2co(i,j) = qq1qq2co(i,j) + (qq1(i,j,k)-qq1m(i,j))*(qq2(i,j,k)-qq2m(i,j))
       enddo
       qq1qq2co(i,j) = qq1qq2co(i,j)/real(kx-2*margin)
    enddo
    enddo
    
    return
  end subroutine cor

!========================================================================================|
  subroutine pdf2d(qq,jx,kx,qqh,nh,ww)
!========================================================================================|
! qq[jx,kx]: real(kind(0.e0)): variable
! ww[jx,kx]: real(kind(0.e0)): weight
! jx,kx: integer: dimension size
!
! qqh[nh]: real(kind(0.e0)): PDF
! nh: integer: bin number of PDF
! 
!----------------------------------------------------------------------------------------|
    implicit none

    integer :: j,k,jx,kx
    integer :: n
    integer :: nh
    real(KIND=KIND(0.e0)) :: sumq
    real(KIND=KIND(0.e0)) :: minq,maxq,binq
    real(KIND=KIND(0.e0)), dimension(jx,kx) :: qq,ww
    real(KIND=KIND(0.e0)), dimension(nh) :: qqh
!----------------------------------------------------------------------------------------|

    maxq = maxval(qq)
    minq = minval(qq)
    binq = (maxq-minq)/real(nh-3)
    qqh(1:nh) = 0.d0

    if(binq /= 0.d0) then
       do k = 1,kx
       do j = 1,jx
          n = int((qq(j,k)-minq)/binq)+1
          if(n < 0 .or. n > nh) write(*,*) n,binq
          qqh(n) = qqh(n) + ww(j,k)
       enddo
       enddo
       sumq = 0.d0
       do n = 1,nh-2
          sumq = sumq + qqh(n)
       enddo
       qqh = qqh/sumq/binq
       qqh(nh-1) = maxq
       qqh(nh) = minq
    else
       qqh(1:nh) = 0.d0
    endif

    return
  end subroutine pdf2d

!========================================================================================|
#endif
end module remap_def
