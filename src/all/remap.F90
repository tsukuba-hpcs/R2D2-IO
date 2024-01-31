!----------------------------------------------------------------------|
subroutine remap(tu,fr)
!----------------------------------------------------------------------|
! Thie program remaps each MPI threads to Global data
!======================================================================|
  use const_def, only: pi
  use geometry_def, only: ix0,jx0,kx0,nx,ny,nz,nxg,nyg,nzg,margin,mtype
  use implicit_def

  use comm_def, only: npe,myrank,ib,jb,kb
  use mhd_def, only: qq
#ifndef non_IO  
  use remap_def, only: qre&
       & ,remapcount_snd,remapcount_rcv &
       & ,remapijk_snd,remapijk_rcv &
       & ,is0,js0 &
       & ,iixl,jjxl &
       & ,iixlg,jjxlg &
       & ,iix,jjx,kkx &
       & ,mtype_remap
#endif  
  implicit none
  include "mpif.h"  
  integer, dimension(mpi_status_size) :: mstatus_rcv,mstatus_snd

  real(KIND=KIND(0.d0)), dimension(nxg,nyg,nzg) :: tu,fr
#ifndef non_IO  
  integer :: np
  integer :: ii0,ii1,iii0,iii1,iii2,iii
  integer :: jj0,jj1,jjj0,jjj1,jjj2,jjj
  integer :: kk0,kk1,kkk0,kkk1,kkk2,kkk
  real(KIND=KIND(0.e0)), dimension(:,:,:,:), allocatable :: bufsnd_qq,bufrcv_qq
  integer :: mmx
  integer :: ireq_snd_qq,ireq_rcv_qq
  integer :: jst,jet,kst,ket
  
  !======================================================================|
  ! avoid compiler warning
  iii2 = 0
  jjj2 = 0
  kkk2 = 0

  iii0 = 0
  jjj0 = 0
  kkk0 = 0

  !!!----------------------------!
  do np = 0,npe-1
     if(remapcount_snd(np) /= 0) then
        iii0 = remapijk_snd(1,np)
        jjj0 = remapijk_snd(2,np)
        kkk0 = remapijk_snd(3,np)

        iii1 = remapijk_snd(4,np)
        jjj1 = remapijk_snd(5,np)
        kkk1 = remapijk_snd(6,np)

        iii2 = iii1 - iii0 + 1
        jjj2 = jjj1 - jjj0 + 1
        kkk2 = kkk1 - kkk0 + 1

        ii0 = iii0 - ib*nx
        ii1 = iii1 - ib*nx

        jj0 = jjj0 - jb*ny
        jj1 = jjj1 - jb*ny

        kk0 = kkk0 - kb*nz
        kk1 = kkk1 - kb*nz

        if(allocated(bufsnd_qq)) deallocate(bufsnd_qq)
        allocate(bufsnd_qq(iii2,jjj2,kkk2,mtype_remap))
        bufsnd_qq(1:iii2,1:jjj2,1:kkk2,1:mtype) = real(qq(ii0:ii1,jj0:jj1,kk0:kk1,1:mtype))
        
        bufsnd_qq(1:iii2,1:jjj2,1:kkk2,mtype+1) = sngl(fr(        ii0:ii1,jj0:jj1,kk0:kk1))
        bufsnd_qq(1:iii2,1:jjj2,1:kkk2,mtype+2) = sngl(tu(        ii0:ii1,jj0:jj1,kk0:kk1))

        mmx = mtype_remap*iii2*jjj2*kkk2
        call mpi_isend(bufsnd_qq,mmx,mpi_real,np,10,mpi_comm_world,ireq_snd_qq,merr)
     endif
     
     if(remapcount_rcv(np) /= 0) then
        iii0 = remapijk_rcv(1,np)
        jjj0 = remapijk_rcv(2,np)
        kkk0 = remapijk_rcv(3,np)

        iii1 = remapijk_rcv(4,np)
        jjj1 = remapijk_rcv(5,np)
        kkk1 = remapijk_rcv(6,np)

        iii2 = iii1 - iii0 + 1
        jjj2 = jjj1 - jjj0 + 1
        kkk2 = kkk1 - kkk0 + 1

        if(allocated(bufrcv_qq)) deallocate(bufrcv_qq)
        allocate(bufrcv_qq(iii2,jjj2,kkk2,mtype_remap))
        
        mmx = mtype_remap*iii2*jjj2*kkk2
        call mpi_irecv(bufrcv_qq,mmx,mpi_real,np,10,mpi_comm_world,ireq_rcv_qq,merr)        
     endif

     if(remapcount_snd(np) /= 0) then
        call mpi_wait(ireq_snd_qq,mstatus_snd,merr)
     endif

     if(remapcount_rcv(np) /= 0) then 
        call mpi_wait(ireq_rcv_qq,mstatus_rcv,merr)
        
        do k = 1,kkk2
        do j = 1,jjj2
        do i = 1,iii2
        do m = 1,mtype_remap
           iii = i - 1 + iii0 - is0(myrank) + 1
           jjj = j - 1 + jjj0 - js0(myrank) + 1
           kkk = k - 1 + kkk0
           qre(iii,jjj,kkk,m) = bufrcv_qq(i,j,k,m)
        enddo
        enddo
        enddo
        enddo
     endif
  enddo
  
  
  !-----------------------------
  if(iixl(myrank) /= 0 .and. jjxl(myrank) /= 0) then
     ist = 1 + margin
     iet = iixlg(myrank) - margin
     jst = 1 + margin
     jet = jjxlg(myrank) - margin
     kst = 1 + margin
     ket = kkx - margin

     call remap_calc(qre,iixlg(myrank),jjxlg(myrank),kkx,iixl(myrank),jjxl(myrank))
  endif

#endif
  return
end subroutine remap
