!========================================================================================|
!
!  R2D2 radiation MHD code
!  Copyright(C) Hideyuki Hotta <hotta.hieyuki@gmail.com>
!  
!========================================================================================|
!----------------------------------------------------------------------------------------|
subroutine comm_init
!----------------------------------------------------------------------------------------|
!
! This is a subroutine for defining the block information.
! Block means the each MPI blocks
! 
!----------------------------------------------------------------------------------------|
!----------------------------------------------------------------------------------------|
  use implicit_def
  use geometry_def, only: &
       & xdcheck,ydcheck,zdcheck &
       & ,ix0,jx0,kx0 &
       & ,iper,jper,kper,npe0,yinyang
  use comm_def, only: xyz,blck,myrank,npe &
       & ,ib,jb,kb,myrnk_loca,mpi_yinyang

#ifdef YinYang
  use yinyang_def, only: mpi_comm_world_yy,myrank_yy,npe_yy
#endif
  
  implicit none
#ifdef fujitsu
  include "mpif-ext.h"
#endif
  include "mpif.h"

  integer :: ndims ! number of MPI process in each direction
  integer :: cart_comm
  integer, dimension(:), allocatable :: dims,xyz0
  logical, dimension(:), allocatable :: periodic
  logical :: reorder
  integer :: np

  !integer :: id,id0,id1
  integer :: ista,iend

  integer :: key,color
  logical :: origin_rank_flag
  integer :: mpi_comm_world_tmp,myrank_tmp,npe_tmp
!----------------------------------------------------------------------------------------|

#ifdef fujitsu
  integer :: fjmpi_dim
  integer :: fjmpi_ix0,fjmpi_jx0,fjmpi_kx0
  integer, dimension(3) :: fjmpi_coords
  integer, parameter :: fjmpi_maxppn = 4
  integer, dimension(fjmpi_maxppn) :: fjmpi_ranks
  integer :: fjmpi_outppn
  integer :: fjmpi_cmg_order
  integer :: ibl,jbl,kbl
  integer :: cmg_x = 1, cmg_y = 2, cmg_z = 2

#ifdef YinYang
  integer, dimension(0:npe/2-1) :: fjmpi_x,fjmpi_y,fjmpi_z
  integer, dimension(0:npe-1) :: rank2rank_yy
  integer, dimension(0:npe-1) :: mpi_yinyang_global
  integer, dimension(0:npe-1,3) :: tmm
#else
  integer, dimension(0:npe-1) :: fjmpi_x,fjmpi_y,fjmpi_z
#endif


#else ! not fujitsu  

#ifdef YinYang
  integer, dimension(0:npe/2-1) :: fjmpi_x,fjmpi_y,fjmpi_z
  integer, dimension(0:npe-1) :: rank2rank_yy,rank2rank_yy_local
  integer, dimension(0:npe-1) :: mpi_yinyang_global,mpi_yinyang_global_local
  integer, dimension(0:npe-1,3) :: tmm
#else
  integer, dimension(0:npe-1) :: fjmpi_x,fjmpi_y,fjmpi_z
#endif
  
  
#endif
!----------------------------------------------------------------------------|

  ! dimension check
  if(xdcheck*ydcheck*zdcheck == 2) ndims = 1 ! 1D problem
  if(xdcheck*ydcheck*zdcheck == 4) ndims = 2 ! 2D problem
  if(xdcheck*ydcheck*zdcheck == 8) ndims = 3 ! 3D problem

! Fugaku and Flow topology
#ifdef fujitsu
  ! if this is not 3D problem tofu topology is not used
  if(ndims /= 3) goto 1000

  !----------------------------
  ! Check tofu logical dimension 
  call fjmpi_topology_get_dimension(fjmpi_dim,merr)
  if(fjmpi_dim /= 3) then
     if(myrank == 0) then
        write(*,*) 'node specification is not 3D'
        write(*,*) 'thus move to normal MPI coordinate'
     endif
     goto 1000
  endif

  !----------------------------
  ! Check tofu logical shape
  call fjmpi_topology_get_shape(fjmpi_ix0,fjmpi_jx0,fjmpi_kx0,merr)
#ifdef YinYang
  if((fjmpi_ix0*cmg_x /= ix0) .or. (fjmpi_jx0*cmg_y/2 /= jx0) .or. (fjmpi_kx0*cmg_z /= kx0)) then
#else
  if((fjmpi_ix0*cmg_x /= ix0) .or. (fjmpi_jx0*cmg_y /= jx0) .or. (fjmpi_kx0*cmg_z /= kx0)) then
#endif
     if(myrank == 0) then
        write(*,*) '3D tofu coordinate is different from code intenstion'
        write(*,*) 'thus move to normal MPI coordinate'
     endif
     goto 1000
  endif

  !----------------------------
  ! Check coordinate of myrank
  call fjmpi_topology_get_coords(mpi_comm_world,myrank,fjmpi_logical,fjmpi_dim,fjmpi_coords,merr)

  !----------------------------
  ! Check 
  call fjmpi_topology_get_ranks(mpi_comm_world,fjmpi_logical,fjmpi_coords,fjmpi_maxppn,fjmpi_outppn,fjmpi_ranks,merr)
  do m = 1,4
     if(myrank == fjmpi_ranks(m)) fjmpi_cmg_order = m - 1
  enddo

  kbl = fjmpi_cmg_order/cmg_x/cmg_y
  jbl = (fjmpi_cmg_order - cmg_x*cmg_y*kbl)/cmg_x
  ibl = fjmpi_cmg_order - jbl*cmg_x - kbl*cmg_x*cmg_y
    
  ib = fjmpi_coords(1)*cmg_x + ibl 
  jb = fjmpi_coords(2)*cmg_y + jbl
  kb = fjmpi_coords(3)*cmg_z + kbl

#ifdef YinYang
  if (jb < jx0) then
     mpi_yinyang = 0
  else
     mpi_yinyang = 1
     jb = jb - jx0
  endif
!  if (kb < kx0) then
!     mpi_yinyang = 0
!  else
!     mpi_yinyang = 1
!     kb = kb - kx0
!  endif
#endif
  
  if(myrank == 0) then
     write(*,*) 'We use Tofu oriented MPI topology'
  endif

#ifdef YinYang
  ! split mpi_comm_world for Yin-Yang
  key = 0
  if (mpi_yinyang == 0) then
     color = 0
  else
     color = 1
  endif
  call mpi_comm_split(mpi_comm_world,color,key,mpi_comm_world_yy,merr)
  call mpi_comm_size(mpi_comm_world_yy,npe_yy   ,merr)
  call mpi_comm_rank(mpi_comm_world_yy,myrank_yy,merr)  

  fjmpi_x(myrank_yy) = ib
  fjmpi_y(myrank_yy) = jb
  fjmpi_z(myrank_yy) = kb

  call mpi_allgather(fjmpi_x(myrank_yy),1,mpi_integer,fjmpi_x(0),1,mpi_integer,mpi_comm_world_yy,merr)
  call mpi_allgather(fjmpi_y(myrank_yy),1,mpi_integer,fjmpi_y(0),1,mpi_integer,mpi_comm_world_yy,merr)
  call mpi_allgather(fjmpi_z(myrank_yy),1,mpi_integer,fjmpi_z(0),1,mpi_integer,mpi_comm_world_yy,merr)
#else
  fjmpi_x(myrank) = ib
  fjmpi_y(myrank) = jb
  fjmpi_z(myrank) = kb

  call mpi_allgather(fjmpi_x(myrank),1,mpi_integer,fjmpi_x(0),1,mpi_integer,mpi_comm_world,merr)
  call mpi_allgather(fjmpi_y(myrank),1,mpi_integer,fjmpi_y(0),1,mpi_integer,mpi_comm_world,merr)
  call mpi_allgather(fjmpi_z(myrank),1,mpi_integer,fjmpi_z(0),1,mpi_integer,mpi_comm_world,merr)
#endif
  
  xyz(:,1) = fjmpi_x
  xyz(:,2) = fjmpi_y
  xyz(:,3) = fjmpi_z
 
  goto 2000
#endif

1000 continue

  allocate(dims(ndims)) ! number of MPI process in each direction
  allocate(periodic(ndims)) ! periodic information for MPI_CART_CREATE
  allocate(xyz0(ndims)) !

  xyz(:,:) = 0
  reorder = .true.
  periodic(:) = .true.

  ! 1D
  if(ndims == 1) then
     if(xdcheck /= 1) then
        if(iper == 0) periodic(1) = .false.
        dims(1) = ix0
     endif
     if(ydcheck /= 1) then 
        if(jper == 0) periodic(1) = .false.
        dims(1) = jx0
     endif
     if(zdcheck /= 1) then
        if(kper == 0) periodic(1) = .false.
        dims(1) = kx0
     endif
  endif

  ! 2D
  if(ndims == 2) then
     if(xdcheck == 1) then
        if(jper == 0) periodic(1) = .false.
        if(kper == 0) periodic(2) = .false.
        dims(1) = jx0
        dims(2) = kx0
     endif

     if(ydcheck == 1) then
        if(iper == 0) periodic(1) = .false.
        if(kper == 0) periodic(2) = .false.
        dims(1) = ix0
        dims(2) = kx0
     endif

     if(zdcheck == 1) then
        if(iper == 0) periodic(1) = .false.
        if(jper == 0) periodic(2) = .false.
        dims(1) = ix0
        dims(2) = jx0
     endif
  endif

  ! 3D
  if(ndims == 3) then
     if(iper == 0) periodic(1) = .false.
     if(jper == 0) periodic(2) = .false.
     if(kper == 0) periodic(3) = .false.
     dims(1) = ix0
     dims(2) = jx0
     dims(3) = kx0
  endif

  origin_rank_flag = .false.
#ifdef YinYang
     ! split mpi_comm_world for Yin-Yang
     key = 0
     if(myrank < npe/2) then
        mpi_yinyang = 0
        color = 0
     else
        mpi_yinyang = 1
        color = 1
     endif
     call mpi_comm_split(mpi_comm_world,color,key,mpi_comm_world_yy,merr)
     call mpi_comm_size(mpi_comm_world_yy,npe_yy   ,merr)
     call mpi_comm_rank(mpi_comm_world_yy,myrank_yy,merr)
     mpi_comm_world_tmp = mpi_comm_world_yy
     myrank_tmp = myrank_yy
     npe_tmp = npe/2
     if(myrank == 0 .or. myrank == npe/2) origin_rank_flag = .true.     
#else
     mpi_comm_world_tmp = mpi_comm_world
     myrank_tmp = myrank
     npe_tmp = npe
     if(myrank == 0) origin_rank_flag = .true.
#endif
    
  call mpi_cart_create(mpi_comm_world_tmp,ndims,dims,periodic,reorder,cart_comm,merr)

  if(origin_rank_flag) then
     do np = 0,npe_tmp-1
        call mpi_cart_coords(cart_comm,np,ndims,xyz0,merr)

        ! 1D
        if(ndims == 1) then
           if(xdcheck /= 1) xyz(np+1,1) = xyz0(1)
           if(ydcheck /= 1) xyz(np+1,2) = xyz0(1)
           if(zdcheck /= 1) xyz(np+1,3) = xyz0(1)
        endif

        ! 2D
        if(ndims == 2) then
           if(xdcheck == 1) then
              xyz(np+1,2) = xyz0(1)
              xyz(np+1,3) = xyz0(2)
           endif

           if(ydcheck == 1) then
              xyz(np+1,1) = xyz0(1)
              xyz(np+1,3) = xyz0(2)
           endif

           if(zdcheck == 1) then
              xyz(np+1,1) = xyz0(1)
              xyz(np+1,2) = xyz0(2)
           endif
        endif

        ! 3D
        if(ndims == 3) then
           xyz(np+1,:) = xyz0(:)
        endif
     enddo

  endif

  call mpi_bcast(xyz,npe_tmp*3,mpi_integer,0,mpi_comm_world_tmp,merr)

2000 continue

  ! Output background variables
  if(myrank == 0) then
     open(idf,file="./data/param/xyz.dac",form="unformatted")
     write(idf) xyz
     close(idf)     
  endif

#ifdef YinYang

#ifdef fujitsu  
  rank2rank_yy(myrank) = myrank_yy
  mpi_yinyang_global(myrank) = mpi_yinyang
  call mpi_allgather(rank2rank_yy(myrank),1,mpi_integer,rank2rank_yy(0),1,mpi_integer,mpi_comm_world,merr)
  call mpi_allgather(mpi_yinyang_global(myrank),1,mpi_integer,mpi_yinyang_global(0),1,mpi_integer,mpi_comm_world,merr)
#else
  rank2rank_yy_local(myrank) = myrank_yy
  mpi_yinyang_global_local(myrank) = mpi_yinyang
  call mpi_allgather(rank2rank_yy_local(myrank),1,mpi_integer,rank2rank_yy(0),1,mpi_integer,mpi_comm_world,merr)
  call mpi_allgather(mpi_yinyang_global_local(myrank),1,mpi_integer,mpi_yinyang_global(0),1,mpi_integer,mpi_comm_world,merr)
#endif  

  do np = 0,npe-1
     blck%iloca(np,1:3) = xyz(rank2rank_yy(np)+1,1:3)
  enddo

#else  
  blck%iloca(0:npe0-1,1:3) = xyz(1:ix0*jx0*kx0,1:3)
#endif
  
  ib = blck%iloca(myrank,1)
  jb = blck%iloca(myrank,2)
  kb = blck%iloca(myrank,3)
  


!----------------------------------------------------------------------------------------|
  ! when we put the location the function reply ID
  allocate(myrnk_loca(yinyang,-1:ix0,-1:jx0,-1:kx0))
   
  ! When we put the location the function reply ID
  ! -1: physical boundary
  myrnk_loca(:,:,:,:) = -1
  
  ! if yinyang = 1 without Yin-Yang grid
  ! if yinyang = 2 with    Yin-Yang grid
  !id0 = 0
#ifdef YinYang
  do np = 0,npe-1
     if( mpi_yinyang_global(np) == 0) then! for Yin grid
        m = 1
     else
        m = 2
     endif
     i = blck%iloca(np,1)
     j = blck%iloca(np,2)
     k = blck%iloca(np,3)
     myrnk_loca(m,i,j,k) = np
  enddo

#else
  do np = 0,npe-1
     i = blck%iloca(np,1)
     j = blck%iloca(np,2)
     k = blck%iloca(np,3)
     myrnk_loca(1,i,j,k) = np
  enddo
#endif

!  do m = 1,yinyang
!     if(m == 1) then ! for Yin grid
!        id0 = 0
!        id1 = npe0/yinyang-1
!     endif
!     if(m == 2) then ! for Yang grid
!        id0 = npe0/yinyang
!        id1 = npe0-1
!     endif
!     do id = id0,id1
!        i = blck%iloca(id,1)
!        j = blck%iloca(id,2)
!        k = blck%iloca(id,3)
!        myrnk_loca(m,i,j,k) = id
!     enddo

!----------------------------------------------------------------------------------------|
  ! search neighbor block
!  do id = id0,id1
  do np = 0,npe-1
#ifdef YinYang
     if( mpi_yinyang_global(np) == 0 ) then! for Yin grid
        m = 1
     else ! for Yang grid
        m = 2
     endif
#else
     m = 1
#endif

     i = blck%iloca(np,1)
     j = blck%iloca(np,2)
     k = blck%iloca(np,3)

     blck%ineib(np,1) = myrnk_loca(m,i+1,j,k)
     blck%ineib(np,2) = myrnk_loca(m,i-1,j,k)
     blck%ineib(np,3) = myrnk_loca(m,i,j+1,k)
     blck%ineib(np,4) = myrnk_loca(m,i,j-1,k)
     blck%ineib(np,5) = myrnk_loca(m,i,j,k+1)
     blck%ineib(np,6) = myrnk_loca(m,i,j,k-1)
        
     ! for periodic boundary condition
     if(iper == 1) then
        if(blck%ineib(np,1) == -1) blck%ineib(np,1) = myrnk_loca(m,0    ,j,k)
        if(blck%ineib(np,2) == -1) blck%ineib(np,2) = myrnk_loca(m,ix0-1,j,k)
     endif
        
     if(jper == 1) then
        if(blck%ineib(np,3) == -1) blck%ineib(np,3) = myrnk_loca(m,i,0    ,k)
        if(blck%ineib(np,4) == -1) blck%ineib(np,4) = myrnk_loca(m,i,jx0-1,k)
     endif
        
     if(kper == 1) then
        if(blck%ineib(np,5) == -1) blck%ineib(np,5) = myrnk_loca(m,i,j,0    )
        if(blck%ineib(np,6) == -1) blck%ineib(np,6) = myrnk_loca(m,i,j,kx0-1)
     endif
  enddo
!----------------------------------------------------------------------------------------|

  return
end subroutine comm_init
