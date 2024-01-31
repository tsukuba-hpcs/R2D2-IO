!========================================================================================|
!
!  R2D2 radiation MHD code
!  Copyright(C) Hideyuki Hotta <hotta.hieyuki@gmail.com>
!  
!========================================================================================|
module comm_def
  use geometry_def, only: ix0,jx0,kx0,npe0
  ! geometry array for MPI thread
  ! Myrank in MPI thread
  integer, save :: myrank,npe

  integer, save :: mpi_yinyang
  integer, dimension(ix0*jx0*kx0,3), save :: xyz  

  !Structure of the ID information
  type id_info
     sequence
     integer, dimension(0:npe0-1,3) :: iloca
     integer, dimension(0:npe0-1,6) :: ineib
  end type id_info
  type(id_info), save :: blck
  integer, save :: ib,jb,kb
  integer, dimension(:,:,:,:), allocatable, save :: myrnk_loca
  ! enter ID and return the location. Required to search neibor and Yin-Yang connection.

end module comm_def
