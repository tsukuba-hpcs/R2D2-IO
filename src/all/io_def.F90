!========================================================================================|
!
!  R2D2 radiation MHD code
!  Copyright(C) Hideyuki Hotta <hotta.hieyuki@gmail.com>
!  
!========================================================================================|
module io_def
!========================================================================================|

  ! variables for MPI I/O
  logical, parameter, public :: memory_saving = .false.
  integer, dimension(4), public, save :: gsize_io,ssize_io,start_io
  integer, public, save :: global_io,fh_io,mmx_io

  integer, dimension(3), public, save :: gsize_io_memory,ssize_io_memory,start_io_memory
  integer, public, save :: global_io_memory,fh_io_memory,mmx_io_memory

  integer, public, save :: mpi_comm_world_jk ! refined at io_init.F90 (for allreduce for RTE)

  private
  interface io_params_write
     module procedure  io_params_write_d, io_params_write_i &
          & , io_params_write_c, io_params_write_l
  end interface io_params_write

  public :: io_params_write

contains

  !---------------------------------------------------
  subroutine io_params_write_d(idf,param,char_param)
    implicit none
    
    integer, intent(in) :: idf
    real(KIND(0.d0)), intent(in) :: param
    character(*), intent(in) :: char_param

    write(idf,'(e16.10,1x,A,1x,A)') param,char_param,"d"
    
    return
  end subroutine io_params_write_d

  !---------------------------------------------------
  subroutine io_params_write_i(idf,param,char_param)
    implicit none
    
    integer, intent(in) :: idf
    integer, intent(in) :: param
    character(*), intent(in) :: char_param

    write(idf,'(i8,1x,A,1x,A)') param,char_param,"i"
    
    return
  end subroutine io_params_write_i

  !---------------------------------------------------
  subroutine io_params_write_c(idf,param,char_param)
    implicit none
    
    integer, intent(in) :: idf
    character(*), intent(in) :: param
    character(*), intent(in) :: char_param

    write(idf,'(A,1x,A,1x,A)') param,char_param,"c"
    
    return
  end subroutine io_params_write_c

  !---------------------------------------------------
  subroutine io_params_write_l(idf,param,char_param)
    implicit none
    
    integer, intent(in) :: idf
    logical, intent(in) :: param
    character(*), intent(in) :: char_param

    write(idf,'(L1,1x,A,1x,A)') param,char_param,"l"
    
    return
  end subroutine io_params_write_l
  
!========================================================================================|
  
end module io_def
