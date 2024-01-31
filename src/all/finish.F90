subroutine finish
  use implicit_def
  use comm_def, only: myrank
  implicit none
  
  include "mpif.h"  
  if(myrank == 0) then
     write(*,*) "Finish correctly, great!!"
  endif
  call mpi_barrier(mpi_comm_world,merr)
  call mpi_finalize(merr)
  stop
  
  return
end subroutine finish
  
  
