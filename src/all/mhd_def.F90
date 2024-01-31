!========================================================================================|
!
!  R2D2 radiation MHD code
!  Copyright(C) 2019 Hideyuki Hotta <hotta@chiba-u.jp>
!  
!========================================================================================|
module mhd_def
  use geometry_def, only: mtype,nxg,nyg,nzg,ix0

  ! basic variables, ro,vx,vy,vz,bx,by,bz,se( or en)
  real(KIND(0.d0)), dimension(nxg,nyg,nzg,mtype), save :: qq,qqm
  real(KIND(0.d0)), dimension(nxg,nyg,nzg,mtype+3), save :: qqp
end module mhd_def
