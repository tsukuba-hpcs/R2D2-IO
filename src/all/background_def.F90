!========================================================================================|
!
!  R2D2 radiation MHD code
!  Copyright(C) Hideyuki Hotta <hotta.hieyuki@gmail.com>
!  
!========================================================================================|
! Define parameters for 3D full radiation transfer
module background_def
  use const_def, only: pi
  use geometry_def, only: ix00,nxg,nyg,nzg

  ! background variables in global view
  real(KIND(0.d0)), dimension(ix00), save :: ro00,pr00,te00,se00,en00
  real(KIND(0.d0)), dimension(ix00), save :: op00,tu00
  real(KIND(0.d0)), dimension(ix00), save :: dsedr00,dtedr00
  real(KIND(0.d0)), dimension(ix00), save :: dprdro00,dprdse00
  real(KIND(0.d0)), dimension(ix00), save :: dtedro00,dtedse00
  real(KIND(0.d0)), dimension(ix00), save :: dendro00,dendse00
  real(KIND(0.d0)), dimension(ix00), save :: gx00,cp00,ma00
  real(KIND(0.d0)), dimension(ix00), save :: xi00
  ! artificial flux and artificial source term
  real(KIND(0.d0)), dimension(ix00), save :: fa00,sa00,la00

  ! background variables in local view
  real(KIND(0.d0)), dimension(nxg), save :: ro0,pr0,te0,se0,en0,se0x
  real(KIND(0.d0)), dimension(nxg), save :: ro0o,te0o,se0o
  real(KIND(0.d0)), dimension(nxg), save :: op0,tu0
  real(KIND(0.d0)), dimension(nxg), save :: dtedr,dsedr
  real(KIND(0.d0)), dimension(nxg), save :: dprdro,dprdse,dtedro,dtedse,dendro,dendse
  real(KIND(0.d0)), dimension(nxg), save :: gx,cp,xi,xi2o
  real(KIND(0.d0)), dimension(nxg), save :: sa0

  ! This parameter controls the treatment of **0 + **1
  real(KIND(0.d0)), save :: backfac

  ! for artificial viscosity
  real(KIND(0.d0)), dimension(nxg), save :: xi2m

  ! Rotation
  real(KIND(0.d0)), parameter :: omfac = 0.d0 ! ratio to the solar rotation
  real(KIND(0.d0)), parameter :: om00 = 2.d0*pi*413.d-9*omfac
  real(KIND(0.d0)), parameter :: theta = 70.d0/180.d0*pi
#ifdef spherical
  real(KIND(0.d0)), dimension(nyg,nzg), save :: omx,omy,omz
#else
  real(KIND(0.d0)), parameter :: omx = +om00*cos(theta)
  real(KIND(0.d0)), parameter :: omy = -om00*sin(theta)
  real(KIND(0.d0)), parameter :: omz = 0.d0
#endif
  
  ! divb cleaning time scale
  real(KIND(0.d0)), dimension(nxg) :: tauo_divb

contains
!========================================================================================|
  subroutine interpolate_background(qq_int,qq_real,r_int,r_real,ix_int,ix_real)
!========================================================================================|
!----------------------------------------------------------------------------------------|
! This is a subroutine for linear interpolation used in
! background_init.f90
!----------------------------------------------------------------------------------------|
    
    implicit none
    integer, intent(in) :: ix_int
    integer :: i
    integer :: ii,ir
    integer :: ix_real

    real(KIND(0.d0)) :: dxx, dxx1, dxx2
    real(KIND(0.d0)), dimension(ix_int ) :: qq_int , r_int
    real(KIND(0.d0)), dimension(ix_real) :: qq_real, r_real
 !----------------------------------------------------------------------------------------|
  
    ii = 1
    do i = 1,ix_int
       do ir = 1,ix_real
          ii = ir
          dxx = r_int(i) - r_real(ir)
          if(dxx < 0.) goto 1000
       enddo
1000   continue
       if(ii <= 1) ii = 2
       dxx1 = r_real(ii) - r_int(i)
       dxx2 = r_int(i) - r_real(ii-1)
       qq_int(i) = (qq_real(ii)*dxx2 + qq_real(ii-1)*dxx1)/(dxx1 + dxx2)
    enddo    
    return
  end subroutine interpolate_background
!-------------------------------------------------------------|

end module background_def
