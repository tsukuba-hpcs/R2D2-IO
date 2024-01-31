!========================================================================================|
!
!  R2D2 radiation MHD code
!  Copyright(C) Hideyuki Hotta <hotta.hieyuki@gmail.com>
!  
!========================================================================================|
module const_def
  ! define double precision
  real(KIND(0.d0)), parameter :: pi = 3.141592653589d0
  real(KIND(0.d0)), parameter :: pii = 1.d0/pi

  real(KIND(0.d0)), parameter :: pii4 = 1.d0/4.d0/pi
  real(KIND(0.d0)), parameter :: pii8 = 1.d0/8.d0/pi
  real(KIND(0.d0)), parameter :: sb = 5.67d-5 ! Stefan-Boltzman constant

#ifdef ideal
!  real(KIND(0.d0)), parameter :: gm = 1.4d0 ! ratio of heat capacities
!  real(KIND(0.d0)), parameter :: gm = 2.d0 ! ratio of heat capacities
  real(KIND(0.d0)), parameter :: gm = 5.d0/3.d0 ! ratio of heat capacities
  real(KIND(0.d0)), parameter :: rr = 1.d0  ! gas constant
  real(KIND(0.d0)), parameter :: cv = rr/(gm - 1.d0)  ! heat capcity at constant volume
  real(KIND(0.d0)), parameter :: theta = 10.d0 ! temperature gradient
  real(KIND(0.d0)), parameter :: m_ad = 1.d0 ! adiabatic index
  real(KIND(0.d0)), parameter :: k_cond = 0.07d0 ! thermal conduction coefficient
  real(KIND(0.d0)), parameter :: sigma = 1.d0 ! Prandtl number
  real(KIND(0.d0)), parameter :: m_visc = (gm-1)/gm/rr*sigma*k_cond ! viscosity coefficient
#endif  
  
end module const_def
