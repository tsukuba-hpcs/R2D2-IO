!----------------------------------------------------------------------|
! prepare Yin-Yang boundary
subroutine yinyang_init
#ifdef YinYang
  use geometry_def, only: nxg, nyg ,nzg, nx, ny, nz &
       & , xdcheck, ydcheck, zdcheck, mtype &
       & , marginx, marginy, marginz,margin &
       & , xmin,xmax,ymin,ymax,zmin,zmax &
       & , ix0,jx0,kx0 &
       & ,i1,j1,k1,i2,j2,k2 &
       & ,yinyang_t,yinyang,n_int &
       & ,y,z &
       & ,dy,dz
  use const_def, only: pi
  use yinyang_def, only: &
       &  yinyangcount_rcv,yinyangcount_snd &
       & ,yinyang_yz_rcv,yinyang_yz_snd &
       & ,yinyang_jk_snd,yinyang_jk_rcv &
       & ,yinyang_dyz,yzmask,yzmask_act,yo,zo,coss,sins
  
  use comm_def, only: myrnk_loca,myrank,npe,blck,ib,jb,kb,mpi_yinyang
  use implicit_def
  !$ use omp_lib
  implicit none

  include "mpif.h"

  real(KIND(0.d0)), dimension(:,:), allocatable :: yinyang_yz_snd_l,yinyang_yz_rcv_l
  !dir$ attributes align:64 :: yinyang_yz_snd_l,yinyang_yz_rcv_l
  integer, dimension(mpi_status_size) :: mstatus
  integer :: ireq_rcv0,ireq_snd0
  logical, dimension(nyg,nzg) :: flag_yz

  real(KIND(0.d0)) :: ycmin,ycmax,zcmin,zcmax

  real(KIND(0.d0)) :: sct,sco
  
  integer :: ibnd,id_neibor
  integer :: kk0,kk1,jj0,jj1
  integer :: ii,jj,kk
  integer :: npo
  integer :: in_out0,in_out1
  integer :: iyinyango
  integer :: jbo,kbo,mmx,np

  real(KIND(0.d0)) :: dyl,dzl
  real(KIND(0.d0)) :: flag0,flag0g,flag1,flag1g
  real(KIND(0.d0)) :: yt,zt

  integer :: jst,jen
  integer :: kst,ken
!----------------------------------------------------------------------|


  ! Yin-Yang convert
  do k = 1,nzg
  do j = 1,nyg
     yo(j,k) = acos(sin(y(j))*sin(z(k)))
     zo(j,k) = asin(cos(y(j))/sin(yo(j,k)))
     sct =  sin(y(j)   )*cos(z(k))
     sco = -sin(yo(j,k))*cos(zo(j,k))
     if(sign(1.d0,sct) /= sign(1.d0,sco)) then
        zo(j,k) = sign(1.d0,zo(j,k))*pi - zo(j,k)
     endif
     coss(j,k) = -sin(zo(j,k))*sin(z(k))
     sins(j,k) =  cos(zo(j,k))/sin(y(j))
  enddo
  enddo
  
  dyl = (ymax-ymin)/dble(jx0) ! The width of "grid" in y-direction
  dzl = (zmax-zmin)/dble(kx0) ! The width of "grid" in z-direction
  
  flag_yz(:,:) = .true.
  ycmin = ymin + dy(1+marginy)*n_int
  ycmax = ymax - dy(1+marginz)*n_int
  zcmin = zmin + dz(1+marginy)*n_int
  zcmax = zmax - dz(1+marginz)*n_int

  ! the number of receiving data
  yinyangcount_rcv(:) = 0
  yinyangcount_snd(:) = 0
  ! searching recieve data

  yinyang_yz_snd(:,:,:) = 1.d0
  yinyang_yz_rcv(:,:,:) = 0.d0

  yinyang_jk_snd(:,:,:) = -1
  yinyang_jk_rcv(:,:,:) = -1

  if(mpi_yinyang == 0) then ! Judge the Yin-Yang
  !if(1 == 1) then ! Judge the Yin-Yang
     yzmask_act(:,:) = 1.d0

     ! indicator for myrnk_loca
     ! iyinyango = 1: Yin grid
     ! iyinyango = 2: Yang grid
     iyinyango = 2
     do ibnd = 3,6
        id_neibor = blck%ineib(myrank,ibnd)
        ! locate the boundary grid
        kk0 = 1
        kk1 = nzg
        jj0 = 1
        jj1 = nyg
        
        if(ibnd == 3) then
           jj0 = nyg - margin + 1
           jj1 = nyg
        endif
        if(ibnd == 4) then
           jj0 = 1
           jj1 = margin
        endif
        
        if(ibnd == 5) then
           kk0 = nzg - margin + 1
           kk1 = nzg
        endif
        if(ibnd == 6) then
           kk0 = 1
           kk1 = margin
        endif
        
        if(id_neibor == -1) then
           do k = kk0,kk1
           do j = jj0,jj1
              ! locate my location in the other coordinate
              ! MPI rank is investigated
              jbo = int((yo(j,k) - ymin)/dyl)
              kbo = int((zo(j,k) - zmin)/dzl)
              if(jbo > jx0-1) write(*,*) "bad j"
              if(kbo > kx0-1) write(*,*) "bad k"
              npo = myrnk_loca(iyinyango,ib,jbo,kbo) ! specify the location in Yang grid
              yinyangcount_rcv(npo) = yinyangcount_rcv(npo) + 1

              yinyang_yz_rcv(1,yinyangcount_rcv(npo),npo) = yo(j,k) ! location in Yang
              yinyang_yz_rcv(2,yinyangcount_rcv(npo),npo) = zo(j,k) ! location in Yang

              yinyang_yz_rcv(3,yinyangcount_rcv(npo),npo) = y(j) ! location in Yin
              yinyang_yz_rcv(4,yinyangcount_rcv(npo),npo) = z(k) ! location in Yin
              
              ! warning ! rcv is correct
              yinyang_jk_rcv(1,yinyangcount_rcv(npo),npo) = j
              yinyang_jk_rcv(2,yinyangcount_rcv(npo),npo) = k
           enddo
           enddo
        endif
     enddo
  else ! for Yang-grid
     iyinyango = 1     
     ! in_out0 : judgement for previous grid
     ! in_out1 : judgement for current grid
     ! in_out* = 0: out of Yin-grid (precisely speaking out of criterion)
     ! in_out* = 1: in Yin-grid (precisely speaking in criterion)

     yzmask_act(:,:) = 1.d0
     do k = 1,nzg
     do j = 1,nyg
        if(    ycmin < yo(j,k) .and. ycmax > yo(j,k) .and. &
             & zcmin < zo(j,k) .and. zcmax > zo(j,k)) then
           yzmask_act(j,k) = 0.d0
        endif        
     enddo
     enddo
     ! scaninng y-direction
     !if(iyzcase /= 2) then
     !if(1 /= 1) then
     do k = 1,nzg
        do j = 2,nyg
           in_out0 = 0 ! 0 means outside Yin grid
           in_out1 = 0 ! 1 means inside  Yin grid
           if(    ycmin < yo(j  ,k) .and. ycmax > yo(j  ,k) .and. &
                & zcmin < zo(j  ,k) .and. zcmax > zo(j  ,k)) then
              in_out0 = 1
           endif
           
           if(    ycmin < yo(j-1,k) .and. ycmax > yo(j-1,k) .and. &
                & zcmin < zo(j-1,k) .and. zcmax > zo(j-1,k)) then
              in_out1 = 1
           endif

           
           
           if(in_out0 /= in_out1) then
              if(in_out0 == 1) then! current position is inside
                 ! upper boundary
                 jst = j
                 jen = min(j + marginy - 1,nyg)
                 if(jen /= j + marginy - 1) then
                    ! this is allowed only where the MPI boundary
                    if(jb == jx0 - 1) then
                       write(*,*) "Bad Yin-Yang yloc01"
                       call finish
                    endif
                 endif
              else ! previous position is inside
                 ! lower boundary
                 jst = max(j - marginy,1)
                 jen = j - 1
                 if(jst /= j - marginy) then
                    ! this is allowed only where the MPI boundary
                    if(jb == 0) then
                       write(*,*) "Bad Yin-Yang yloc02"
                       call finish
                    endif
                 endif
              endif
              do jj = jst,jen
                 flag_yz(jj,k) = .false.
                 jbo = int((yo(jj,k) - ymin)/dyl)
                 kbo = int((zo(jj,k) - zmin)/dzl)
                 npo = myrnk_loca(iyinyango,ib,jbo,kbo)
                 if(npo == -1) then
                    write(*,*) npo,iyinyango,ib,jbo,kbo
                 endif
                 ! specify the location in Yang grid
                 yinyangcount_rcv(npo) = yinyangcount_rcv(npo) + 1
                 
                 yinyang_yz_rcv(1,yinyangcount_rcv(npo),npo) = yo(jj,k)
                 yinyang_yz_rcv(2,yinyangcount_rcv(npo),npo) = zo(jj,k)
                 
                 yinyang_yz_rcv(3,yinyangcount_rcv(npo),npo) = y(jj)
                 yinyang_yz_rcv(4,yinyangcount_rcv(npo),npo) = z(k)
                 
                 ! warning ! rcv is corr2ect
                 yinyang_jk_rcv(1,yinyangcount_rcv(npo),npo) = jj
                 yinyang_jk_rcv(2,yinyangcount_rcv(npo),npo) = k
              enddo
           endif
        enddo
     enddo
     
     do j = 1,nyg
        do k = 2,nzg
           in_out0 = 0 ! 0 means outside Yin grid
           in_out1 = 0 ! 1 means inside  Yin grid
           if(    ycmin < yo(j,k  ) .and. ycmax > yo(j,k  ) .and. &
                & zcmin < zo(j,k  ) .and. zcmax > zo(j,k  )) then
              in_out0 = 1
           endif
           
           if(    ycmin < yo(j,k-1) .and. ycmax > yo(j,k-1) .and. &
                & zcmin < zo(j,k-1) .and. zcmax > zo(j,k-1)) then
              in_out1 = 1
           endif
           
           if(in_out0 /= in_out1) then
              if(in_out0 == 1) then! current position is inside
                 ! upper boundary
                 kst = k
                 ken = min(k + marginz - 1,nzg)
                 if(ken /= k + marginz - 1) then
                    ! this is allowed only where the MPI boundary
                    if(kb == kx0 - 1) then
                       write(*,*) "Bad Yin-Yang yloc01"
                       call finish
                    endif
                 endif
              else ! previous position is inside
                 ! lower boundary
                 kst = max(k - marginz,1)
                 ken = k - 1
                 if(kst /= k - marginz) then
                    ! this is allowed only where the MPI boundary
                    if(kb == 0) then
                       write(*,*) "Bad Yin-Yang yloc02"
                       call finish
                    endif
                 endif
              endif
        
              do kk = kst,ken
                 if(flag_yz(j,kk)) then
                    if(    ycmin < yo(j,kk) .and. ycmax > yo(j,kk) .and. &
                         & zcmin < zo(j,kk) .and. zcmax > zo(j,kk)) then
                    
                       jbo = int((yo(j,kk) - ymin)/dyl)
                       kbo = int((zo(j,kk) - zmin)/dzl)
                       npo = myrnk_loca(iyinyango,ib,jbo,kbo)
                       ! specify the location in Yang grid
                       yinyangcount_rcv(npo) = yinyangcount_rcv(npo) + 1
                       
                       yinyang_yz_rcv(1,yinyangcount_rcv(npo),npo) = yo(j,kk)
                       yinyang_yz_rcv(2,yinyangcount_rcv(npo),npo) = zo(j,kk)
                       
                       yinyang_yz_rcv(3,yinyangcount_rcv(npo),npo) = y(j)
                       yinyang_yz_rcv(4,yinyangcount_rcv(npo),npo) = z(kk)
                       
                       ! warning ! rcv is correct
                       yinyang_jk_rcv(1,yinyangcount_rcv(npo),npo) = j
                       yinyang_jk_rcv(2,yinyangcount_rcv(npo),npo) = kk
                    endif
                 endif
              enddo
           endif
        enddo
     enddo
  endif
  
  call mpi_alltoall(yinyangcount_rcv(0),1,mpi_integer,yinyangcount_snd(0) &
       & ,1,mpi_integer,mpi_comm_world,merr)
  
  ! sending "information"
  do np = 0,npe-1
     if(yinyangcount_rcv(np) /= 0) then
        if(allocated(yinyang_yz_rcv_l)) deallocate(yinyang_yz_rcv_l)
        allocate(yinyang_yz_rcv_l(4,yinyangcount_rcv(np)))
        yinyang_yz_rcv_l(1:4,1:yinyangcount_rcv(np)) = yinyang_yz_rcv(1:4,1:yinyangcount_rcv(np),np)
        mmx = yinyangcount_rcv(np)*4
        call mpi_isend(yinyang_yz_rcv_l,mmx,mpi_real8,np,10,mpi_comm_world,ireq_rcv0,merr)
     endif

  ! sending "information"
     if(yinyangcount_snd(np) /= 0) then
        mmx = yinyangcount_snd(np)
        if(allocated(yinyang_yz_snd_l)) deallocate(yinyang_yz_snd_l)!
        allocate(yinyang_yz_snd_l(4,yinyangcount_snd(np)))
        mmx = yinyangcount_snd(np)*4
        call mpi_irecv(yinyang_yz_snd_l,mmx,mpi_real8,np,10,mpi_comm_world,ireq_snd0,merr)
        call mpi_waitall(1,ireq_snd0,mstatus,merr)
        yinyang_yz_snd(1:4,1:yinyangcount_snd(np),np) = yinyang_yz_snd_l(1:4,1:yinyangcount_snd(np))
     endif
     
     if(yinyangcount_rcv(np) /= 0) call mpi_waitall(1,ireq_rcv0,mstatus,merr)
  enddo

  
  flag0 = 0.d0
  flag1 = 0.d0
  do np = 0,npe-1
     if(yinyangcount_snd(np) /=0) then
        do ii = 1,yinyangcount_snd(np)
           yt = yinyang_yz_snd(1,ii,np)
           zt = yinyang_yz_snd(2,ii,np)

           j = int((yt - y(1))/dy(1+marginy)) + 1
           k = int((zt - z(1))/dz(1+marginz)) + 1
                     
           if(j == 1) write(*,*) yt,y(1),yt-y(1),ymin,dy
           yinyang_jk_snd(1,ii,np) = j
           yinyang_jk_snd(2,ii,np) = k


           if(j+1 > nyg .or. j < 1 .or. k+1 > nzg .or. k < 1 ) then
              write(*,*) yo*180./pi,zo*180./pi,jj,kk
              flag0 = 1.d0
           endif

           if(y(j+1) > ymax .or. y(j) < ymin .or. z(k+1) > zmax .or. z(k) < zmin ) then
              write(*,*) yo*180./pi,zo*180./pi,jj,kk,myrank,np
              flag1 = 1.d0
           endif

goto 1000           
           ! Lagrange interpolation coefficient 3rd order (j-1, j, j+1, j+2)
           ! for j - 1
           yinyang_dyz(1 ,ii,np) = &
                &  (yt    -y(j  ))*(yt    -y(j+1))*(yt    -y(j+2))/ &
                & ((y(j-1)-y(j  ))*(y(j-1)-y(j+1))*(y(j-1)-y(j+2)))

           ! for j
           yinyang_dyz(2 ,ii,np) = &
                &  (yt    -y(j-1))*(yt    -y(j+1))*(yt    -y(j+2))/ &
                & ((y(j  )-y(j-1))*(y(j  )-y(j+1))*(y(j  )-y(j+2)))

           ! for j + 1
           yinyang_dyz(3 ,ii,np) = &
                &  (yt    -y(j-1))*(yt    -y(j  ))*(yt    -y(j+2))/ &
                & ((y(j+1)-y(j-1))*(y(j+1)-y(j  ))*(y(j+1)-y(j+2)))

           ! for j + 2
           yinyang_dyz(4 ,ii,np) = &
                &  (yt    -y(j-1))*(yt    -y(j  ))*(yt    -y(j+1))/ &
                & ((y(j+2)-y(j-1))*(y(j+2)-y(j  ))*(y(j+2)-y(j+1)))
           
           ! for k - 1
           yinyang_dyz(5 ,ii,np) = &
                &  (zt    -z(k  ))*(zt    -z(k+1))*(zt    -z(k+2))/ &
                & ((z(k-1)-z(k  ))*(z(k-1)-z(k+1))*(z(k-1)-z(k+2)))

           ! for k
           yinyang_dyz(6 ,ii,np) = &
                &  (zt    -z(k-1))*(zt    -z(k+1))*(zt    -z(k+2))/ &
                & ((z(k  )-z(k-1))*(z(k  )-z(k+1))*(z(k  )-z(k+2)))

           ! for k + 1
           yinyang_dyz(7 ,ii,np) = &
                &  (zt    -z(k-1))*(zt    -z(k  ))*(zt    -z(k+2))/ &
                & ((z(k+1)-z(k-1))*(z(k+1)-z(k  ))*(z(k+1)-z(k+2)))

           ! for k + 2
           yinyang_dyz(8 ,ii,np) = &
                &  (zt    -z(k-1))*(zt    -z(k  ))*(zt    -z(k+1))/ &
                & ((z(k+2)-z(k-1))*(z(k+2)-z(k  ))*(z(k+2)-z(k+1)))
1000 continue
           
           do n = 1,2*margin
              yinyang_dyz(n,ii,np) = 1.d0
              do l = 1,2*margin
                 if (l /= n) then
                    yinyang_dyz(n,ii,np) = &
                         & yinyang_dyz(n,ii,np) &
                         & *(yt - y(j-margin+l))/(y(j-margin+n) - y(j-margin+l))
                 endif
              enddo
           enddo

           do n = 1,2*margin
              yinyang_dyz(n+2*margin,ii,np) = 1.d0
              do l = 1,2*margin
                 if (l /= n) then
                    yinyang_dyz(n+2*margin,ii,np) = &
                         & yinyang_dyz(n+2*margin,ii,np) &
                         & *(zt - z(k-margin+l))/(z(k-margin+n) - z(k-margin+l))
                 endif
              enddo
           enddo           
        enddo
     endif
  enddo
  
  call mpi_allreduce(flag0,flag0g,1,mpi_real8,mpi_max,mpi_comm_world,merr)

  
  if(flag0g == 1.d0) then
     if(myrank == 0) then
        write(*,*) ""
        write(*,*) "### Coennection is bad between Yin-Yang.                   ###"
        write(*,*) "### Please check the subroutine yinyang_boundary_intro.f90 ###"
        write(*,*) "### Finishing..."
        write(*,*) ""
     endif             
  endif

  call mpi_allreduce(flag1,flag1g,1,mpi_real8,mpi_max,mpi_comm_world,merr)
  if(flag1g == 1.d0) then
     if(myrank == 0) then
        write(*,*) ""
        write(*,*) "### Yin-Yang cannot covers each other                      ###"
        write(*,*) "### Please check the subroutine yinyang_boundary_intro.f90 ###"
        write(*,*) "### Finishing..."
        write(*,*) ""
     endif             
  endif

  if(flag0g == 1.d0 .or.flag1g == 1.d0) then
     call mpi_finalize(merr)
     stop     
  endif

  yzmask = 1.d0
  
  if(mpi_yinyang == 1) then
     do k = 1,nzg
     do j = 1,nyg
!        yo = acos(sin(y(j))*sin(z(k)))
!        zo = asin(cos(y(j))/sin(yo))
!        sct =  sin(y(j))*cos(z(k))
!        sco = -sin(yo  )*cos(zo  )
!        if(sign(1.d0,sct) /= sign(1.d0,sco)) then
!           zo = sign(1.d0,zo)*pi - zo
!        endif
        if(    ymin <= yo(j,k) .and. ymax >= yo(j,k) .and. &
             & zmin <= zo(j,k) .and. zmax >= zo(j,k)) then
           yzmask(j,k) = 0.d0
        endif
     enddo
     enddo
  endif
  
  ! YinYang
#endif  
  return
end subroutine yinyang_init
