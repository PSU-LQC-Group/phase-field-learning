MODULE VARIABLES
IMPLICIT NONE

!GLOBAL VARIABLES
  integer               :: i,j, k, Nx, Ny , Rad, t_loop, il
  integer               :: file_skip, tmax, restart, restart_time
  integer               :: B0l_stages, B0l_jump, B0l_change_time, B0l_quench_time
  integer, parameter    :: Nmax=1200, idum=23020230
  real*8,  parameter    :: pi = 3.14159265, theta=pi/5
  real*8                :: dt, dx, tt,uu,vv,ww,Box,B2l,B0l,RKK,rno,co,eta,Mn,Mc
  real*8                :: t4,mid_sys1,mid_sys2,qo,cs,summ,DB0l,B0l_range,xr,yr,d1,d2
  real*8                :: F(0:Nmax,0:Nmax)
  real*8                :: G(0:Nmax,0:Nmax), Bl(0:Nmax,0:Nmax)
  real*8                :: grad2_n(0:Nmax,0:Nmax), grad2_c(0:Nmax,0:Nmax)
  real*8                :: grad4_n(0:Nmax,0:NMax), grad2_nc(0:Nmax,0:Nmax)
  real*8                :: grad4_nc(0:Nmax,0:Nmax), RHS(0:Nmax,0:Nmax)
  real*8                :: grad2_mn(0:Nmax,0:Nmax), grad2_mc(0:Nmax,0:Nmax)
  character(len=9)      :: cn

CONTAINS

subroutine read_globals
  

  !read in initial parameters from file (see inut for definitions)
  open(1,file='input')
    read(1,*) dx
    read(1,*) dt
    read(1,*) tt
    read(1,*) uu
    read(1,*) vv
    read(1,*) ww
    read(1,*) Box
    read(1,*) B2l
    read(1,*) B0l
    read(1,*) B0l_range
    read(1,*) B0l_stages
    read(1,*) B0l_quench_time
    read(1,*) RKK
    read(1,*) rno
    read(1,*) co
    read(1,*) cs
    read(1,*) eta
    read(1,*) file_skip
    read(1,*) tmax
    read(1,*) Nx
    read(1,*) Ny
    read(1,*) Rad
	read(1,*) restart
	read(1,*) restart_time
  close(1)

  !initialize specific parameters 
  F=0.0d0; G=0.0d0

  !set movilities to be the same and=1 (like rescaling time by mobility)
  Mn=1.0d0 ;Mc=1.0d0

  !define inverse square mesh size
  t4 = 1 / dx**2

  !test that system sizes do not exceed declared array dimensions (avoid this using dynamic alloc)
  if(Nx.gt.Nmax.or.Ny.gt.Nmax)then 
   print*, 'One of the system dimensions exceeds array dimensions'
   stop
  endif

  DB0l=B0l_range/(B0l_stages-1)
  B0l_change_time=B0l_quench_time/B0l_stages
 
  
 print*,'Global data:'
 print*,'_________________________________'
 print 10, B0l-Box, DB0l, B0l_change_time, ww, dt, dx, tmax, file_skip, Nx, Ny, rno, co, cs, Rad, theta
 10 FORMAT(' DB0(init) =',F8.4,/ ' DB0 change=', F8.4, / ' DB0 change rate (steps)=', I6, / &
           & ' ww=',F8.4, / ' dt =',F6.4,/' dx = ',F6.4,/ ' tmax=',I8,/ ' file_skip=',I6,/ ' Nx=',I6,/ &
           & ' Ny=',I6,/ ' rno =',F5.3,/ ' co =',F6.3,/ ' cs =' ,F8.4, / ' seed diameter=',I6, / ' misorientation=', F8.4)
           
 print*,'_________________________________'

end subroutine read_globals


END MODULE VARIABLES
