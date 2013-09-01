MODULE VARIABLES
IMPLICIT NONE

!GLOBAL VARIABLES
  integer               :: i,j, k, Nx, Ny , t_loop, il
  integer               :: file_skip, tmax, Rad
  integer, parameter    :: Nxmax=1200, Nymax=1200, idum=92923932

  real*8                :: dt, dx, dx2_in, dx_in, dt_in, ave_psi, a_12
  real*8                :: MAG_sq2, a_s, epsilon, epsilon_4,  partition
  real*8                :: OMEGA, dist, derx, dery, Atheta, Aptheta
  real*8                :: cipjp, cipjm, cimjm, cimjp, JR, JL, JT, JB
  real*8                :: MAG, JR_a, JL_a, JT_a, JB_a, x_cross
  real*8                :: lambda, W2, D_bar, tau, x1, x2, y1, y2
  real*8,  parameter    :: eps=1.0d-8, a_1=0.8839, a_2=0.6267, a_t=1.0/(2*sqrt(2.0d0))

  real*8                :: PSI(0:Nxmax,0:Nymax), C(0:Nxmax,0:Nymax)
  real*8                :: grad2(0:Nxmax,0:Nymax), EU(0:Nxmax,0:NyMax) 
  real*8                :: DPSI(0:Nxmax,0:Nymax)


  character(len=9)      :: cn

CONTAINS

subroutine read_globals
  
 !read in initial parameters from file (see input file for definitions)
  open(1,file='input')
    read(1,*) dx
    read(1,*) dt
    read(1,*) file_skip
    read(1,*) tmax
    read(1,*) Nx
    read(1,*) Ny
    read(1,*) lambda
    read(1,*) OMEGA
    read(1,*) epsilon_4
    read(1,*) Rad
    read(1,*) partition
  close(1)

  !initialize specific parameters
  PSI=0.0d0

  C=0.0d0

  EU=0.0d0

  !test that system sizes do not exceed declared array dimensions (avoid this using dynamic alloc)
  if(Nx+1.gt.Nxmax.or.Ny+1.gt.Nymax)then 
   print*, 'One of the system dimensions exceeds array dimensions'
   stop
  endif
	
  dx2_in = 1.0d0 / dx**2
  dx_in = 1.0d0 / dx
  dt_in = 1.0d0 / dt
  D_bar=a_2*lambda
  a_s = 1.0d0-3*epsilon_4
  epsilon = 4.0d0*epsilon_4/a_s 
  a_12 = 4.0d0*a_s*epsilon

 print*,'Global data:'
 print*,'_________________________________'
 print 10, OMEGA, epsilon_4, lambda, dt, dx, tmax, file_skip, Nx, Ny
 10 FORMAT(' DELTA=',F8.4,/ ' epsilon_4=',F8.4,/ ' lambda=',F8.4,/ ' dt =',F6.4,/' dx = ',F6.4,/ & 
       &    ' tmax=',I6,/ ' file_skip=',I6,/ ' Nx=',I6,/ ' Ny=',I6)
 print*,'_________________________________'

end subroutine read_globals


END MODULE VARIABLES
