MODULE VARIABLES
IMPLICIT NONE

!GLOBAL VARIABLES
  integer               :: i,j, k, Nx, Ny , t_loop, il
  integer               :: file_skip, tmax
  integer, parameter    :: Nmax=2200,idum=92923932
  real*8                :: dt, dx, a_2, a_4, M, W2, dx2_in, ave_psi
  real*8                :: PSI(0:Nmax,0:Nmax)
  real*8                :: grad2(0:Nmax,0:Nmax)
  real*8                :: RHS(0:Nmax,0:NMax)
  character(len=9)      :: cn

CONTAINS

 subroutine read_globals
  
 !read initial parameters from file 
  open(1,file='input')
    read(1,*) dx          !space step
    read(1,*) dt          !time step
    read(1,*) a_2         !coefficient of the phi^2 term in f(phi)
    read(1,*) a_4         !coefficient of the phi^4 term in f(phi)
    read(1,*) W2          !gradient squared coefficient squared W^2
    read(1,*) file_skip   ! how often to print
    read(1,*) tmax        !how many total time steps to take
    read(1,*) Nx          !system size in the x direction
    read(1,*) Ny          !system size in the y direction 
  close(1)

  !initialize specific parameters 
  PSI=0.0d0

  M=1.0d0

  dx2_in = 1 / dx**2

  !test that system sizes do not exceed declared dimensions (avoid this using dynamic allocation)
  if(Nx+1.gt.Nmax.or.Ny+1.gt.Nmax)then 
   print*, 'One of the system dimensions exceeds array dimensions'
   stop
  endif
  
 print*,'Global data:'
 print*,'_________________________________'
 print 10, a_2, dt, dx, tmax, file_skip, Nx, Ny
 10 FORMAT(' a_2 =',F8.4,/ 'dt =',F6.4,/' dx = ',F6.4,/ ' & 
       &     tmax=',I6,/ ' file_skip=',I6,/ ' Nx=',I6,/ ' Ny=',I6)
 print*,'_________________________________'

end subroutine read_globals


END MODULE VARIABLES
