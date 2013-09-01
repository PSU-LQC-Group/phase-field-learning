!!
!!Module contains routines for numerical integration of Model B from text by Provatas and Elder
! using central difference (spherical laplacian) formula used for nabla^2
!
MODULE SOLVER 
USE VARIABLES !!where variables are defined
USE UTIL      !!For printig and random number generation
IMPLICIT NONE

CONTAINS !-----------

 !Main routine for time evolution of model B
 subroutine calculate


  !set initial conditions. NOTE: PSI-->order parameter 
  !initialize random seed for random number generation in gasdev()
  call srand(idum)
   do i=1,Nx
    do j=1,Ny
	  !set PSI=0 + fluctuations 
     PSI(i,j)=0.000001*gasdev2()
   end do  
  end do  

  !initialize to zero the avergae of the ordr parameter
  ave_psi=0.0d0
  do i=1,Nx
   do j=1,Ny
     ave_psi=ave_psi+PSI(i,j)
   end do
  end do 
    
  !write the initial average value of the order parameter
  open(2,file='averga_psi',status='unknown')
  write(2,*) ave_psi/((Nx)*(Ny))

  !print intial PSI field
  open(1,file='out_0',status='unknown')
    call print_2Dfield(PSI, Nx, Ny, 1)
  close(1)

 !!!START TIME MARCHING
 do t_loop=1, tmax

    !!!!!!!!!!!!!!! order parameter equation update (see text) !!!!!!!!!!!!!!!!!!!!

    !1. periodic BC for PSI and compute grad^2(PSI) array
    call PERIODIC(PSI)
    call NABLA2(PSI,grad2)

    !2. calculate right hand side of model A
    do i=1, Nx
     do j=1,Ny
      RHS(i,j)=-W2*grad2(i,j)+a_2*PSI(i,j)+a_4*PSI(i,j)**3
     end do 
    end do

    !3. re-enforce periodic BC for PSI
    call PERIODIC(RHS)
    call NABLA2(RHS,grad2mu)

    !4. take one time step forward
    do i=1, Nx
     do j=1,Ny
       PSI(i,j)=PSI(i,j)+dt*M*grad2mu(i,j)
     end do
    end do

    !initialize to zero average for this time step
    ave_psi=0.0d0
    do i=1, Nx
     do j=1,Ny
       ave_psi=ave_psi+PSI(i,j)
     end do 
   end do 

    !!!!!!!!!!!!!!!!!!! I/O  (print out average & PSI)!!!!!!!!!!!!!!!!!!!!!!!!
    write(2,*) ave_psi/((Nx)*(Ny))
    if(mod(t_loop,file_skip)+1 == 1)then
       call chari(t_loop,cn,il)
       open(1,file='out_'//cn(1:il),status='unknown')
         call print_2Dfield( PSI, Nx, Ny, 1)
       close(1)
    end if

 end do !time loop

 !close file opened to store the average of the order parameter
 close(2)
 end subroutine calculate 


  !-------------------------------------------------------------------
  !finite difference calculation of Laplacian of A, stores answer in GR2. (See appendix)
   subroutine NABLA2(A,GR2)
   real*8, dimension(0:,0:) ::A,GR2
    GR2=0.0d0
    do i=1,Nx
     do j=1,Ny
        GR2(i,j) = ( (A(i+1,j) + A(i-1,j) + A(i,j+1) +A(i,j-1))/2.0d0  &
              & +   (A(i+1,j+1)+A(i-1,j+1)+A(i+1,j-1) +A(i-1,j-1))/4.0d0 &
              & -   3.0d0*A(i,j) )*dx2_in
     end do
    end do
   end subroutine NABLA2

 !-----------------------
 !enforceis periodic boundary conditions (see text)
  subroutine PERIODIC(A)
  real*8, dimension(0:,0:) ::A
    A(Nx+1,:)    = A(1,:)
    A(0,:)      = A(Nx,:)
    A(:,Ny+1)    = A(:,1)
    A(:,0)      = A(:,Ny)
  end subroutine periodic
  !!!!!!!!!!!!!!!!!!!!!!!

END MODULE SOLVER
