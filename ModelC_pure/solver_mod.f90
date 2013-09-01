!!
!!Numerical integration of Model C for a pure material from text by Provatas and Elder
!
! using central difference (spherical laplacian) formula used for nabla^2
!
MODULE SOLVER
USE VARIABLES
USE UTIL
IMPLICIT NONE

CONTAINS !-----------

 subroutine calculate

  !initialize random seed for random number generation (if needed)
  call srand(idum)

  !!!!!!set initial conditions. NOTE: PSI-->order parameter !!!!!!
  do i=1,Nx
   do j=1,Ny
     dist=sqrt( ( (i-1)*dx )**2 + ( (j-1)*dx )**2 )-Rad*dx
     PSI(i,j)=-tanh(dist/sqrt(2.0))
   end do  
  end do  
  U=DELTA

  !Enfrce zero flux  boundary conditions for next iteration
  call ZERO_FLUX(PSI)
  call ZERO_FLUX(U)

  !print intial PSI/U fields
  open(1,file='out_0',status='unknown')
    call print_2Dfields(PSI, U, Nx, Ny, 1)
  close(1)

  !open file to store interface position
  open(12,file='Int_position',status='unknown')

 !!!START TIME MARCHING!!!!!!
 do t_loop=1, tmax

    !!!!!!!!!!!!!!! order parameter (non-conserved) equation update !!!!!
    DO i=1,Nx
     DO j=1,Ny
      !!compute corner [ghost] node values around (i,j) (see finite volume text in appendinx)
        cipjp=( PSI(i+1,j+1) + PSI(i,j+1) + PSI(i,j) + PSI(i+1,j) )/4.0d0
        cipjm=( PSI(i+1,j) + PSI(i,j) + PSI(i,j-1) + PSI(i+1,j-1) )/4.0d0
        cimjp=( PSI(i,j+1) + PSI(i-1,j+1) + PSI(i-1,j) + PSI(i,j) )/4.0d0
        cimjm=( PSI(i,j) + PSI(i-1,j) + PSI(i-1,j-1) + PSI(i,j-1) )/4.0d0

      !! Calculate right edge flux
        DERX= PSI(i+1,j)-PSI(i,j) 
        DERY= cipjp - cipjm 
        call set_aniso_parameters  !!Computes local A & A' in text
        JR = Atheta * ( Atheta*DERX + Aptheta*DERY )
       
      !! Calculate left edge flux
        DERX= PSI(i,j)-PSI(i-1,j) 
        DERY= cimjp - cimjm
        call set_aniso_parameters  !!Computes local A & A' in text
        JL = Atheta * ( Atheta*DERX + Aptheta*DERY )

      !! Calculate top edge flux
        DERX= cipjp - cimjp
        DERY= PSI(i,j+1)-PSI(i,j) 
        call set_aniso_parameters  !!Computes local A & A' in text
        JT = Atheta * ( Atheta*DERY - Aptheta*DERX )

      !! Calculate bottom edge flux
        DERX= cipjm - cimjm
        DERY= PSI(i,j)-PSI(i,j-1) 
        call set_aniso_parameters  !!Computes local A & A' in text
        JB = Atheta * ( Atheta*DERY - Aptheta*DERX )

      !! Calculate local a(n) term for tau
        DERX= PSI(i+1,j)-PSI(i-1,j) 
        DERY= PSI(i,j+1)-PSI(i,j-1) 
        call set_aniso_parameters  !!Computes local A & A' in text

      !!Update PSI-equation (see text)
        dPSI(i,j) = (dt / Atheta**2) * ( dx2_in * ( (JR - JL) + (JT -JB ) ) &
           &     +  PSI(i,j)-PSI(i,j)**3  &
           &     - lambda*U(i,j)*( 1.0d0 - PSI(i,j)**2 )**2 )

        PSI(i,j)=PSI(i,j) + dPSI(i,j)

     END DO  
    END DO  
    

    !!!!!!!!!!!!!!! diffusion (conserved) equation update !!!!!!!!!!!!!
    call NABLA2(U,grad2)
    DO i=1,Nx
     DO j=1,Ny
      U(i,j) = U(i,j) + alpha *dt*grad2(i,j) + 0.5d0*dPSI(i,j) 
     END DO
    END DO

  
    !Enfrce zero flux  boundary conditions for next iteration
    call ZERO_FLUX(PSI)
    call ZERO_FLUX(U)


    !!!!!!!!!!!!!!!!locate and print interface  position !!!!!!!!!!
    DO i=1,Nx
     !check for phi-crossing y-axis 
     IF( PSI(i,1) * PSI(i+1,1) .lt. 0.0d0 )THEN
       !from (x,PSI(x)) points below and above PSI=0 axis, construct a line
       y1=PSI(i,1)
       x1=(i-1)*dx
       y2=PSI(i+1,1)
       x2=i*dx
       exit
     ENDIF
    END DO
   !interface defiend as where line crosses PSI=0
    x_cross=x1 - ( ( x2-x1 ) / ( y2-y1) ) * y1
    write(12,*)  (t_loop-1)*dt,x_cross


    !!!!!!!!!!!!!!!!!!! I/O  (print out PSI)!!!!!!!!!!!!!!!!!!!!!!!!
    if(mod(t_loop,file_skip)+1 == 1)then
       call chari(t_loop,cn,il)
       open(1,file='out_'//cn(1:il),status='unknown')
         call print_2Dfields( PSI, U, Nx, Ny, 1)
       close(1)
    end if


 end do !time loop
 
  !close  file that stores interface position
  close(12)
 end subroutine calculate 

  !-------------------------------------------------------------------
  !for given derivatives in the x (DERX) and y (DERY) directions, computes anisotropy (see text)
  subroutine set_aniso_parameters

  MAG_sq2 = (DERX**2 + DERY**2)**2
  IF(MAG_sq2.gt.eps)THEN

   Atheta = a_s*( 1 + epsilon*(DERX**4 + DERY**4) / MAG_sq2   )
   Aptheta = a_12*DERX*DERY*(DERX**2 - DERY**2) /  MAG_sq2

  ELSE 
  
     Atheta = a_s
     Aptheta = 0.0d0

  ENDIF

  end subroutine set_aniso_parameters

  !-------------------------------------------------------------------
  !finite difference formulation of laplacian (see text)
  subroutine NABLA2(A,GR2)
  real*8, dimension(0:,0:) ::A,GR2
   GR2=0.0d0
   do i=1,Nx
    do j=1,Ny
       GR2(i,j) = ( (A(i+1,j) + A(i-1,j) + A(i,j+1) +A(i,j-1))/2.0  &
              & +   (A(i+1,j+1)+A(i-1,j+1)+A(i+1,j-1) +A(i-1,j-1))/4.0 &
              & -   3*A(i,j) )* dx2_in
 
    end do
   end do
  end subroutine NABLA2

 !-----------------------
 !eforces periodic boundary conditions on array A (see text)
  subroutine PERIODIC(A)
  real*8, dimension(0:,0:) ::A
    A(Nx+1,:)    = A(1,:)
    A(0,:)      = A(Nx,:)
    A(:,Ny+1)    = A(:,1)
    A(:,0)      = A(:,Ny)
  end subroutine periodic

 !-----------------------
 !enforces zero flux boundary conditions on array A (see text)
  subroutine ZERO_FLUX(A)
  real*8, dimension(0:,0:) ::A
    A(Nx+1,:)    = A(Nx-1,:)
    A(0,:)      = A(2,:)
    A(:,Ny+1)    = A(:,Ny-1)
    A(:,0)      = A(:,2)
  end subroutine ZERO_flux
  !!!!!!!!!!!!!!!!!!!!!!!




END MODULE SOLVER
