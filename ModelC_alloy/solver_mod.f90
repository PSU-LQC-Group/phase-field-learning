!!
!!Numerical integration of Model C for a dilute alloy from text by Provatas and Elder
!
! using central difference (spherical laplacian) formula used for nabla^2
!
MODULE SOLVER
USE VARIABLES
USE UTIL
IMPLICIT NONE

CONTAINS !-----------

 subroutine calculate

 !initialize random number seed for random number genertion (if needed)
 call srand(idum)

 !set initial conditions. NOTE: PSI-->order parameter 
  do i=1,Nx
   do j=1,Ny
     dist=sqrt( ( (i-1)*dx )**2 + ( (j-1)*dx )**2 )-Rad*dx
     PSI(i,j)=-tanh(dist/sqrt(2.0))
   end do  
  end do  
  !initialize the reduced chemical potential u and concentration C
  EU=1.0d0 - (1-partition) * OMEGA
  C=0.5d0 * ( (1+partition) - (1-partition)*PSI )*EU

  !Enfrce zero flux  boundary conditions for next iteration
  call ZERO_FLUX(PSI)
  call ZERO_FLUX(C)
  call ZERO_FLUX(EU)

  !print intial PSI/C fields
  open(1,file='out_0',status='unknown')
    call print_2Dfields(PSI, C, Nx, Ny, 1)
  close(1)

  !open file to store interface position
  open(12,file='Int_position',status='unknown')

 !!!START TIME MARCHING
 do t_loop=1, tmax


    !!!!!!!!!!!!!!! order parameter (non-conserved) equation update (see text)!!!!!

    DO i=1,Nx
     DO j=1,Ny
      !!set corner [ghost] node values around (i,j) (see finite volume text in appendix)
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

      !!Update PSI-equation. U of pure model replaced by reduced chem pot EU (see text)
        dPSI(i,j) = (dt / Atheta**2) * ( dx2_in * ( (JR - JL) + (JT -JB ) ) &
           &     +  PSI(i,j)-PSI(i,j)**3 - ( lambda/(1.0d0-partition) ) & 
           &     * (EU(i,j) - 1.0d0)*( 1.0d0 - PSI(i,j)**2 )**2 )

        PSI(i,j)=PSI(i,j) + dPSI(i,j)

     END DO  
    END DO  

    !re-enforc zero flux on phi field, which will be needed again this iteration
    call ZERO_FLUX(PSI)

    !!!!!!!!!!!!!!! diffusion (conserved) equation update (see text) !!!!!!!!!!!!!
    DO i=1,Nx
     DO j=1,Ny
      !!set corner [ghost] node values around (i,j) (see finite volume text in appendix)
        cipjp=( PSI(i+1,j+1) + PSI(i,j+1) + PSI(i,j) + PSI(i+1,j) )/4.0d0
        cipjm=( PSI(i+1,j) + PSI(i,j) + PSI(i,j-1) + PSI(i+1,j-1) )/4.0d0
        cimjp=( PSI(i,j+1) + PSI(i-1,j+1) + PSI(i-1,j) + PSI(i,j) )/4.0d0
        cimjm=( PSI(i,j) + PSI(i-1,j) + PSI(i-1,j-1) + PSI(i,j-1) )/4.0d0

      !! Calculate right edge fluxe for anti-trapping term 
        DERX= PSI(i+1,j)-PSI(i,j) 
        DERY= cipjp - cipjm 
        MAG= DERX**2 + DERY**2
        JR= - D_bar*( C(i+1,j) + C(i,j) )*Q( 0.5d0*( PSI(i+1,j)+PSI(i,j) ) ) &
          & * ( EU(i+1,j) - EU(i,j) ) / ( EU(i+1,j) + EU(i,j) )  
        if(MAG.gt.eps)then
         JR_a= - (a_t/4.0)*(1.0d0-partition)*( EU(i+1,j) + EU(i,j)  ) &
          &  * ( DPSI(i+1,j) + DPSI(i,j) ) * DERX/sqrt( MAG )
        else
          JR_a=0.0d0
        endif

      !! Calculate left edge flux for anti-trapping term 
        DERX= PSI(i,j)-PSI(i-1,j) 
        DERY= cimjp - cimjm
        MAG= DERX**2 + DERY**2
        JL= - D_bar*( C(i-1,j) + C(i,j) )*Q( 0.5d0*( PSI(i-1,j)+PSI(i,j) ) ) &
          & * ( EU(i,j) - EU(i-1,j) ) / ( EU(i-1,j) + EU(i,j) ) 
        if(MAG.gt.eps)then
         JL_a= - (a_t/4.0)*(1.0d0-partition)*( EU(i-1,j) + EU(i,j)  ) &
           & * ( DPSI(i-1,j) + DPSI(i,j) ) * DERX/sqrt( MAG )
        else
         JL_a=0.0d0
        endif

      !! Calculate top edge flux for anti-trapping term 
        DERX= cipjp - cimjp
        DERY= PSI(i,j+1)-PSI(i,j) 
        MAG= DERX**2 + DERY**2
        JT= - D_bar*( C(i,j+1) + C(i,j) )*Q( 0.5d0*( PSI(i,j+1)+PSI(i,j) ) ) &
          & * ( EU(i,j+1) - EU(i,j) ) / ( EU(i,j+1) + EU(i,j) )
        if(MAG.gt.eps)then
         JT_a= - (a_t/4.0)*(1.0d0-partition)*( EU(i,j+1) + EU(i,j)  ) &
           & * ( DPSI(i,j+1) + DPSI(i,j) ) * DERY/sqrt( MAG )
        else
         JT_a=0.0d0
        endif

      !! Calculate bottom edge flux for anti-trapping term 
        DERX= cipjm - cimjm
        DERY= PSI(i,j)-PSI(i,j-1) 
        MAG= DERX**2 + DERY**2
        JB= - D_bar*( C(i,j-1) + C(i,j) )*Q( 0.5d0*( PSI(i,j-1)+PSI(i,j) ) ) &
          & * ( EU(i,j) - EU(i,j-1) ) / ( EU(i,j-1) + EU(i,j) ) 
        if(MAG.gt.eps)then
         JB_a= - (a_t/4.0)*(1.0d0-partition)*( EU(i,j-1) + EU(i,j)  ) &
           & * ( DPSI(i,j-1) + DPSI(i,j) ) * DERY/sqrt( MAG )
        else
         JB_a=0.0d0
        endif  

      !!Update C-equation (uses finite volume)
        C(i,j) = C(i,j) - dt*dx2_in*( (JR - JL) + (JT -JB ) ) &
               &        - dx_in*( (JR_a - JL_a) + (JT_a -JB_a ) ) 
     END DO
    END DO

  
    !!!!!!!!!!!!!!! update reduced chem pot, EU field, for next itertion!!!!!!!!!!!!!!
    DO i=1,Nx
     DO j=1,Ny
       EU(i,j)=2*C(i,j) / ( (1.0d0+partition) - (1.0d0-partition)*PSI(i,j) )
     END DO 
   END DO


    !Enfrce zero flux  boundary conditions for next iteration
    call ZERO_FLUX(PSI)
    call ZERO_FLUX(C)
    call ZERO_FLUX(EU)



    !!!!!!!!!!!!!!!!locate and print interface  position !!!!!!!!!!
    DO i=1,Nx
     !check for phi-crossing y-axis
     IF( PSI(i,1) * PSI(i+1,1) .lt. 0.0d0 )THEN
       !from (x,PSI(x)) point below and above PSI=0 axis, construct a line
       y1=PSI(i,1)
       x1=(i-1)*dx
       y2=PSI(i+1,1)
       x2=i*dx
       exit
     ENDIF
    END DO
    !interface defined as point where line crosses PSI=0
    x_cross=x1 - ( ( x2-x1 ) / ( y2-y1) ) * y1
    write(12,*)  (t_loop-1)*dt,x_cross

    !!!!!!!!!!!!!!!!!!! I/O  (print out PSI)!!!!!!!!!!!!!!!!!!!!!!!!
    if(mod(t_loop,file_skip)+1 == 1)then
       call chari(t_loop,cn,il)
       !open file to store phase and concentration fields
        open(1,file='out_'//cn(1:il),status='unknown')
        call print_2Dfields( PSI, C, Nx, Ny, 1)
        close(1)

       !open file to store centreline solute concentration 
        open(22,file='centre_conc_'//cn(1:il),status='unknown')
        call print_1Dfield( C(:,1), Nx, 22 )
        close(22)
    end if

 end do !time loop
 
  !close file that stores interface position
  close(12)
 end subroutine calculate 

  !-------------------
  !interpolated phase-dependent portion of diffusivity
  function Q(x)
  real*8 ::Q
  real*8 ::x
   Q=(1.0d0-x) / (1.0d0+partition-(1.0d0-partition )*x)
  end function Q

  !-------------------------------------------------------------------
  subroutine set_aniso_parameters
  !unsing derivatives in the x (DERX) and y (DERY) directions, compute anisotropy (see text)
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
  !finite differenced laplacian (see text)
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
 !enforce periodic poinday conditions
  subroutine PERIODIC(A)
  real*8, dimension(0:,0:) ::A
    A(Nx+1,:)    = A(1,:)
    A(0,:)      = A(Nx,:)
    A(:,Ny+1)    = A(:,1)
    A(:,0)      = A(:,Ny)
  end subroutine periodic

 !-----------------------
 !enforce zero flux bounary conditions
  subroutine ZERO_FLUX(A)
  real*8, dimension(0:,0:) ::A
    A(Nx+1,:)    = A(Nx-1,:)
    A(0,:)      = A(2,:)
    A(:,Ny+1)    = A(:,Ny-1)
    A(:,0)      = A(:,2)
  end subroutine ZERO_flux
  !!!!!!!!!!!!!!!!!!!!!!!




END MODULE SOLVER
