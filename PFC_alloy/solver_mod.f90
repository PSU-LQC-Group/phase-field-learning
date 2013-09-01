!!
!!Numerical integration of simple alloy PFC model from text by Provatas and Elder 
! Free energy in Section 9.2. and Eqs. 9.50 for the dynamcis (without noise!)
! using central difference (spherical laplacian) formula for all derivatives
!
MODULE SOLVER
USE VARIABLES
USE UTIL
IMPLICIT NONE

CONTAINS

 subroutine calculate
 
 !set initial conditions:
 !!NOTE: F-->density & G-->concentration
  
 if(restart==0)then
 !!!!!!!!!!case 1: set initial condition 

  !!initialize density 
  call srand(idum)
  mid_sys1=Nx/4
  mid_sys2=3*Nx/4
  qo=0.8660d0
  do i=1,Nx
   do j=1,Ny
    d1=sqrt((i*dx-mid_sys1*dx)**2+(j*dx-mid_sys1*dx)**2)-Rad*dx
	d2=sqrt((i*dx-mid_sys2*dx)**2+(j*dx-mid_sys2*dx)**2)-Rad*dx
	 !if(i.gt.mid_sys1-Rad.and.i.lt.mid_sys1+Rad)then !case 1: two long strips of solid
	 if(d1.le.0.0d0) then  !case 2: two circular seeds
      ! small sinusoidal seed at left
      F(i,j) = 0.5*( cos(2*qo*(j-1)*dx/sqrt(3.0))/2.0 &
      &      - cos(qo*(i-1)*dx)*cos(qo*(j-1)*dx/sqrt(3.0d0)) )+rno
  	 !else if(i.gt.mid_sys2-Rad.and.i.lt.mid_sys2+Rad)then !case 1
	 else if(d2.le.0.0d0)then !case 2:
      ! small sinusoidal seed at right, misoreinted by theta
      xr=  (i-1)*dx*cos(theta)+(j-1)*dx*sin(theta)
      yr= -(i-1)*dx*sin(theta)+(j-1)*dx*cos(theta)
      F(i,j) = 0.5*( cos(2*qo*yr/sqrt(3.0))/2.0 &
      &      - cos(qo*xr)*cos(qo*yr/sqrt(3.0d0)) )+rno
     else
      !  everywhere else
      F(i,j)=0.001*gasdev2()+rno
    endif
   end do  
  end do  

  !!initialize concentration
  do i=1,Nx
   do j=1,Ny
    d1=sqrt((i*dx-mid_sys1*dx)**2+(j*dx-mid_sys1*dx)**2)-Rad*dx
	d2=sqrt((i*dx-mid_sys2*dx)**2+(j*dx-mid_sys2*dx)**2)-Rad*dx
    !if(i.gt.mid_sys1-Rad.and.i.lt.mid_sys1+Rad)then !case 1: two long strips of solid
    if(d1.le.0.0d0) then  !case 2: two circular seeds 
	   G(i,j)=0.001*gasdev2()+cs
    !else if(i.gt.mid_sys2-Rad.and.i.lt.mid_sys2+Rad)then !case 1
	else if(d2.le.0.0d0)then !case 2:
      G(i,j)=0.001*gasdev2()+cs
    else
      G(i,j)=0.001*gasdev2()+co
    endif
   end do
  end do

 else if(restart==1)then
 !!!!!!!case 2: srtart from previous file --which you arrange to call "restart_file"
  open(1,file='restart_file',status='unknown')
  do i=1,Nx
   do j=1,Ny
    read(1,12) F(i,j),G(i,j)
	12 FORMAT(E23.16,1x,E23.16)
   end do
  end do 
  close(1)
 
 endif

  !!initialize periodi buffer arrays (see text)
  call PERIODIC(F)
  call PERIODIC(G)

  !print intial ic into file
  if(restart==0)
  open(1,file='out_0',status='unknown')
    call print_2Dfields(F, G, Nx, Ny, 1)
  close(1)
  endif
  
  !initialize Bol be in the "firt" jump time sequence 
  B0l_jump=1

  !start keeping record of temperature changes (i.e. changes in B0l)
  !open(16,file='B0l_record',status='unknown')
  print*, 'DB0l for next ',B0l_change_time,' steps=',B0l-Box

 !!!START TIME STEPPING
 do t_loop=1, tmax

    !!!!!!!!!!!!!!! adjust temperature !!!!!!!!!!!!!!!!!!!!!!!!
    if(mod(t_loop,B0l_change_time+1)+1 == 1)then
     B0l_jump=B0l_jump+1

     if(B0l_jump.le.B0l_stages)then
      B0l=B0l-DB0l
	  print*, 'At time step=',t_loop
      print*, 'DB0l for next ',B0l_change_time,' steps=',B0l-Box
     endif

    endif


    !!!!!!!!!!!!!!! density equation update !!!!!!!!!!!!!!!!!!!!

    !1a. periodic BC for F & G and grad2 arrays for density equation
    call PERIODIC(F)
    call PERIODIC(G)
    call NABLA2(F,grad2_n)
    call NABLA2(F*G,grad2_nc)


    !1b. periodic BC for gra2_n & grad2_nc and  grad4 arrays for density equation
    call PERIODIC(grad2_n)
    call PERIODIC(grad2_nc)
    call NABLA2(grad2_n,grad4_n)
    call NABLA2(grad2_nc,grad4_nc)
	
    !1c. calculate the solute dependent coefficient Bl
    do i=1,Nx
     do j=1,Ny
      Bl(i,j)=B0l+B2l*G(i,j)**2
     end do 
    end do 

    !1d. update right hand side of density chemical potential
    do i=1, Nx
     do j=1,Ny
      RHS(i,j)= Bl(i,j)*F(i,j)+Box*(2*grad2_n(i,j)+grad4_n(i,j)) & 
  &            +2*eta*Box*( G(i,j)*(grad2_n(i,j)+grad4_n(i,j)) &
  &            +grad2_nc(i,j)+grad4_nc(i,j) ) &
  &            - tt*F(i,j)**2 + vv*F(i,j)**3
     end do 
    end do

    !1e. periodic BC for grad2 of RHS and grad2 of RHS
    call PERIODIC(RHS)
    call NABLA2(RHS,grad2_mn)

    !1f. update rho-equation
    do i=1, Nx
     do j=1,Ny
       F(i,j) = F(i,j) + dt*Mn*grad2_mn(i,j)
     end do
    end do

    !!!!!!!!!!!!!!!!! concentration eq update !!!!!!!!!!!!!!!!!!!!!!!

    !2a. re-set BC for F (G still periodic) & grad2 arrays for conc equation
    call PERIODIC(F)
	!G still remains periodic from above update
    call NABLA2(F,grad2_n)
    call NABLA2(G,grad2_c)

    !2b. periodic BC for gra2_n & grad4 arrays for concentration equation
    call PERIODIC(grad2_n)
    call NABLA2(grad2_n,grad4_n)

    !2c. calculate right hand side of density eqn
    do i=1, Nx
     do j=1,Ny
      RHS(i,j) = 2*Box*eta*F(i,j)*(grad2_n(i,j)+grad4_n(i,j)) &
   &           + (ww+B2l*F(i,j)**2)*G(i,j) &
   &           + uu*G(i,j)**3 - RKK*grad2_c(i,j)
     end do
    end do 

    !2d. periodic BC for RHS and  grad2 of RHS
    call PERIODIC(RHS)
    call NABLA2(RHS,grad2_mc)

    !2e. integrate ofconcentration (\delta N) equation
    do i=1,Nx
     do j=1,Ny
       G(i,j) = G(i,j) + dt*Mn*grad2_mc(i,j)
     end do
    end do


    !!!!!!!!!!!!!!!!!!! I/O !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if(mod(t_loop,file_skip)+1 == 1)then
       call chari(t_loop+restart_time,cn,il)
       open(1,file='out_'//cn(1:il),status='unknown')
         call print_2Dfields( F, G, Nx, Ny, 1)
       close(1)
    end if


 end do !time loop

 !close(16)
 end subroutine calculate 


  !!!!!!!!!!!!!!!!!!!!!!!!
  subroutine NABLA2(A,GR2)
  real*8, dimension(0:,0:) ::A,GR2
   GR2=0.0d0
   do i=1,Nx
    do j=1,Ny
       GR2(i,j) = ( (A(i+1,j) + A(i-1,j) + A(i,j+1) +A(i,j-1))/2.0  &
              & +   (A(i+1,j+1)+A(i-1,j+1)+A(i+1,j-1) +A(i-1,j-1))/4.0 &
              & -   3*A(i,j) )* t4   
 
    end do
   end do
  end subroutine NABLA2

  subroutine PERIODIC(A)
  real*8, dimension(0:,0:) ::A
    A(Nx+1,:)    = A(1,:)
    A(0,:)      = A(Nx,:)
    A(:,Ny+1)    = A(:,1)
    A(:,0)      = A(:,Ny)
  end subroutine periodic
  !!!!!!!!!!!!!!!!!!!!!!!

END MODULE SOLVER
