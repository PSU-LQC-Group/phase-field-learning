MODULE UTIL
USE VARIABLES
implicit none ;save

CONTAINS

!----------------------------------------------------------------------------
!print 2D array signal1 & signal2 of dimensions Nx by Ny into file numbered=file_no
subroutine print_2Dfields(signal1, signal2, Nx, Ny, file_no)
  real*8, dimension(0:,0:)   ::signal1
  real*8, dimension(0:,0:)   ::signal2
  integer		     ::file_no, Nx,Ny !Nx, Ny -number of gridpoints
  integer		     ::i,j
 
  ! print 2D field              !!!!!!!!!!!
  do i=1, Nx
   do j=1, Ny
     write(file_no,10) signal1(i,j), signal2(i,j)
     10 FORMAT(E23.16,1x,E23.16)
   end do
  end do
end subroutine print_2Dfields

!----------------------------------------------------------------------------
!1D analogue of previous routine 
subroutine print_1Dfield(signal1, Nx, file_no)
  real*8, dimension(0:)      ::signal1
  integer		     ::file_no, Nx !Nx -number of gridpoints
  integer		     ::i,j
 
  ! print 1D field1              !!!!!!!!!!!
  do i=1, Nx
     write(file_no,10) (i-1)*dx , signal1(i)
     10 FORMAT(E23.16,1x,E23.16)
  end do
end subroutine print_1Dfield

 !--------------------------------------------------
 !converts integer "i" into a character string, referenced as "cn(1:il)" in commands
 subroutine chari(i,ci,il)
    integer ::i,ii,i3,kk,k,j,j1,il
        real    ::ri
        character(len=9) ci
        character(len=10) str

        ii=i
4       if(ii.gt.999999999) then
                ri=ii
                ii=nint(ri/10)
                goto 4
        end if
        i3=ii
        str='0123456789'
        do 11 k=1,9
                j=10**k
                j1=10**(k-1)
                if((i3.ge.j1).and.(i3.lt.j)) il=k
11      continue
        do 22 k=il,1,-1
                kk=mod(ii,10)+1
                ci(k:k)=str(kk:kk)
                ii=ii/10
22      continue
        return
        end subroutine chari

!-------------------------------------------------
!generates uniform random numbers between [0,1]
FUNCTION ran(idum) 
IMPLICIT NONE 
INTEGER, PARAMETER :: K4B=selected_int_kind(9) 
INTEGER(K4B), INTENT(INOUT) :: idum 
REAL :: ran  
INTEGER(K4B), PARAMETER :: IA=16807,IM=2147483647,IQ=127773,IR=2836 
REAL, SAVE :: am
INTEGER(K4B), SAVE :: ix=-1,iy=-1,k 

if (idum <= 0 .or. iy < 0) then 
   am=nearest(1.0,-1.0)/IM 
   iy=ior(ieor(888889999,abs(idum)),1) 
   ix=ieor(777755555,abs(idum)) 
   idum=abs(idum)+1 
end if 
  ix=ieor(ix,ishft(ix,13)) 
  ix=ieor(ix,ishft(ix,-17)) 
  ix=ieor(ix,ishft(ix,5)) 
  k=iy/IQ 
  iy=IA*(iy-k*IQ)-IR*k 
  if (iy < 0) iy=iy+IM 
  ran=am*ior(iand(IM,ieor(ix,iy)),1) 
END FUNCTION ran

!-----------------------------------------
!generates gaussian distributed deviates with zero mean and unit STD, uses "ran"
FUNCTION gasdev(idum) 
INTEGER idum 
REAL gasdev  
INTEGER iset 
REAL fac,gset,rsq,v1,v2,ran1 
SAVE iset,gset 
DATA iset/0/ 

if (idum.lt.0) iset=0 
if (iset.eq.0) then 
1  v1=2.0*ran(idum)-1.0 
   v2=2.0*ran(idum)-1.0 
  rsq=v1**2+v2**2 
  if(rsq.ge.1..or.rsq.eq.0.)goto 1 
  fac=sqrt(-2.*log(rsq)/rsq) 
  gset=v1*fac 
  gasdev=v2*fac 
  iset=1 
else 
  gasdev=gset 
  iset=0 
endif 
return 
END FUNCTION gasdev

!-----------------------------------------
!gaussian deviates as above, may not work on all compilers
FUNCTION gasdev2() 
INTEGER idum 
REAL gasdev2
INTEGER iset 
REAL fac,gset,rsq,v1,v2,ran1 
SAVE iset,gset 
DATA iset/0/ 

if (idum.lt.0) iset=0 
if (iset.eq.0) then 
1  v1=2.0*rand()-1.0 
   v2=2.0*rand()-1.0 
  rsq=v1**2+v2**2 
  if(rsq.ge.1..or.rsq.eq.0.)goto 1 
  fac=sqrt(-2.*log(rsq)/rsq) 
  gset=v1*fac 
  gasdev2=v2*fac 
  iset=1 
else 
  gasdev2=gset 
  iset=0 
endif 
return 
END FUNCTION gasdev2

!--------------------------------------------------------------------
!generates random seed
  FUNCTION random_seed(NN, psi_01, init_seeds, seed_size, idum)
    INTEGER                   :: NN, asizex
    real*8, dimension(NN)      ::random_seed    
    INTEGER                   :: i, j, k, i_size, init_seeds, seed_size
    INTEGER                   :: i_seed, i_position, idum
    REAL*8                    :: psi_01
   asizex = NN
   random_seed = psi_01
   do k=1, init_seeds
     !random size generation for k-th seed
     i_size = 1+int(seed_size*ran(idum))
       if(i_size<1) i_size=1 
     ! print*, 'Seed, seed size ixj:', k, i_size, j_size
     !random position generation for k-th seed
      i_seed = 1+int(asizex * ran(idum))
       if(i_seed<1) i_seed=1 
     ! print*, 'Seed, seed position ixj:', k, i_seed, j_seed
     !boundary conditions 
     do i=1,i_size
        if(i_seed+i>asizex)then 
          i_position = ((i_seed+i)-asizex)
        else 
          i_position = i_seed+i
        end if
        ! print*, 'i_pos :', i_position
        random_seed(i_position) = 0.1*gasdev(idum)
     end do
    end do

    END FUNCTION random_seed

 !----------------------------------------------------------------------


End module UTIL
