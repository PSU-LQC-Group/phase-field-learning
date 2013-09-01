This codes simulate the simple binary Alloy Phase Field Crystal Model described in the text in section 9.2. It is written in C and uses fftw to run in parallel on multiple processors. It handles the situation eta=0.


_______________________INFORMATION--------------------------
Binary Alloy Phase Field Crystal Model Readme File
The input file alloy.in contains the following parameters (whose variable names in the PFC code are in the right hand column of the file):
densf	- name of files containing density data
concf	- name of files containing concentration data
Nx - number of grid points in the x-direction
Nx - number of grid points in the y-direction
dx - grid spacing in the x-direction
dt - size of time step
BL - difference BL-BX gives non-dimensional temperature
BL2 - 
BX - difference BL-BX gives non-dimensional temperature
RK - prefactor of gradient^2 concentration term in chemical potential
w - prefactor of quadratic concentration term in chemical potential
t - prefactor of cubic density term in chemical potential
u - prefactor of quartic concentration term in chemical potential
v - prefactor of quartic density term in chemical potential
rno - average mean density
co - average concentration in solid
cl - average concentration in liquid
nnend - last time step to be computed
nout - output file when time steps are an integer multiple of nout
nstart - first time step to be computed/read in
ntype - determines the initial condition to be used
qo - magnitude of reciprocal lattice vector
dLx - size of seed in x-direction
dLy - size of seed in y-direction
noise - amplitude of noise in concentration field


