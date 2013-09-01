PROGRAM LOADER
USE VARIABLES
USE SOLVER
IMPLICIT NONE

!initilialize and read in model/material parameters
call read_globals

!run time evolutio of alloy PFC model
call calculate

END PROGRAM LOADER
