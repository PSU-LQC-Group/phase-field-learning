PROGRAM MANAGER
USE VARIABLES
USE SOLVER
IMPLICIT NONE

!initialize variables and read in model/material parameters
call read_globals

!run time evolution of Model A
call calculate

END PROGRAM MANAGER
