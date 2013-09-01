PROGRAM MANAGER
USE VARIABLES
USE SOLVER
IMPLICIT NONE

!intialize variables and read in model/material parameters
call read_globals

!run time evolution of model C for alloys (dilute alloy model in text)
call calculate

END PROGRAM MANAGER
