F90 = /sw/bin/g95
F77 = /usr/bin/f77
FFLAGS =
########################################################################
MOD_OBJ =\
          variables_mod.o util_mod.o solver_mod.o

F90_OBJS = manager.o
##################################################################### Libraries

ALL_OBJS = $(MOD_OBJ) $(F90_OBJS) 

NDS_3D: $(ALL_OBJS) 
	$(F90) $(FFLAGS) $(ALL_OBJS)   -o modelC_pure

clean:  
	rm *.mod *.o *~ signal_* *.png *.csv modelC_pure

clean2:
	rm out_*
##################################################################### Rules
variables_mod.o        : variables_mod.f90
		         $(F90) $(FFLAGS)  -c variables_mod.f90
util_mod.o             :variables_mod.o util_mod.f90
		         $(F90) $(FFLAGS)  -c util_mod.f90
solver_mod.o           :variables_mod.o util_mod.o solver_mod.f90
		         $(F90) $(FFLAGS)  -c solver_mod.f90
manager.o 	       : manager.f90 variables_mod.o util_mod.o solver_mod.o
		         $(F90) $(FFLAGS)  -c manager.f90

######################################################################cleanup 




