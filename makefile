LIBS="-lblas -llapack"
demo: driver.F rowmap.F 
	gfortran -o $@ $^ $(LIBS)
