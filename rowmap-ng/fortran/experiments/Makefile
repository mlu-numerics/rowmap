# Builds the experiment drivers.
#   driver      : original ../rowmap.F (unchanged)          -> sft 2, 3 and autonomous runs
#   driver_sft  : patched copy rowmap_sft.F (see PATCH(sft)) -> sft 4, 5
FC      = gfortran
FFLAGS  = -O2 -g -fno-fast-math
LEGACY  = -std=legacy -w
LIBS    = -lopenblas
B       = build

all: $(B)/driver $(B)/driver_sft $(B)/sweep

$(B)/sweep: sweep.f90 $(B)/rowmap.o | $(B)
	$(FC) $(FFLAGS) -J$(B) -o $@ sweep.f90 $(B)/rowmap.o $(LIBS)

$(B):
	mkdir -p $(B)

$(B)/problems.o: problems.f90 | $(B)
	$(FC) $(FFLAGS) -J$(B) -c $< -o $@

$(B)/driver.o: driver.f90 $(B)/problems.o
	$(FC) $(FFLAGS) -J$(B) -c $< -o $@

$(B)/rowmap.o: ../rowmap.F | $(B)
	$(FC) $(FFLAGS) $(LEGACY) -c $< -o $@

$(B)/rowmap_sft.o: rowmap_sft.F | $(B)
	$(FC) $(FFLAGS) $(LEGACY) -c $< -o $@

$(B)/driver: $(B)/driver.o $(B)/problems.o $(B)/rowmap.o
	$(FC) -o $@ $^ $(LIBS)

$(B)/driver_sft: $(B)/driver.o $(B)/problems.o $(B)/rowmap_sft.o
	$(FC) -o $@ $^ $(LIBS)

clean:
	rm -rf $(B)

.PHONY: all clean
