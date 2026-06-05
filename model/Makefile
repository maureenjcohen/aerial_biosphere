FC      = gfortran
FFLAGS  = -O2 -Wall -Wextra -fcheck=all -fbacktrace
FFLAGS_FAST = -O3 -march=native

EXE     = bio_1d
SRCS    = bio_params.F90 bio_atmosphere.F90 bio_model.F90 bio_driver.F90
OBJS    = $(SRCS:.F90=.o)
MODS    = bio_params.mod bio_atmosphere.mod bio_model.mod

.PHONY: all clean fast

all: $(EXE)

fast: FFLAGS = $(FFLAGS_FAST)
fast: $(EXE)

$(EXE): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^

bio_params.o: bio_params.F90
	$(FC) $(FFLAGS) -c $<

bio_atmosphere.o: bio_atmosphere.F90 bio_params.o
	$(FC) $(FFLAGS) -c $<

bio_model.o: bio_model.F90 bio_params.o bio_atmosphere.o
	$(FC) $(FFLAGS) -c $<

bio_driver.o: bio_driver.F90 bio_params.o bio_atmosphere.o bio_model.o
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f $(OBJS) $(MODS) $(EXE)
