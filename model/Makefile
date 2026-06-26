FC      = gfortran
FFLAGS  = -O2 -Wall -Wextra -fcheck=all -fbacktrace
FFLAGS_FAST = -O3 -march=native

EXE     = bio_1d
SRCS    = bio_params.F90 bio_atmosphere.F90 bio_model.F90 bio_driver.F90
OBJS    = $(SRCS:.F90=.o)
MODS    = bio_params.mod bio_atmosphere.mod bio_model.mod

# Venus cloud biosphere model: environment module + biology module + driver
EXE_CLOUD   = bio_cloud
CLOUD_SRCS  = bio_cloud_model.F90 bio_venus.F90 bio_venus_driver.F90
CLOUD_OBJS  = $(CLOUD_SRCS:.F90=.o)
CLOUD_MODS  = bio_cloud_model.mod bio_venus.mod

.PHONY: all clean fast cloud

all: $(EXE) $(EXE_CLOUD)

fast: FFLAGS = $(FFLAGS_FAST)
fast: $(EXE) $(EXE_CLOUD)

cloud: $(EXE_CLOUD)

$(EXE): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^

$(EXE_CLOUD): $(CLOUD_OBJS)
	$(FC) $(FFLAGS) -o $@ $^

bio_cloud_model.o: bio_cloud_model.F90
	$(FC) $(FFLAGS) -c $<

bio_venus.o: bio_venus.F90 bio_cloud_model.o
	$(FC) $(FFLAGS) -c $<

bio_venus_driver.o: bio_venus_driver.F90 bio_cloud_model.o bio_venus.o
	$(FC) $(FFLAGS) -c $<

bio_params.o: bio_params.F90
	$(FC) $(FFLAGS) -c $<

bio_atmosphere.o: bio_atmosphere.F90 bio_params.o
	$(FC) $(FFLAGS) -c $<

bio_model.o: bio_model.F90 bio_params.o bio_atmosphere.o
	$(FC) $(FFLAGS) -c $<

bio_driver.o: bio_driver.F90 bio_params.o bio_atmosphere.o bio_model.o
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f $(OBJS) $(MODS) $(EXE) $(EXE_CLOUD) $(CLOUD_OBJS) $(CLOUD_MODS)
