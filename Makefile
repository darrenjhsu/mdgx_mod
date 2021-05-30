
include ../config.h

MDGXWRAP_CFILES = Wrappers.c

MDGXWRAP_OBJS = Wrappers.o

yes:  install

no:
	@echo "Skipping installation of mdgx"

cuda: $(BINDIR)/mdgx.cuda$(SFX)

noparallel:
	@echo "Skipping installation of mdgx"

install: $(BINDIR)/mdgx$(SFX) $(MDGXWRAP_OBJS) $(LIBDIR)/libmdgx$(SHARED_SUFFIX)

parallel: $(BINDIR)/mdgx.MPI$(SFX)

uninstall:
	/bin/rm -f $(BINDIR)/mdgx$(SFX) $(BINDIR)/mdgx.MPI$(SFX) $(BINDIR)/mdgx.cuda$(SFX)
	/bin/rm -f $(LIBDIR)/libmdgx.* $(LIBDIR)/libmdgxMPI.*
	/bin/rm -f $(INCDIR)/mdgx.h

MDGX_CFILES = \
	    CompFrc.c \
	    BSpline.c \
	    Grid.c \
	    Random.c \
	    mdgxVector.c \
	    Matrix.c \
	    ChargeMap.c \
	    pmeRecip.c \
	    mleRecip.c \
	    pmeDirect.c \
	    CrdManip.c \
	    CellManip.c \
	    Topology.c \
	    Trajectory.c \
	    SpecialMath.c \
	    Nonbonded.c \
	    Bonded.c \
	    Parse.c \
	    Command.c \
	    Constraints.c \
	    Thermostats.c \
	    Barostats.c \
	    Integrator.c \
	    Timings.c \
	    Manual.c \
	    VirtualSites.c \
	    Buckingham.c \
	    ThermoDyn.c \
	    ChargeFit.c \
	    ParamOut.c \
	    ParamRead.c \
	    ParamFit.c \
	    IPolQ.c \
            ConfigSamp.c \
            SinglePointEval.c \
	    LoopBuilder.c \
	    Peptide.c \
	    Gpu.c \
	    Restraints.c \
	    Debug.c \
	    BroadcastCommand.c \
	    MPIMap.c \
	    MPITypeCast.c \
	    ptrajmask.c \
	    AmberNetcdf.c \
	    mdgx.c

MDGX_OBJS = \
	    CompFrc.o \
	    BSpline.o \
	    Grid.o \
	    Random.o \
	    mdgxVector.o \
	    Matrix.o \
	    ChargeMap.o \
	    pmeRecip.o \
	    mleRecip.o \
	    pmeDirect.o \
	    CrdManip.o \
	    CellManip.o \
	    Topology.o \
	    Trajectory.o \
	    SpecialMath.o \
	    Nonbonded.o \
	    Bonded.o \
	    Parse.o \
	    Command.o \
	    Constraints.o \
	    Thermostats.o \
	    Barostats.o \
	    Integrator.o \
	    Timings.o \
	    Manual.o \
	    VirtualSites.o \
	    Buckingham.o \
	    ThermoDyn.o \
	    ChargeFit.o \
	    ParamOut.o \
	    ParamRead.o \
	    ParamFit.o \
	    IPolQ.o \
	    ConfigSamp.o \
            SinglePointEval.o \
	    LoopBuilder.o \
	    Peptide.o \
	    Gpu.o \
	    Restraints.o \
	    Debug.o \
	    BroadcastCommand.o \
	    MPIMap.o \
	    MPITypeCast.o \
	    ptrajmask.o \
	    AmberNetcdf.o \
	    mdgx.o

MDGX_HEADERS = \
	       CompFrc.h \
	       BSpline.h \
	       Grid.h \
	       Random.h \
	       mdgxVector.h \
	       Matrix.h \
	       ChargeMap.h \
	       pmeRecip.h \
	       mleRecip.h \
	       pmeDirect.h \
	       CrdManip.h \
	       CellManip.h \
	       Topology.h \
	       Trajectory.h \
	       SpecialMath.h \
	       Nonbonded.h \
	       Bonded.h \
	       Parse.h \
	       Command.h \
	       Constraints.h \
	       Thermostats.h \
	       Barostats.h \
	       Integrator.h \
	       Timings.h \
	       Manual.h \
	       VirtualSites.h \
	       Buckingham.h \
	       ThermoDyn.h \
	       ChargeFit.h \
	       ParamOut.h \
	       ParamRead.h \
	       ParamFit.h \
	       IPolQ.h \
               ConfigSamp.h \
               SinglePointEval.h \
               Peptide.h \
	       Gpu.h \
	       Restraints.h \
	       Debug.h \
	       BroadcastCommand.h \
	       MPIMap.h \
	       MPITypeCast.h \
	       ptrajmask.h \
	       AmberNetcdf.h 

MDGX_CUDA_FILES = \
	          ArraySimulator.cu

MDGX_CUDA_OBJS = \
	         ArraySimulator.o \

MDGX_CUDA_HEADERS = \
	            ArraySimulator.h

$(BINDIR)/mdgx$(SFX) : $(FFTW3) $(MDGX_OBJS)
	@echo "[MDGX]  CC  $@"
	$(VB)$(CC) $(COPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) $(LDFLAGS) \
	-o $@ $(MDGX_OBJS) -L$(LIBDIR) $(FLIBS_FFTW3) $(NETCDFLIB) $(LM)

$(BINDIR)/mdgx.cuda$(SFX) : $(FFTW3) $(MDGX_CUDA_OBJS) $(MDGX_OBJS)
	@echo "[MDGX]  CC  $@"
	$(VB)$(CC) $(COPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) $(LDFLAGS) \
	-o $@ $(MDGX_CUDA_OBJS) $(MDGX_OBJS) $(MDGX_CU_LIBS) -L$(LIBDIR) $(FLIBS_FFTW3) $(NETCDFLIB) $(LM)

$(BINDIR)/mdgx.MPI$(SFX) : $(MDGX_OBJS)
	@echo "[MDGX]  CC  $@"
	$(VB)$(CC) $(COPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) $(LDFLAGS) \
	-o $@ $(MDGX_OBJS) -L$(LIBDIR) $(FLIBSF) $(FLIBS_FFTW3) $(NETCDFLIB) $(LM)

$(LIBDIR)/libmdgx$(SHARED_SUFFIX): $(MDGX_OBJS) $(MDGXWRAP_OBJS) $(BINDIR)/ucpp
	@echo "[MDGX]  CC  $@"
	$(VB)$(CC) $(MAKE_SHARED) -o $@ $(CFLAGS) $(COPTFLAGS) $(MDGX_OBJS) $(MDGXWRAP_OBJS) \
	-L$(LIBDIR) $(FLIBSF) $(FLIBS_FFTW3) $(NETCDFLIB) $(LM)
	$(VB)cp -p mdgxapi.h $(INCDIR)/mdgx.h

# mdgxapi.h is a pruned version of the mdgxhcat.sh output;
# for problem cases try using the mdgxhcat.sh output via this target.
mdgx.h:
	./mdgxhcat.sh $(CPP)
	mv mdgx.h $(INCDIR)
	#cp -p mdgxapi.h $(INCDIR)/mdgx.h

$(FFTW3):
	cd ../fftw-3.3 && $(MAKE) && $(MAKE) -j 1 install

$(BINDIR)/ucpp:
	cd ../ucpp-1.3 && $(MAKE) install

.c.o:
	@echo "[MDGX]  CC  $<"
	$(VB)$(CC) -c $(COPTFLAGS) $(MDGX_CFLAGS) $(CFLAGS) $(AMBERCFLAGS) $(MDGX_CU_DEFINES) $(MDGX_CU_INCLUDES) -I$(INCDIR) $(NETCDFINC) -o $@ $<

# This is puerile, but simply enumerate the CUDA units here.
ArraySimulator.o : ArraySimulator.cu
	@echo "[MDGX]  NVCC  $<"
	$(VB)$(NVCC) $(MDGX_CU_DEFINES) $(MDGX_CU_INCLUDES) -c -I$(INCDIR) $(NETCDFINC) -o $@ $<

# Dependencies
ArraySimulator.cu : KernelMacros.h kDynamics.h

# This is puerile, but simply enumerate the CUDA units here.
ArraySimulator.o : ArraySimulator.cu
	@echo "[MDGX]  NVCC  $<"
	$(VB)$(NVCC) -DCUDA -c -I$(INCDIR) $(NETCDFINC) -o $@ $<

clean:
	/bin/rm -f $(MDGX_OBJS) $(MDGX_CUDA_OBJS) $(MDGXWRAP_OBJS) 
