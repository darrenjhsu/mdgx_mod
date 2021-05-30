project(mdgx C)

#Sources common to mdgx and libmdgx
set(MDGX_SOURCES 	    
	CompFrc.c 
	BSpline.c 
	Grid.c 
	Random.c 
	mdgxVector.c 
	Matrix.c 
	ChargeMap.c 
	pmeRecip.c 
	mleRecip.c 
	pmeDirect.c 
	CrdManip.c 
	CellManip.c 
	Topology.c 
	Trajectory.c 
	SpecialMath.c 
	Nonbonded.c 
	Bonded.c 
	Parse.c 
	Command.c 
	Constraints.c 
	Thermostats.c 
	Barostats.c 
	Integrator.c 
	Timings.c 
	Manual.c 
	VirtualSites.c 
	Buckingham.c 
	ThermoDyn.c 
	ChargeFit.c 
	ParamOut.c 
	ParamRead.c 
	ParamFit.c 
	IPolQ.c 
	ConfigSamp.c 
	SinglePointEval.c 
	LoopBuilder.c 
	Peptide.c 
	Gpu.c 
	Restraints.c 
	Debug.c 
	BroadcastCommand.c 
	MPIMap.c 
	MPITypeCast.c 
	ptrajmask.c 
	AmberNetcdf.c 
	mdgx.c)


set(MDGX_CUDA_SOURCES ArraySimulator.cu)
	  
# Sources included only in libmdgx    
set(MDGX_WRAPPER_SOURCES Wrappers.c)

#------------------------------------------------------------------------------------------

message("Let's compile some mdgx things! zzzzzzzzzzzzzz")
add_compile_options(${OPT_CFLAGS})
include_directories($<TARGET_PROPERTY:fftw,INTERFACE_INCLUDE_DIRECTORIES> $<TARGET_PROPERTY:netcdf,INTERFACE_INCLUDE_DIRECTORIES>)

#common objects in the executable and library
add_library(mdgx_common_obj OBJECT ${MDGX_SOURCES})
make_pic_if_needed(mdgx_common_obj)

#mdgx executable
add_executable(mdgx $<TARGET_OBJECTS:mdgx_common_obj>)
target_link_libraries(mdgx fftw netcdf netlib)
install(TARGETS mdgx DESTINATION ${BINDIR}) 

#libmdgx
add_library(libmdgx $<TARGET_OBJECTS:mdgx_common_obj> ${MDGX_WRAPPER_SOURCES})
target_link_libraries(libmdgx netcdf fftw netlib)
remove_prefix(libmdgx)
install_libraries(libmdgx)

	
#install pruned version of generated header
install(FILES mdgxapi.h DESTINATION ${INCDIR} RENAME mdgx.h) 

# --------------------------------------------------------------------
# Parallel versions
# --------------------------------------------------------------------
if(MPI)
	make_mpi_version(mdgx_common_obj mdgx_common_obj_mpi LANGUAGES C)
	 
	make_mpi_version(mdgx mdgx.MPI LANGUAGES C SWAP_SOURCES $<TARGET_OBJECTS:mdgx_common_obj> TO $<TARGET_OBJECTS:mdgx_common_obj_mpi> INSTALL)
	make_mpi_version(libmdgx libmdgx_mpi LANGUAGES C SWAP_SOURCES $<TARGET_OBJECTS:mdgx_common_obj> TO $<TARGET_OBJECTS:mdgx_common_obj_mpi> INSTALL)
	
	if(TARGET_WINDOWS)
		# MPI MDGX uses gethostname to identify machines
		targets_link_libraries(libmdgx_mpi mdgx.MPI LIBRARIES ws2_32)
	endif()
endif()

if(OPENMP)
	make_openmp_version(mdgx_common_obj mdgx_common_obj_omp LANGUAGES C)
	 
	make_openmp_version(mdgx mdgx.OMP LANGUAGES C SWAP_SOURCES $<TARGET_OBJECTS:mdgx_common_obj> TO $<TARGET_OBJECTS:mdgx_common_obj_omp> INSTALL)
	make_openmp_version(libmdgx libmdgx_omp LANGUAGES C SWAP_SOURCES $<TARGET_OBJECTS:mdgx_common_obj> TO $<TARGET_OBJECTS:mdgx_common_obj_omp> INSTALL)
endif()

if (CUDA AND EXISTS "${CUDA_nvidia-ml_LIBRARY}") # make sure we have the extra nvidia-ml library we need
	#set(CMAKE_C_FLAGS -DCUDA)
	message("Let's compile some mdgx.cuda things! xxxxxxxxxxxxxx")
	# no need to make a common object library since there is no shared library version for cuda
	cuda_add_executable(mdgx.cuda ${MDGX_SOURCES} ${MDGX_CUDA_SOURCES} OPTIONS -DCUDA)
	target_compile_options(mdgx.cuda PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-DCUDA> $<$<COMPILE_LANGUAGE:C>:-DCUDA>)
	
	#cuda_add_executable(mdgx.cuda ${MDGX_CUDA_SOURCES} OPTIONS -DCUDA)

	target_link_libraries(mdgx.cuda curand cudadevrt nvidia-ml fftw netlib amber_common netcdff netcdf)

	install(TARGETS mdgx.cuda DESTINATION ${BINDIR})

endif()

#------------------------------------------------------------------------------------------
