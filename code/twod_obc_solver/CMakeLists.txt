cmake_minimum_required(VERSION 3.18.4)

################################################################################################################################
##DEVELOPER SECTION#############################################################################################################
##in this section software developers need to make changes when debugging#######################################################
################################################################################################################################
# please name your project accordingly
set(MYPROJECTNAME "gragles_twod_obc_solver")
message([STATUS] "Projectname is ${MYPROJECTNAME}")
project(${MYPROJECTNAME} LANGUAGES C CXX)
set(CMAKE_BUILD_DIR "build")

set(MYPROJECTPATH "${CMAKE_CURRENT_LIST_DIR}/../..")
set(MYEIGENPATH "${MYPROJECTPATH}/code/thirdparty/mandatory/eigen/eigen")
set(MYFFTWPATH "${MYPROJECTPATH}/code/thirdparty/mandatory/fftw")
set(MYIMKLPATH "/opt/intel/oneapi/mkl/2024.0")  # aka MKLROOT, if installed in default locatio, but version may change
set(MYHDFPATH "${MYPROJECTPATH}/code/thirdparty/mandatory/hdf5/CMake-hdf5-1.14.2/build/HDF5-1.14.2-Linux/HDF_Group/HDF5/1.14.2")
set(MYVOROXXPATH "${MYPROJECTPATH}/code/thirdparty/mandatory/voroxx/voro")
set(MYJEMALLOCPATH "${MYPROJECTPATH}/code/thirdparty/optional/jemalloc/jemalloc")
set(MYTOOLSRCPATH "${MYPROJECTPATH}/code/twod_obc_solver/src/cxx")
# pull general information about external libraries and paths related to paraprobe
message([STATUS] "MYPROJECTPATH: ${MYPROJECTPATH}")
message([STATUS] "MYEIGENPATH: ${MYEIGENPATH}")
message([STATUS] "MYFFTWPATH: ${MYFFTWPATH}")
message([STATUS] "MYIMKLPATH: ${MYIMKLPATH}")
message([STATUS] "MYHDF5PATH: ${MYHDFPATH}")
message([STATUS] "MYVOROXXPATH: ${MYVOROXXPATH}")
message([STATUS] "MYJEMALLOCPATH: ${MYJEMALLOCPATH}")
message([STATUS] "MYTOOLPATH: ${MYTOOLPATH}")

# define compiler and Fourier transform library to use
set(EMPLOY_GNUCOMPILER OFF)  # use GNU C/C++ compiler for general purpose if ON, otherwise using Intel compiler (best for Intel CPU)
set(EMPLOY_FFTW_LIBRARY OFF)  # use FFTW library for DFT computations for general purpose architecture but that may be slower than using the MKL library

message([STATUS] "MYCCC_COMPILER: __${CMAKE_C_COMPILER}__")
message([STATUS] "MYCXX_COMPILER: __${CMAKE_CXX_COMPILER}__")

# which FFT library to use
if(EMPLOY_FFTW_LIBRARY)
	ADD_DEFINITIONS("-DUSE_FFTW")
	message([STATUS] "Using Fast Fourier Transform in the West (FFTW) library")
else()
	ADD_DEFINITIONS("-DUSE_MKL")
	message([STATUS] "Using Intel Math Kernel Library Fast Fourier transformations")
endif()
# DUSE_FFTW approximately an order of magnitude overhead compared to DUSE_MKL

#choose optimization level
##-O0 nothing, debug purposes, -O1 moderate optimization, -O2 -O3 for production level up to aggressive architecture specific non-portable optimization
if(EMPLOY_GNUCOMPILER)
	set(MYOPTLEVEL "-O2")
else()
	set(MYOPTLEVEL "-O2 -DMKL_ILP64 -m64")
endif()
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MYOPTLEVEL}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MYOPTLEVEL}")
message([STATUS] "We utilize optimization level ${MYOPTLEVEL}")


################################################################################################################################
##END OF INTERACTION FOR NON PRO USERS##########################################################################################
##here advanced users might want/need to make modifications if they use non default places for thirdparty libraries#############
################################################################################################################################
#HDF5 local installation for advanced I/O, collecting metadata and bundle analysis results together
include_directories("${MYHDFPATH}/include")
set(MYHDFLINKFLAGS "-L${MYHDFPATH}/lib/ ${MYHDFPATH}/lib/libhdf5_hl.a ${MYHDFPATH}/lib/libhdf5.a ${MYHDFPATH}/lib/libz.a -ldl")
if(EMPLOY_FFTW_LIBRARY)
	set(MYFFTW_FLOAT_LINKFLAGS "-L${MYFFTWPATH}/f32/fftw-3.3.10/lib/ ${MYFFTWPATH}/f32/fftw-3.3.10/lib/libfftw3f.a")
	set(MYFFTW_DOUBLE_LINKFLAGS "-L${MYFFTWPATH}/f64/fftw-3.3.10/lib/ ${MYFFTWPATH}/f64/fftw-3.3.10/lib/libfftw3.a")
else()
	# check https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-link-line-advisor.html for more details
	set(MYIMKL_LINKFLAGS "-L${MYIMKLPATH}/lib -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl")
endif()
set(MYJEMALLOC_LINKFLAGS "-L${MYJEMALLOCPATH}/local/lib/ ${MYJEMALLOCPATH}/local/lib/libjemalloc.a")
include_directories("${MYEIGENPATH}")
include_directories("${MYVOROXXPATH}/src")

# https://gitlab.kitware.com/cmake/cmake/-/issues/22831
#setting up compiler-specifics
if(EMPLOY_GNUCOMPILER)
	set(MYVERBOSE "-Wall") # all warnings, -fopt-info for optimization details only
	add_definitions("${MYOPTLEVEL}")
	add_definitions("${MYVERBOSE}")
	set(MYOMP "-fopenmp")
	add_definitions(${MYOMP})
	# add_definitions("-std=c++0x")
	add_definitions("-std=c++17")
else()
	set(MYVERBOSE "-Wall") # all warnings
	add_definitions("${MYOPTLEVEL}")
	add_definitions("${MYVERBOSE}")
	set(MYOMP "-fiopenmp")
	add_definitions(${MYOMP})
	# add_definitions("-std=c++0x")
	add_definitions("-std=c++17")
endif()

# list firstly the precompiled shared aka utils, secondly the tool-specific components, lastly the tool-specific main
add_executable(${MYPROJECTNAME}
	# ${MYTOOLSRCPATH}/allocator.cpp
	${MYTOOLSRCPATH}/../../../utils/src/cxx/PARAPROBE_HDF5Core.cpp
	${MYTOOLSRCPATH}/ggLS.cpp
	${MYTOOLSRCPATH}/mymath.cpp
	${MYTOOLSRCPATH}/SCORE_Random.cpp
	${MYTOOLSRCPATH}/grainhdl.cpp
	${MYTOOLSRCPATH}/box.cpp
	${MYTOOLSRCPATH}/marchingSquares.cpp
	${MYTOOLSRCPATH}/Settings.cpp
	${MYTOOLSRCPATH}/contourSector.cpp
	${MYTOOLSRCPATH}/grahamScan.cpp
	${MYTOOLSRCPATH}/junction.cpp
	${MYTOOLSRCPATH}/minimalisticBoundary.cpp
	${MYTOOLSRCPATH}/grainBoundary.cpp
	${MYTOOLSRCPATH}/misorientation.cpp
	${MYTOOLSRCPATH}/Quaternion.cpp
	${MYTOOLSRCPATH}/IterativeGrainScheduler.cpp
	${MYTOOLSRCPATH}/SquaresGrainScheduler.cpp

	${MYVOROXXPATH}/src/c_loops.cc
	${MYVOROXXPATH}/src/cell.cc
	# ${MYVOROXXPATH}/src/cmd_line.cc
	${MYVOROXXPATH}/src/common.cc
	${MYVOROXXPATH}/src/container.cc
	${MYVOROXXPATH}/src/container_prd.cc
	${MYVOROXXPATH}/src/pre_container.cc
	${MYVOROXXPATH}/src/unitcell.cc
	${MYVOROXXPATH}/src/v_base.cc
	# ${MYVOROXXPATH}/src/v_base_wl.cc
	${MYVOROXXPATH}/src/v_compute.cc
	# ${MYVOROXXPATH}/src/voro++.cc
	${MYVOROXXPATH}/src/wall.cc
)

# linking process, the target link libraries command is specific for each tool of the toolbox
if(EMPLOY_FFTW_LIBRARY)
	target_link_libraries(${MYPROJECTNAME} ${MYOMP} ${MYHDFLINKFLAGS} ${MYFFTW_FLOAT_LINKFLAGS} ${MYFFTW_DOUBLE_LINKFLAGS} -lnuma ${MYJEMALLOC_LINKFLAGS})
else()
	target_link_libraries(${MYPROJECTNAME} ${MYOMP} ${MYHDFLINKFLAGS} ${MYIMKL_LINKFLAGS} -lnuma ${MYJEMALLOC_LINKFLAGS})
endif()
