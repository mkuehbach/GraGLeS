cmake_minimum_required(VERSION 3.18.4)

################################################################################################################################
##DEVELOPER SECTION#############################################################################################################
##in this section software developers need to make changes when debugging#######################################################
################################################################################################################################
#please name your project accordingly
set(MYPROJECTNAME "gragles_structure_generator")
message([STATUS] "Projectname is ${MYPROJECTNAME}")
project(${MYPROJECTNAME} LANGUAGES C CXX)
set(CMAKE_BUILD_DIR "build")


set(MYPROJECTPATH "${CMAKE_CURRENT_LIST_DIR}/../..")
set(MYEIGENPATH "${MYPROJECTPATH}/code/thirdparty/mandatory/eigen/eigen")
set(MYHDFPATH "${MYPROJECTPATH}/code/thirdparty/mandatory/hdf5/CMake-hdf5-1.14.2/build/HDF5-1.14.2-Linux/HDF_Group/HDF5/1.14.2")
set(MYVOROXXPATH "${MYPROJECTPATH}/code/thirdparty/mandatory/voroxx/voro")
set(MYTOOLSRCPATH "${MYPROJECTPATH}/code/structure_generator/src/cxx")

#pull general information about external libraries and paths related to paraprobe
message([STATUS] "MYPROJECTPATH: ${MYPROJECTPATH}")
message([STATUS] "MYEIGENPATH: ${MYEIGENPATH}")
message([STATUS] "MYHDF5PATH: ${MYHDFPATH}")
message([STATUS] "MYVOROXXPATH: ${MYVOROXXPATH}")
message([STATUS] "MYTOOLPATH: ${MYTOOLPATH}")
set(EMPLOY_GNUCOMPILER ON)
message([STATUS] "MYCCC_COMPILER: __${CMAKE_C_COMPILER}__")
message([STATUS] "MYCXX_COMPILER: __${CMAKE_CXX_COMPILER}__")


#choose optimization level
##-O0 nothing, debug purposes, -O1 moderate optimization, -O2 -O3 for production level up to aggressive architecture specific non-portable optimization
if(EMPLOY_GNUCOMPILER)
	set(MYOPTLEVEL "-O2")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MYOPTLEVEL}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MYOPTLEVEL}")
	message([STATUS] "We utilize optimization level ${MYOPTLEVEL}")
endif()

################################################################################################################################
##END OF INTERACTION FOR NON PRO USERS##########################################################################################
##here advanced users might want/need to make modifications if they use non default places for thirdparty libraries#############
################################################################################################################################
#HDF5 local installation for advanced I/O, collecting metadata and bundle analysis results together
include_directories("${MYHDFPATH}/include")
link_directories("${MYHDFPATH}/lib")
set(MYHDFLINKFLAGS "-L${MYHDFPATH}/lib/ ${MYHDFPATH}/lib/libhdf5_hl.a ${MYHDFPATH}/lib/libhdf5.a ${MYHDFPATH}/lib/libz.a -ldl")
include_directories("${MYEIGENPATH}")
include_directories("${MYVOROXXPATH}/src")

#setting up compiler-specifics
if(EMPLOY_GNUCOMPILER)
	# set(MYVERBOSE "-fopt-info") # optimization details only
	set(MYVERBOSE "-Wall") # all warnings
	add_definitions("${MYOPTLEVEL}")
	add_definitions("${MYVERBOSE}")
	set(MYOMP "-fopenmp")  # -lpthread")
	add_definitions(${MYOMP})
	add_definitions("-std=c++17")
endif()

#list firstly the precompiled shared aka utils, secondly the tool-specific components, lastly the tool-specific main
add_executable(${MYPROJECTNAME}
	${MYTOOLSRCPATH}/allocator.cpp
	${MYTOOLSRCPATH}/lodepng.cpp
	${MYTOOLSRCPATH}/newMicrostructureGenerator.cpp
	${MYTOOLSRCPATH}/microStructureHdl.cpp
	${MYTOOLSRCPATH}/Grains.cpp
	${MYTOOLSRCPATH}/SubGrain.cpp
	${MYTOOLSRCPATH}/myQuaternion.cpp
	${MYTOOLSRCPATH}/Settings.cpp
	${MYTOOLSRCPATH}/mymath.cpp
	${MYTOOLSRCPATH}/random.cpp
	${MYTOOLSRCPATH}/IterativeGrainScheduler.cpp

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

#linking process, the target link libraries command is specific for each tool of the toolbox
target_link_libraries(${MYPROJECTNAME} ${MYOMP} ${MYHDFLINKFLAGS} -lnuma)
