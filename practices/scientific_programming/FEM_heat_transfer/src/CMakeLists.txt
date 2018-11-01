cmake_minimum_required(VERSION 2.8)
PROJECT(2d_Unsteady)

set(VTK_DIR "/usr/lib64/vtk")
set(DCO_DIR "$ENV{HOME}/SiSc2/SiScLabWS16-Gr2/dco_cpp_v3.1.4_trial_lin64_gcc")
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

INCLUDE_DIRECTORIES(. ${DCO_DIR}/include)
LINK_DIRECTORIES(${DCO_DIR}/lib)

file(GLOB SOURCES "*.cpp")

add_executable(2d_Unsteady ${SOURCES})

target_link_libraries(2d_Unsteady ${VTK_LIBRARIES})
target_link_libraries(2d_Unsteady libdcoc.a)

set(CMAKE_CXX_FLAGS "-O2")        ## Optimize
set(CMAKE_EXE_LINKER_FLAGS "-s")  ## Strip binary



