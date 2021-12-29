cmake_minimum_required(VERSION 3.0.0)
project(BSMJ VERSION 0.1.0)

# Specified compilers to make sure there is no mistake
set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER "g++")
include(CTest)
enable_testing()

# Wildcard to include all the source files in the directory and the sub-ones.
file(GLOB autopilot_sources ./**.cpp ./**.c) 
ADD_EXECUTABLE(BSMJ ${autopilot_sources})
# Use when you are linking libraries in your source file's directory
# ADD_LIBRARY ( autopilot  ${autopilot_sources})
# TARGET_LINK_LIBRARIES ( autopilot)  

# add_executable(BSMJ "*.cpp")
# Find libraries to include and link
find_package(Boost REQUIRED) # find Boost
include_directories(${Boost_INCLUDE_DIRS})
# message(STATUS "Boost include dir: ${Boost_INCLUDE_DIRS}")
find_library(QUANTLIB_LIBRARY QuantLIB) # find a path of quantlib library
target_link_libraries(BSMJ PUBLIC ${QUANTLIB_LIBRARY})

# Cmake flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++17")
set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)