cmake_minimum_required(VERSION 3.0.0)
project(hagan VERSION 0.1.0)

# Specified compilers to make sure there is no mistake
set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER "g++")
# include(CTest)
# enable_testing()

# Wildcard to include all the source files in the directory and the sub-ones.
file(GLOB autopilot_sources ./**.cpp ./**.c) 
ADD_EXECUTABLE(hagan ${autopilot_sources})
# Use when you are linking libraries in your source file's directory
# ADD_LIBRARY ( autopilot  ${autopilot_sources})
# TARGET_LINK_LIBRARIES ( autopilot)  

# add_executable(hagan "*.cpp")
# Find libraries to include and link
find_package(Boost REQUIRED) # find Boost
include_directories(${Boost_INCLUDE_DIRS})
# message(STATUS "Boost include dir: ${Boost_INCLUDE_DIRS}")
find_library(QUANTLIB_LIBRARY QuantLIB) # find a path of quantlib library
target_link_libraries(hagan PUBLIC ${QUANTLIB_LIBRARY})
# Add Eigen
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/eigen/) # Source codes are there

# Cmake flags
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++17") # faster runtime
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g3 -std=c++17") # faster compilation for debugging
set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
# include(CPack)