cmake_minimum_required(VERSION 3.0.0)
project(qtest VERSION 0.1.0)

# Specified compilers to make sure there is no mistake
set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER "g++")

include(CTest)
enable_testing()

# Header files to include
# If you want to manually specify the path to Boost, you can use the followings
# include_directories(/usr/local/include/)
# include_directories(/usr/local/include/boost/)

add_executable(qtest main.cpp)
find_package(Boost REQUIRED) # find Boost
include_directories(${Boost_INCLUDE_DIRS})
message(STATUS "Boost include dir: ${Boost_INCLUDE_DIRS}")

find_library(QUANTLIB_LIBRARY QuantLIB) # find a path of quantlib library
target_link_libraries(qtest PUBLIC ${QUANTLIB_LIBRARY})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++17")

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
