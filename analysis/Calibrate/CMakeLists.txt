cmake_minimum_required(VERSION 2.8)
project(Calibrate)

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_SOURCE_DIR}/cmake")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
SET (CMAKE_CXX_COMPILER "/usr/bin/clang++")

#----------------------------------------------------------------------------
# Dependencies
#

# Variable for holding dependant libraries
set(LIBRARIES "")

# Include AUSAlib. Look in cmake/FindAUSALIB.cmake for how the specify location
find_package(AUSALIB REQUIRED)
include_directories(${AUSALIB_INCLUDES})
list(APPEND LIBRARIES ${AUSALIB_LIBRARIES})

# Include ROOT
find_package(ROOT REQUIRED)
include_directories(${ROOT_INCLUDE_DIR})
link_directories(${ROOT_LIBRARY_DIR})
list(APPEND LIBRARIES ${ROOT_LIBRARIES} Spectrum MathMore)

link_libraries(${LIBRARIES})

#add_executable(data_fix data_fix.cpp)
add_executable(calibrate calibrate.cpp)
add_executable(calibrate_constructed calibrate_constructed.cpp)
add_executable(mul_compare mul_compare.cpp)
