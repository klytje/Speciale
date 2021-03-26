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

# Include ROOT
find_package(ROOT REQUIRED)
include_directories(${ROOT_INCLUDE_DIR})
link_directories(${ROOT_LIBRARY_DIR})
list(APPEND LIBRARIES ${ROOT_LIBRARIES} Spectrum MathMore)

link_libraries(${LIBRARIES})

add_executable(dalitz_raw dalitz_raw.cpp)
add_executable(dalitz_slice dalitz_slice.cpp)
add_executable(dalitz dalitz.cpp)
add_executable(eCM eCM.cpp)
add_executable(ExpT ExpT.cpp)
add_executable(QpT QpT.cpp)
add_executable(fynbo_plot fynbo_plot.cpp)
add_executable(kinematic_lab kinematic_lab.cpp)
add_executable(kinematic kinematic.cpp)
add_executable(angular_projection angular_projection.cpp)
add_executable(radial_projection radial_projection.cpp)
