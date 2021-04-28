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

add_executable(sim_compare sim_compare.cpp)
add_executable(Dalitz_3D Dalitz_3D.cpp)
add_executable(tdc_plot tdc_plot.cpp)
add_executable(depend_mul depend_mul.cpp)
add_executable(depend_E depend_E.cpp)
add_executable(width_converter width_converter.cpp)
add_executable(depend_width depend_width.cpp)
add_executable(ang_compare ang_compare.cpp)
add_executable(ang_compare_single ang_compare_single.cpp)
add_executable(theta theta.cpp)
add_executable(ang_cor_fit ang_cor_fit.cpp)
add_executable(corr_funcs corr_funcs.cpp)
add_executable(theta_data theta_data.cpp)
