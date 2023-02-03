cmake_minimum_required(VERSION 2.8)
project(EventAnalysis)

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_SOURCE_DIR}/cmake")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
SET (CMAKE_CXX_COMPILER "/usr/bin/g++")

find_package(AusaEventBuilder REQUIRED)
find_package(AUSALIB REQUIRED)
find_package(GSL REQUIRED)
find_package(ROOT REQUIRED)
find_package(CURL REQUIRED)
link_libraries(AUSA::AUSA AUSA::EventBuilder ROOT::Core ROOT::ROOTDataFrame CURL::libcurl)

add_executable(EventAnalysis main.cpp StandardAnalyzer.cpp)
add_executable(tdc_calibration_true tdc_calibration_true.cpp)
add_executable(tdc_calibration_reconstructed tdc_calibration_reconstructed.cpp)
