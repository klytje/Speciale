cmake_minimum_required(VERSION 2.8)
project(EventAnalysis)

# Tell cmake where to look for cmake modules.
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_SOURCE_DIR}/cmake")

#Lets write modern C++
SET (CMAKE_CXX_COMPILER "/usr/bin/clang++")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")

#----------------------------------------------------------------------------
# Dependencies
#

# Variable for holding dependant libraries
set(LIBRARIES "")

find_package(EventBuilder REQUIRED)
include_directories(${EVENT_BUILDER_INCLUDES})
list(APPEND LIBRARIES ${EVENT_BUILDER_LIBRARIES})

# Include AUSAlib. Look in cmake/FindAUSALIB.cmake for how the specify location
find_package(AUSALIB REQUIRED)
include_directories(${AUSALIB_INCLUDES})
list(APPEND LIBRARIES ${AUSALIB_LIBRARIES})

# Include ROOT
find_package(ROOT REQUIRED)
include_directories(${ROOT_INCLUDE_DIR})
link_directories(${ROOT_LIBRARY_DIR})
list(APPEND LIBRARIES ${ROOT_LIBRARIES} Spectrum MathMore)

find_package(CURL)
if (CURL_FOUND)
    include_directories(${CURL_INCLUDE_DIRS})
    list(APPEND LIBRARIES ${CURL_LIBRARIES})
    add_definitions("-DAUSA_HAS_LIBCURL")
endif(CURL_FOUND)

link_libraries(${LIBRARIES})

add_executable(EventAnalysis main.cpp IFA006Analyzer.cpp MulSpectrum.cpp AntiProtonCutter.cpp)
#add_executable(SinglesAnalysis SinglesAnalysis.cpp)
#add_executable(SolidCalc SolidCalc.cpp)
#add_executable(SDSolidCalc SDSolidCalc.cpp)
#add_executable(SingleCrossSection SingleCrossSection.cpp)
#add_executable(DtAnalysis DtAnalysis.cpp)
#add_executable(ThicknessAnalysis ThicknessAnalysis.cpp)
#add_executable(TdcEffeciency TdcEffeciency.cpp)
#add_executable(TdcCum TdcCum.cpp)
#add_executable(PhysAnalysis PhysAnalysis.cpp)
#add_executable(DalitzEff DalitzEff.cpp)
