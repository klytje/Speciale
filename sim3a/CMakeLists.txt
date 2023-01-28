cmake_minimum_required(VERSION 2.8)

#Name of project
project (jonas)

# Tell cmake where to look for cmake modules.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Add CPC
add_subdirectory(external/cpc)



#Use modern c++ standards
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")

############################################
##              Dependencies              ##
############################################

# Variable for holding dependent libraries.
set(LIBRARIES "")

#Find the various modules, namely, ROOT, GSL, Boost and Armadillo
find_package (GSL REQUIRED)
if (GSL_FOUND)
  message(STATUS "Found GSL")
  include_directories(${GSL_INCLUDE_DIRS})
  add_definitions("-DHAVE_GSL")
  list(APPEND LIBRARIES ${GSL_LIBRARIES})
endif (GSL_FOUND)

SET(AUSALIB_LIBRARIES "/home/krelle/tools/ausalib/build/lib/libAUSA.so")
find_package(AUSALIB REQUIRED)
if (AUSALIB_FOUND)
  include_directories(${AUSALIB_INCLUDES})
  link_directories(${AUSALIB_LIB_DIR})
  list(APPEND LIBRARIES ${AUSALIB_LIBRARIES})
endif(AUSALIB_FOUND)

find_package(ROOT REQUIRED)
if (ROOT_FOUND)
  include_directories(${ROOT_INCLUDE_DIR})
  link_directories(${ROOT_LIBRARY_DIR}) 
  list(APPEND LIBRARIES ${ROOT_LIBRARIES} "-lMathMore -lSpectrum")
endif(ROOT_FOUND)

find_package(Armadillo REQUIRED)
if (ARMADILLO_FOUND)
  include_directories(${ARMADILLO_INCLUDE_DIRS})
  #link_directories(${ARMADILLO_LIBRARY_DIRS})
  list(APPEND LIBRARIES ${ARMADILLO_LIBRARIES})
endif(ARMADILLO_FOUND)

find_package(LOGFT REQUIRED)
include_directories(${LOGFT_INCLUDES})
link_directories(${LOGFT_LIB_DIR})
list(APPEND LIBRARIES ${LOGFT_LIBRARIES})

# Includes
include_directories(include)
include_directories(external/cpc/include)

# Jlib
file(GLOB_RECURSE HEADERS "include/*")
file(GLOB_RECURSE SOURCE    "src/*.cpp")
add_library(jonas ${HEADERS} ${SOURCE} $<TARGET_OBJECTS:CPC>)

target_link_libraries(jonas ${LIBRARIES})

# Balamuth simulator
file(GLOB_RECURSE BALAMUTH_SOURCE "balamuth/*")
add_executable(balamuth ${BALAMUTH_SOURCE})
target_link_libraries(balamuth jonas ${LIBRARIES})


# Installation target
INSTALL(TARGETS jonas balamuth
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)
install(DIRECTORY include/jonas DESTINATION include)
