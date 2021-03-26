cmake_minimum_required(VERSION 2.8.4)
project(simX)

#set(CMAKE_CXX_COMPILER "/usr/bin/clang")

# Tell cmake where to look for cmake modules.
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_SOURCE_DIR}/cmake")

############################################
##                 git                    ##
############################################
# Get git information
include(GetGitRevisionDescription REQUIRED)
get_git_head_revision(GIT_BRANCH GIT_SHA1)
string(SUBSTRING ${GIT_SHA1} 0 8 GIT_SHORT_SHA1)

# Write information into cpp file.
set(GIT_FILE "${CMAKE_CURRENT_BINARY_DIR}/git.cpp")
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/source/git.cpp.in" ${GIT_FILE} @ONLY)

############################################
##          Define simX target            ##
############################################
include_directories(include)
include_directories(external/coulomb)
include_directories(external/include)
include_directories(external/include/gte)

# Find files to compile
include(GlobExclude)

file(GLOB_RECURSE HEADERS "include/*")
GLOB_EXCLUDE(SOURCE "source/*.cpp" "/parser/") # All source files except parser
file(GLOB_RECURSE SOURCE_PARSER "source/parser/*") # Parser source files


add_subdirectory(external/coulomb)
add_library(simX_parser OBJECT ${SOURCE_PARSER}) # Make object library of parsers
add_library(simX ${HEADERS}
  			${SOURCE}
  			$<TARGET_OBJECTS:simX_parser>
			$<TARGET_OBJECTS:Coulomb>
			${GIT_FILE}
		)

add_executable(simX_runner runner/main.cpp)
set_target_properties(simX_runner PROPERTIES OUTPUT_NAME simX)

############################################
##              Dependencies              ##
############################################

# Variable for holding dependant libraries
set(LIBRARIES "")

# Include AUSAlib
find_package(AUSALIB REQUIRED)
if (AUSALIB_FOUND)
  include_directories(${AUSALIB_INCLUDES})
  list(APPEND LIBRARIES ${AUSALIB_LIBRARIES})
endif (AUSALIB_FOUND)

# Include ROOT
find_package(ROOT REQUIRED)
if (ROOT_FOUND)
  include_directories(${ROOT_INCLUDE_DIR})
  link_directories(${ROOT_LIBRARY_DIR})
  list(APPEND LIBRARIES ${ROOT_LIBRARIES} Spectrum MathMore)
endif(ROOT_FOUND)

# Include GSL
find_package(GSL)
if (GSL_FOUND)
  include_directories(${GSL_INCLUDE_DIRS})
  add_definitions("-DHAVE_GSL")
  list(APPEND LIBRARIES ${GSL_LIBRARIES})
endif (GSL_FOUND)

#find_package(Boost REQUIRED)
#if (Boost_FOUND)
#  include_directories(${Boost_INCLUDE_DIRS})
#endif (Boost_FOUND)

include_directories("boost_1_58_0")

find_package(CURL)
if (CURL_FOUND)
    include_directories(${CURL_INCLUDE_DIRS})
    list(APPEND LIBRARIES ${CURL_LIBRARIES})
    add_definitions("-DAUSA_HAS_LIBCURL")
endif(CURL_FOUND)


# Add libraries to simX target
target_link_libraries(simX ${LIBRARIES})
target_link_libraries(simX_runner simX ${LIBRARIES})

############################################
##          Instalation target            ##
############################################

# Installation target
INSTALL(TARGETS simX
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

INSTALL(TARGETS simX_runner
        RUNTIME DESTINATION bin
)

install(DIRECTORY include/simX    DESTINATION include)
install(FILES     include/tree.hh DESTINATION include)

############################################
##                Doxygen                 ##
############################################
find_package(Doxygen)
if(DOXYGEN_FOUND)
  set(simX_Doxygen "${CMAKE_CURRENT_BINARY_DIR}/simX.dox")
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/simX.dox.in ${simX_Doxygen} @ONLY)
  add_custom_target(doc
    ${DOXYGEN_EXECUTABLE} ${simX_Doxygen}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen" VERBATIM
)
endif(DOXYGEN_FOUND)
############################################
##         Check compiler version         ##
############################################

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    # require at least gnu 4.7 - earlier version do not have c++11
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.7)
        message(FATAL_ERROR "G++ version must be at least 4.7!")
    endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    # require at least clang 3.5 - earlier version crashes
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.5)
        message(FATAL_ERROR "Clang version must be at least 3.5!")
    endif()
else()
    message(WARNING "You are using an unsupported compiler! Compilation has only been tested with Clang and GCC.")
endif()

############################################
##             Compiler flags             ##
############################################

# We use C++11
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")

if(MSVC)
  # Force to always compile with W4
  if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
    string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
  endif()
elseif(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
  # Update if necessary
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-sign-compare -Wno-reorder")
endif()

# Pick a nice default build type
if (NOT CMAKE_BUILD_TYPE)
	message(STATUS "No build type selected, default to RelWithDebInfo")
	set(CMAKE_BUILD_TYPE "RelWithDebInfo")
endif()

# Use precompiled headers:
option(PCH "If true enable pre compiled headers" OFF)
if (PCH)
	MESSAGE(STATUS "Using precompiled headers (PCH=${PCH})")
    include(cotire)
    set_target_properties(simX PROPERTIES COTIRE_ADD_UNITY_BUILD FALSE)
    cotire(simX)

    set_target_properties(simX_parser
      PROPERTIES COTIRE_CXX_PREFIX_HEADER_INIT "internal/parser_headers.h")

    cotire(simX_parser)
else()
	MESSAGE(STATUS "Not using precompiled headers (PCH=${PCH})")
endif()

############################################
##             Unit Testing               ##
############################################
find_package(UnitTest++)

if(UNITTEST_FOUND)
  message(STATUS "Compiling unit test")
  file(GLOB_RECURSE TEST_SRC
          "test/*.h"
          "test/*.cpp")
  add_executable(simXtest ${TEST_SRC})
  target_include_directories(simXtest PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
  target_link_libraries(simXtest debug ${UNITTEST_DEBUG_LIBRARIES} general ${UNITTEST_LIBRARIES})
  target_link_libraries(simXtest simX ${AUSA} ${LIBRARIES})
else(UNITTEST_FOUND)
  message(STATUS "Not compilin unit tests. Add libunittest++-dev to your system if you want this")
endif(UNITTEST_FOUND)

############################################
##             Code coverage              ##
############################################


IF (CMAKE_BUILD_TYPE STREQUAL "Coverage")
  include(CodeCoverage)
  SETUP_TARGET_FOR_COVERAGE(coverage_report simXtest coverage)
ENDIF()
