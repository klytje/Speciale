# Tries to determine the location of CPC-library (Computer Physics Communications).
# It will try to look in the systems include/lib paths and define
#
#  CPC_FOUND          CPC is installed
#  CPC_INCLUDE_DIRS   Include directory
#  CPC_LIBRARIES      Path to library
#  CPC_LIB_DIRS       Directory of the library

message(STATUS "Looking for Computer Physics Communications library...")

find_path(CPC_INCLUDE_DIR cpc/cwfcomp.h)

find_library(CPC_LIBRARY NAMES libCPC CPC)
find_path(CPC_LIB_DIR lib/libCPC.so)
set(CPC_LIB_DIR ${CPC_LIB_DIR} "/lib")  #Cheating...

# Handle the QUIETLY and REQUIRED arguments and set CPC_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (CPC DEFAULT_MSG CPC_LIBRARY CPC_INCLUDE_DIR)

set(CPC_INCLUDE_DIRS ${CPC_INCLUDE_DIR})
set(CPC_LIBRARIES ${CPC_LIBRARY})
set(CPC_LIB_DIRS ${CPC_LIB_DIR})
