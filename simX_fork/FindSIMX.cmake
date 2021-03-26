# Tries to determine the location of simX.
# It will try to look in the systems include/lib paths
#
# You have two options if you want to have a local build of AUSAlib.
#
# 1) Set these two environment variables
#  SIMX_INC_DIR  Include directory
#  SIMX_LIB_DIR  Directory containing the build library
#
# 2) Set this single environment variable
#  SIMXPATH      Path that should contain libsimX.a and simX/git.h directly or inside include

if (SIMX_INCLUDES)
    # Already in cache, be silent
    set (SIMX_FIND_QUIETLY TRUE)
endif (SIMX_INCLUDES)

find_path (SIMX_INCLUDES simX/git.h
        HINTS $ENV{SIMXPATH}/include ENV SIMX_INC_DIR ENV SIMXPATH)

find_library(SIMX_LIBRARIES simX
        HINTS ENV SIMX_LIB_DIR ENV SIMXPATH)

# Handle the QUIETLY and REQUIRED arguments and set SIMX_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (SIMX DEFAULT_MSG SIMX_LIBRARIES SIMX_INCLUDES)