# Tries to determine the location of AUSAlib.
# It will try to look in the systems include/lib paths
#
# You have two options if you want to have a local build of AUSAlib.
#
# 1) Set these two environment variables
#  AUSALIB_INC_DIR  Include directory
#  AUSALIB_LIB_DIR  Directory containing the build library
#
# 2) Set this single environment variable
#  AUSALIBPATH      Path that should contain libAUSA.a and ausa/AUSA.h directly or inside include

if (AUSALIB_INCLUDES)
    # Already in cache, be silent
    set (AUSALIB_FIND_QUIETLY TRUE)
endif (AUSALIB_INCLUDES)

find_path (AUSALIB_INCLUDES ausa/AUSA.h
        HINTS $ENV{AUSALIBPATH}/include ENV AUSALIB_INC_DIR ENV AUSALIBPATH)

find_library(AUSALIB_LIBRARIES AUSA
        HINTS ENV AUSALIB_LIB_DIR ENV AUSALIBPATH)

# Handle the QUIETLY and REQUIRED arguments and set AUSALIB_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (AUSALIB DEFAULT_MSG AUSALIB_LIBRARIES AUSALIB_INCLUDES)