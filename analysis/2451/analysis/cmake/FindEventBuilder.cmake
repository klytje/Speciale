# Tries to determine the location of EventBuilder.
# It will try to look in the systems include/lib paths
#
# You have two options if you want to have a local build of EventBuilder.
#
# 1) Set these two environment variables
#  EVENT_BUILDER_INC_DIR  Include directory
#  EVENT_BUILDER_LIB_DIR  Directory containing the build library
#
# 2) Set this single environment variable
#  EVENT_BUILDERPATH      Path that should contain libAUSA.a and ausa/AUSA.h directly or inside include

if (EVENT_BUILDER_INCLUDES)
    # Already in cache, be silent
    set (EVENT_BUILDER_FIND_QUIETLY TRUE)
endif (EVENT_BUILDER_INCLUDES)

find_path (EVENT_BUILDER_INCLUDES ausa/event/EventRunner.h
        ENV EVENT_BUILDER_INC_DIR ENV EVENT_BUILDERPATH)

find_library(EVENT_BUILDER_LIBRARIES EventBuilder
        HINTS ENV EVENT_BUILDER_LIB_DIR ENV EVENT_BUILDERPATH)

# Handle the QUIETLY and REQUIRED arguments and set EVENT_BUILDER_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (EventBuilder DEFAULT_MSG EVENT_BUILDER_LIBRARIES EVENT_BUILDER_INCLUDES)