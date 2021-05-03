find_path (AUSADRAW_INCLUDES ausadraw/SetupRender.h
		HINTS ENV AUSADRAW_INC_DIR)
find_library(AUSADRAW_LIBRARIES DrawSetup
		HINTS ENV AUSADRAW_LIB_DIR)
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (AUSADRAW DEFAULT_MSG AUSADRAW_LIBRARIES AUSADRAW_INCLUDES)
