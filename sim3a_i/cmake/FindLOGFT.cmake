find_path (LOGFT_INCLUDES logft/logft.h
		HINTS ENV LOGFT_INC_DIR)
find_library(LOGFT_LIBRARIES logft
		HINTS ENV LOGFT_LIB_DIR)
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (LOGFT DEFAULT_MSG LOGFT_LIBRARIES LOGFT_INCLUDES)
