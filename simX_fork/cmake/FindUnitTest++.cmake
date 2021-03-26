#
# Find the UnitTest++ includes and library
#
# This module defines
# UNITTEST_INCLUDE_DIR, where to find tiff.h, etc.
# UNITTEST_LIBRARIES, the libraries to link against to use UnitTest++.
# UNITTEST_FOUND, If false, do not try to use UnitTest++.

# also defined, but not for general use are
# UNITTEST_LIBRARY, where to find the UnitTest++ library.
# UNITTEST_DEBUG_LIBRARY, where to find the UnitTest++ library in debug mode.

FIND_PATH(UNITTEST_INCLUDE_DIR unittest++/UnitTest++.h
  /usr/local/include
  /usr/include
)

FIND_PATH(UNITTEST_ALT_INCLUDE_DIR UnitTest++/UnitTest++.h
        /usr/local/include
        /usr/include
)

# With Win32, important to have both
IF(WIN32)
  FIND_LIBRARY(UNITTEST_LIBRARY UnitTest++
               ${UNITTEST_INCLUDE_DIR}/../lib
               /usr/local/lib
               /usr/lib)
  FIND_LIBRARY(UNITTEST_DEBUG_LIBRARY UnitTest++d
               ${UNITTEST_INCLUDE_DIR}/../lib
               /usr/local/lib
               /usr/lib)
ELSE(WIN32)
  # On unix system, debug and release have the same name
  FIND_LIBRARY(UNITTEST_LIBRARY UnitTest++
               ${UNITTEST_INCLUDE_DIR}/../lib
               /usr/local/lib
               /usr/lib)
  FIND_LIBRARY(UNITTEST_DEBUG_LIBRARY UnitTest++
               ${UNITTEST_INCLUDE_DIR}/../lib
               /usr/local/lib
               /usr/lib)
ENDIF(WIN32)

IF(UNITTEST_ALT_INCLUDE_DIR)
  SET(UNITTEST_INCLUDE_DIR ${UNITTEST_ALT_INCLUDE_DIR})
ENDIF(UNITTEST_ALT_INCLUDE_DIR)

IF(UNITTEST_INCLUDE_DIR)
  IF(UNITTEST_LIBRARY)
    SET(UNITTEST_FOUND "YES")
    SET(UNITTEST_LIBRARIES ${UNITTEST_LIBRARY} ${CMAKE_DL_LIBS})
    SET(UNITTEST_DEBUG_LIBRARIES ${UNITTEST_DEBUG_LIBRARY} ${CMAKE_DL_LIBS})
  ENDIF(UNITTEST_LIBRARY)
ENDIF(UNITTEST_INCLUDE_DIR)

IF(UNITTEST_ALT_INCLUDE_DIR)
  file(GLOB_RECURSE unitfiles
          "${UNITTEST_INCLUDE_DIR}/UnitTest++/*.h")
  foreach(unitfile ${unitfiles})
    GET_FILENAME_COMPONENT(unitname ${unitfile} NAME)
    string(TOUPPER ${unitname} UNITFILE)
    string(REPLACE "+" "P" UNITFILE ${UNITFILE})
    string(REPLACE "." "_" UNITFILE ${UNITFILE})
    message(${UNITFILE})
    FILE(WRITE ${CMAKE_CURRENT_BINARY_DIR}/unittest++/${unitname}
"#ifndef AUSALIB_${UNITFILE}
#define AUSALIB_${UNITFILE}

#include <${unitfile}>


#endif //AUSALIB_${UNITFILE}
")
  endforeach()
ENDIF(UNITTEST_ALT_INCLUDE_DIR)