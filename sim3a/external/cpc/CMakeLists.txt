# Find files to compile
include_directories(include)

file(GLOB_RECURSE HEADERS "include/*.h")
file(GLOB_RECURSE SOURCE    "src/*.cpp")
add_library(CPC OBJECT ${HEADERS} ${SOURCE})
