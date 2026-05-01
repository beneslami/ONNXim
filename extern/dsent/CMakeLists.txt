cmake_minimum_required(VERSION 3.5)
project(DSENT)

file(GLOB_RECURSE SRCS "*.h" "*.cc")

# Remove Python interface - not needed for ONNXim integration
list(FILTER SRCS EXCLUDE REGEX ".*interface\\.cc$")

add_library(dsent STATIC ${SRCS})

include_directories("${CMAKE_SOURCE_DIR}")

set_target_properties(dsent PROPERTIES
    COMPILE_FLAGS "-std=c++11 -O2"
)
