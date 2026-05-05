CMAKE_MINIMUM_REQUIRED(VERSION 3.15)

set(LIB_NAME "Simulator")

file(GLOB_RECURSE SRC_FILES
  "${CMAKE_SOURCE_DIR}/src/*.h"
  "${CMAKE_SOURCE_DIR}/src/*.cc"
)

add_executable(${LIB_NAME} ${SRC_FILES})
add_library(${LIB_NAME}_lib ${SRC_FILES})

# SuperLU SRC only for ThreeDICEWrapper.cc
# Keeps 'queue' file from shadowing C++ standard <queue> everywhere else
set_source_files_properties(
    "${CMAKE_SOURCE_DIR}/src/ThreeDICEWrapper.cc"
    PROPERTIES COMPILE_OPTIONS
        "-I${CMAKE_SOURCE_DIR}/extern/superlu_mt-4.0.0/SRC"
)