cmake_minimum_required(VERSION 3.15.0)
set(project_name "AiFrameworkSim")
project(${project_name})

set(CMAKE_POLICY_VERSION_MINIMUM 3.5)

include(ExternalProject)
include(ProcessorCount)
ProcessorCount(NPROC)
if(NPROC EQUAL 0)
    set(NPROC 4)
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
include("${CMAKE_SOURCE_DIR}/build/conanbuildinfo.cmake")
conan_basic_setup()

option(USE_RAMULATOR "USE_RAMULATOR" ON)
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/extern)

# C++ settings
set(CMAKE_CXX_STANDARD 20)
set(ONNX_ML 1)
set(JSON_BuildTests OFF CACHE INTERNAL "")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/build/bin")
set(LIBRARY_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/build/lib")
add_definitions(-D_GLIBCXX_USE_CXX11_ABI=0)
message("BINARY DIR ${CMAKE_BINARY_DIR}")
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(-fsanitize=address)
    add_link_options(-fsanitize=address)
endif()

###############################################################################
# SuperLU_MT (shared by DRAMsim3 and 3D-ICE)
###############################################################################
ExternalProject_Add(superlu_mt_build
    SOURCE_DIR "${PROJECT_SOURCE_DIR}/extern/superlu_mt-4.0.0"
    BUILD_IN_SOURCE TRUE
    CONFIGURE_COMMAND ""
    BUILD_COMMAND make lib
    INSTALL_COMMAND ""
    BUILD_BYPRODUCTS
        "${PROJECT_SOURCE_DIR}/extern/superlu_mt-4.0.0/lib/libsuperlu_mt_OPENMP.a"
)

add_library(superlu_mt STATIC IMPORTED GLOBAL)
set_target_properties(superlu_mt PROPERTIES
    IMPORTED_LOCATION
        "${PROJECT_SOURCE_DIR}/extern/superlu_mt-4.0.0/lib/libsuperlu_mt_OPENMP.a"
)
add_dependencies(superlu_mt superlu_mt_build)

# Point DRAMsim3's find_library(SUPERLU) to our shared SuperLU_MT
set(SUPERLU
    "${PROJECT_SOURCE_DIR}/extern/superlu_mt-4.0.0/lib/libsuperlu_mt_OPENMP.a"
    CACHE FILEPATH "SuperLU_MT library" FORCE
)

# Create symlink so DRAMsim3 finds SuperLU_MT headers at its hardcoded path
execute_process(
    COMMAND ${CMAKE_COMMAND} -E create_symlink
        "${PROJECT_SOURCE_DIR}/extern/superlu_mt-4.0.0"
        "${PROJECT_SOURCE_DIR}/extern/DRAMsim3/ext/SuperLU_MT_3.1"
)

###############################################################################
# 3D-ICE
###############################################################################
ExternalProject_Add(threed_ice_build
    SOURCE_DIR "${PROJECT_SOURCE_DIR}/extern/3d-ice"
    BUILD_IN_SOURCE TRUE
    CONFIGURE_COMMAND ""
    BUILD_COMMAND make -j${NPROC} lib "CFLAGS=-MD -O3 -fopenmp -fPIC -w -I${PROJECT_SOURCE_DIR}/extern/superlu_mt-4.0.0/SRC"
    INSTALL_COMMAND ""
    DEPENDS superlu_mt_build
    BUILD_BYPRODUCTS
        "${PROJECT_SOURCE_DIR}/extern/3d-ice/lib/libthreed-ice-3.1.0.a"
)

add_library(threed_ice STATIC IMPORTED GLOBAL)
set_target_properties(threed_ice PROPERTIES
    IMPORTED_LOCATION
        "${PROJECT_SOURCE_DIR}/extern/3d-ice/lib/libthreed-ice-3.1.0.a"
)
add_dependencies(threed_ice threed_ice_build)

find_package(OpenMP REQUIRED)
find_library(OPENBLAS_LIB openblas REQUIRED)

###############################################################################
# ONNXim source
###############################################################################
add_subdirectory("${PROJECT_SOURCE_DIR}/src")

###############################################################################
# External libraries
###############################################################################
add_subdirectory("${PROJECT_SOURCE_DIR}/extern/ramulator_custom")

add_subdirectory("${PROJECT_SOURCE_DIR}/extern/ramulator2")
include_directories("${PROJECT_SOURCE_DIR}/extern/ramulator2/src")
include_directories("${PROJECT_SOURCE_DIR}/extern/ramulator2/resources/ndp_wrappers")

add_subdirectory("${PROJECT_SOURCE_DIR}/extern/booksim")

set(THERMAL OFF CACHE BOOL "Enable DRAMsim3 thermal simulation")
add_subdirectory("${PROJECT_SOURCE_DIR}/extern/DRAMsim3")
add_dependencies(dramsim3 superlu_mt_build)
target_include_directories(dramsim3 PUBLIC
    "${PROJECT_SOURCE_DIR}/extern/DRAMsim3/src"
    "${PROJECT_SOURCE_DIR}/extern/DRAMsim3/ext/headers"
)

add_subdirectory("${PROJECT_SOURCE_DIR}/extern/dsent")
target_include_directories(dsent PUBLIC "${PROJECT_SOURCE_DIR}/extern/dsent")

add_subdirectory("${PROJECT_SOURCE_DIR}/extern/protobuf/cmake" EXCLUDE_FROM_ALL)
set_target_properties(libprotoc PROPERTIES FOLDER "external/protobuf")
set_target_properties(protoc PROPERTIES FOLDER "external/protobuf")

add_definitions("-DONNX_NAMESPACE=onnx")
add_subdirectory("${PROJECT_SOURCE_DIR}/extern/onnx" EXCLUDE_FROM_ALL)
set_target_properties(onnx PROPERTIES FOLDER "extern/onnx")
set_target_properties(onnx_proto PROPERTIES FOLDER "extern/onnx")

###############################################################################
# Simulator include directories
###############################################################################
target_include_directories(Simulator PUBLIC
    "${PROJECT_SOURCE_DIR}/extern/DRAMsim3/src"
    "${PROJECT_SOURCE_DIR}/extern/dsent"
    "${PROJECT_SOURCE_DIR}/extern/3d-ice/include"
    "${CMAKE_BINARY_DIR}/extern/onnx"
    ${ONNX_INCLUDE_DIRS}
)
target_include_directories(Simulator_lib PUBLIC
    "${PROJECT_SOURCE_DIR}/extern/DRAMsim3/src"
    "${PROJECT_SOURCE_DIR}/extern/dsent"
    "${PROJECT_SOURCE_DIR}/extern/3d-ice/include"
    "${CMAKE_BINARY_DIR}/extern/onnx"
    ${ONNX_INCLUDE_DIRS}
)

# SuperLU SRC only for ThreeDICEWrapper.cc — prevents 'queue' file conflict
set_source_files_properties(
    "${PROJECT_SOURCE_DIR}/src/ThreeDICEWrapper.cc"
    PROPERTIES COMPILE_OPTIONS
        "-I${PROJECT_SOURCE_DIR}/extern/superlu_mt-4.0.0/SRC"
)

###############################################################################
# Build ordering
###############################################################################
add_dependencies(Simulator     onnx_proto threed_ice_build superlu_mt_build)
add_dependencies(Simulator_lib onnx_proto threed_ice_build superlu_mt_build)

###############################################################################
# Linking
###############################################################################
target_link_libraries(Simulator
    ramulator1 booksim2 ramulator dramsim3 dsent
    threed_ice superlu_mt ${OPENBLAS_LIB} OpenMP::OpenMP_CXX
    ${PROTOBUF_LIB} onnx_proto ${CONAN_LIBS} stdc++fs
)

target_link_libraries(Simulator_lib
    ramulator1 booksim2 ramulator dramsim3 dsent
    threed_ice superlu_mt ${OPENBLAS_LIB} OpenMP::OpenMP_CXX
    ${PROTOBUF_LIB} onnx_proto ${CONAN_LIBS} stdc++fs
)