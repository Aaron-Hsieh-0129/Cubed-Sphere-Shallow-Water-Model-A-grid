cmake_minimum_required(VERSION 3.10)

# Set default build type to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build" FORCE)
endif()

project(csswm VERSION 0.1.0)

# Set default prefix directory if not specified
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
    set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install" CACHE PATH "Installation directory prefix" FORCE)
endif()

# Provide option for users to specify custom compiler paths
set(CMAKE_C_COMPILER "$ENV{CC}" CACHE FILEPATH "C compiler")
set(CMAKE_CXX_COMPILER "$ENV{CXX}" CACHE FILEPATH "C++ compiler")

# Set compiler flags
set(CMAKE_CXX_FLAGS "-Wall -Wextra" CACHE STRING "C++ compiler flags")
set(CMAKE_CXX_FLAGS_DEBUG "-g" CACHE STRING "Debug flags")
set(CMAKE_CXX_FLAGS_RELEASE "-O3" CACHE STRING "Release flags")

# Enable testing and export compile commands
include(CTest)
enable_testing()
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE BOOL "Export compile commands" FORCE)

# Include ExternalProject for dependency management
include(ExternalProject)

# Set where dependencies will be downloaded and installed
set(DEP_INSTALL_DIR "${CMAKE_SOURCE_DIR}/_deps/install")
file(MAKE_DIRECTORY ${DEP_INSTALL_DIR})

# Get number of available cores
include(ProcessorCount)
ProcessorCount(8)
if(N_CORES GREATER 1)
    set(CMAKE_BUILD_PARALLEL_LEVEL ${N_CORES})
    message(STATUS "Using ${N_CORES} cores for building")
else()
    message(STATUS "Using single core for building")
endif()

# Local include directories
include_directories(
    ${CMAKE_SOURCE_DIR}/include
    ${DEP_INSTALL_DIR}/include
)

# HDF5
find_library(HDF5_LIB NAMES hdf5 hdf5_serial HINTS ${DEP_INSTALL_DIR}/lib NO_DEFAULT_PATH)
if(NOT HDF5_LIB)
    message(STATUS "HDF5 not found in ${DEP_INSTALL_DIR}, installing...")
    ExternalProject_Add(
        hdf5_ext
        URL            https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.14/src/hdf5-1.8.14.tar.gz
        PREFIX         ${CMAKE_SOURCE_DIR}/_deps/hdf5
        CONFIGURE_COMMAND ${CMAKE_COMMAND} -E env
            CC=${CMAKE_C_COMPILER}
            CXX=${CMAKE_CXX_COMPILER}
            <SOURCE_DIR>/configure
            --prefix=${DEP_INSTALL_DIR}
            --enable-cxx
        BUILD_COMMAND  make -j${N_CORES}
        INSTALL_COMMAND make install
        UPDATE_COMMAND ""
        BUILD_IN_SOURCE 1  # Ensures build happens in source dir, like your manual test
        BUILD_ALWAYS   FALSE
    )
    set(HDF5_DIR ${DEP_INSTALL_DIR})
    set(HDF5_NEEDS_INSTALL TRUE)
    # Set expected library path after install
    set(HDF5_LIB ${DEP_INSTALL_DIR}/lib/libhdf5.so)
else()
    message(STATUS "HDF5 found in ${DEP_INSTALL_DIR}/lib")
endif()

# NetCDF-C (depends on HDF5)
find_library(NETCDF_LIB NAMES netcdf HINTS ${DEP_INSTALL_DIR}/lib NO_DEFAULT_PATH)
if(NOT NETCDF_LIB)
    message(STATUS "NetCDF-C not found in ${DEP_INSTALL_DIR}, installing...")
    ExternalProject_Add(
        netcdf-c_ext
        GIT_REPOSITORY https://github.com/Unidata/netcdf-c.git
        GIT_TAG        v4.3.3.1
        GIT_SHALLOW    TRUE
        PREFIX         ${CMAKE_SOURCE_DIR}/_deps/netcdf-c
        CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=${DEP_INSTALL_DIR} --with-hdf5=${DEP_INSTALL_DIR} CPPFLAGS=-I${DEP_INSTALL_DIR}/include LDFLAGS=-L${DEP_INSTALL_DIR}/lib
        BUILD_COMMAND  make -j${N_CORES}
        INSTALL_COMMAND make install
        DEPENDS        hdf5_ext
        UPDATE_COMMAND ""
        BUILD_ALWAYS   FALSE
    )
    set(NETCDF_C_DIR ${DEP_INSTALL_DIR})
    set(NETCDF_C_NEEDS_INSTALL TRUE)
    # Set expected library path after install
    set(NETCDF_LIB ${DEP_INSTALL_DIR}/lib/libnetcdf.so)
else()
    message(STATUS "NetCDF-C found in ${DEP_INSTALL_DIR}/lib")
endif()

# NetCDF-C++4 (depends on netcdf-c and hdf5)
find_library(NETCDF_CXX4_LIB NAMES netcdf_c++4 HINTS ${DEP_INSTALL_DIR}/lib NO_DEFAULT_PATH)
if(NOT NETCDF_CXX4_LIB)
    message(STATUS "NetCDF-C++4 not found in ${DEP_INSTALL_DIR}, installing...")
    set(NETCDF_CXX4_DEPENDS "")
    if(HDF5_NEEDS_INSTALL)
        list(APPEND NETCDF_CXX4_DEPENDS hdf5_ext)
    endif()
    if(NETCDF_C_NEEDS_INSTALL)
        list(APPEND NETCDF_CXX4_DEPENDS netcdf-c_ext)
    endif()
    ExternalProject_Add(
        netcdf-cxx4_ext
        GIT_REPOSITORY https://github.com/Unidata/netcdf-cxx4.git
        GIT_TAG        v4.2.1
        GIT_SHALLOW    TRUE
        PREFIX         ${CMAKE_SOURCE_DIR}/_deps/netcdf-cxx4
        CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=${DEP_INSTALL_DIR} --with-netcdf=${DEP_INSTALL_DIR} --with-hdf5=${DEP_INSTALL_DIR} CPPFLAGS=-I${DEP_INSTALL_DIR}/include LDFLAGS=-L${DEP_INSTALL_DIR}/lib
        BUILD_COMMAND  make -j${N_CORES}
        INSTALL_COMMAND make install
        DEPENDS        ${NETCDF_CXX4_DEPENDS}
        UPDATE_COMMAND ""
        BUILD_ALWAYS   FALSE
    )
    set(NETCDF_CXX4_NEEDS_INSTALL TRUE)
    # Set expected library path after install
    set(NETCDF_CXX4_LIB ${DEP_INSTALL_DIR}/lib/libnetcdf_c++4.so)
else()
    message(STATUS "NetCDF-C++4 found in ${DEP_INSTALL_DIR}/lib")
endif()

# Source files
file(GLOB_RECURSE SRC
    "${CMAKE_SOURCE_DIR}/src/*.cpp"
    "${CMAKE_SOURCE_DIR}/src/*.hpp"
)

# Add executable
add_executable(${PROJECT_NAME} ${SRC})

# Link libraries (no need to check NOTFOUND since we set paths)
target_link_libraries(${PROJECT_NAME} PRIVATE
    ${HDF5_LIB}
    ${NETCDF_LIB}
    ${NETCDF_CXX4_LIB}
    ${LIBS}
)

# Add dependencies if needed
if(HDF5_NEEDS_INSTALL)
    add_dependencies(${PROJECT_NAME} hdf5_ext)
endif()
if(NETCDF_C_NEEDS_INSTALL)
    add_dependencies(${PROJECT_NAME} netcdf-c_ext)
endif()
if(NETCDF_CXX4_NEEDS_INSTALL)
    add_dependencies(${PROJECT_NAME} netcdf-cxx4_ext)
endif()

# Installation rules
install(TARGETS ${PROJECT_NAME}
    RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
)

# CPack settings
set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

# Print configuration summary
message(STATUS "Configuration Summary:")
message(STATUS "  Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "  Install prefix: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "  C Compiler: ${CMAKE_C_COMPILER}")
message(STATUS "  C++ Compiler: ${CMAKE_CXX_COMPILER}")
