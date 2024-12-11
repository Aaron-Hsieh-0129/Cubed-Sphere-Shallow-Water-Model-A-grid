cmake_minimum_required(VERSION 3.0.0)

# Set default build type to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

project(csswm VERSION 0.1.0)

# Set compiler
set(CMAKE_C_COMPILER "/data/Aaron/gcc13/bin/gcc")
set(CMAKE_CXX_COMPILER "/data/Aaron/gcc13/bin/g++")

# Set compiler flags
set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# Enable testing and export compile commands
include(CTest)
enable_testing()
set(CMAKE_EXPORT_COMPILE_COMMANDS True)

# Add local include directories
include_directories(
    /data/Aaron/local/include
    ./include
)

# Add an option for PETSc support
option(USE_PETSC "Enable PETSc support" OFF)

# Include PETSc if USE_PETSC is ON
if(USE_PETSC)
    message(STATUS "PETSc support enabled.")
    find_package(PETSc REQUIRED)

    include_directories(${PETSC_INCLUDES})
    set(LIBS ${LIBS} ${PETSC_LIBRARIES})
else()
    message(STATUS "PETSc support disabled.")
endif()

# Source files
file(GLOB_RECURSE SRC src/*.cpp src/*.hpp)

# Add the executable target
add_executable(csswm ${SRC})

# Link libraries
find_library(libncxxPath netcdf_c++4 "/data/Aaron/local/lib")
target_link_libraries(csswm ${libncxxPath} ${LIBS})

# CPack settings
set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
