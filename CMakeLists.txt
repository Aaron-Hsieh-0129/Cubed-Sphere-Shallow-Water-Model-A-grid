cmake_minimum_required(VERSION 3.0.0)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

project(csswm VERSION 0.1.0)

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(CTest)
enable_testing()

set(CMAKE_EXPORT_COMPILE_COMMANDS True)

include_directories(
    /home/Aaron/local/include
)

file(GLOB_RECURSE SRC src/*.cpp src/*.hpp)
add_executable(csswm ${SRC})

find_library(libncxxPath netcdf_c++4 "/home/Aaron/local/lib")
target_link_libraries(csswm ${libncxxPath})

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
