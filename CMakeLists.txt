cmake_minimum_required(VERSION 3.0.0)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(csswm VERSION 0.1.0)

include(CTest)
enable_testing()

include_directories(
    include/netcdf
)

file(GLOB_RECURSE SRC src/*.cpp src/*.hpp)
add_executable(csswm ${SRC})

find_library(libncPath netcdf lib/netcdf)
find_library(libncxxPath netcdf_c++4 lib/netcdf)
target_link_libraries(csswm ${libncPath} ${libncxxPath})

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
