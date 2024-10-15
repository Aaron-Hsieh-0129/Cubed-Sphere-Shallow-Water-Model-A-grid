Cubed-Sphere Shallow Water Model (A grid)
=========================================

This is a cubed-sphere shallow water model


Prerequisite
------------

- C++ compiler (higher than C++11)
- CMake (higher than 3.0.0) (You can create your own Makefile by translating the CMakefile.txt if you don't want to use CMake)
- netcdf-cxx4 (hdf5, netcdf-c are needed for netcdf-cxx) [optional]


How to Use
----------

1. Clone the project using:

   .. code-block:: bash

      git clone https://github.com/Aaron-Hsieh-0129/Cubed-Sphere-Shallow-Water-Model-A-grid.git

2. You can change the model settings by changing the macro flags in the `./src/define.hpp` and also the configuration object in the `./src/main.cpp`. In general, you only need to modify the configuration in `./src/main.cpp`. But if you want to change the numerical methods or the running cases, you might need to modify the flags in the `./src/Config.hpp`.

3. You are able to run the model by running the command under the project folder:

   .. code-block:: bash

      sh run.sh

   or you can use your own command by referencing the command in `run.sh`.

Optional for NetCDF output
-------------------------------------------

1. Install netcdf-cxx

   It's a little bit complicated to install libraries for C/C++. I will provide a tutorial for installing C/C++ compiler and the libraries in another file. Here, you don't need to have sudo privilege to install anything.

2. Link the installed libraries.

- If you don't need to use PETSc and netcdf, you might need to turn off the link command in CMakeLists.txt.
- You need to change the libraries path (netcdf, petsc) to your own path.
- Change include path in CMakeLists.txt:

  .. code-block:: cmake

    include_directories(
      </path/to/your/netcdf>/include
    )

- Change library link path:

  .. code-block:: cmake

    find_library(libncxxPath netcdf_c++4 "<path to your netcdf_c++4>/lib")


3. Change the flags in `./src/define.hpp`:

- Turn on `NCOUTPUT` and turn off `TXTOUTPUT` to use netcdf output.
