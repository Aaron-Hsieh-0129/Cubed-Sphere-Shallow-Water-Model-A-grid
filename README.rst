Cubed-Sphere Shallow Water Model (A grid)
=========================================

This is a cubed-sphere shallow water model


Prerequisite
------------

- C++ compiler (higher than C++11)
- CMake (higher than 3.10.0)
- hdf5 1.8.14 (can be automatically downloaded through CMake)
- netcdf-c 4.3.3.1 (can be automatically downloaded through CMake)
- netcdf-cxx4 4.2.1 (can be automatically downloaded through CMake)


How to Use
----------

1. Clone the project using:

   .. code-block:: bash

      git clone https://github.com/Aaron-Hsieh-0129/Cubed-Sphere-Shallow-Water-Model-A-grid.git

2. In `csswm_config.txt`, some common configurations can be set inside. For different test cases, modify the flags inside `src/deine.hpp`. 

3. Compile the project using
   
   .. code-block:: bash

      mkdir build && cd build
      cmake ..
      # If you want to specify your compiler
      cmake .. -DCMAKE_C_COMPILER=/your/compiler/path -DCMAKE_CXX_COMPILER=/your/compiler/path
      make

   If the dependencies are not detected, they will be installed in `project root/_deps`. 

4. You are able to run the model by running the command under the project folder:

   .. code-block:: bash

      sh run.sh

   or you can use your own command by referencing the command in `run.sh`.


