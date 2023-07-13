# This shell script is to execute the whole program, plot all graphs, and make videos 
#!/bin/bash
rm -rf build
mkdir build
cd build/ && cmake ../ && make && ./csswm