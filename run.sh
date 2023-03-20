# This shell script is to execute the whole program, plot all graphs, and make videos 
#!/bin/bash
rm -rf build outputs
rm -rf /data/Aaron/CSSWM-A-grid/Cases/Barotropic/180/outputs/
mkdir build
cd build/ && cmake ../ && make && ./csswm