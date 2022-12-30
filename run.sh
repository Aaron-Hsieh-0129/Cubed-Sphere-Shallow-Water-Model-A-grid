# This shell script is to execute the whole program, plot all graphs, and make videos 
#!/bin/bash
mkdir bin outputs graphs
make clean
make
cd bin && ./csswm
cd ../src && python plot_h_wind.py
cd ../graphs/h/curvilinear
pwd
ffmpeg -r 3 -i %d.png -pix_fmt yuv420p h.mov
cd ../sphere
ffmpeg -r 3 -i %d.png -pix_fmt yuv420p h.mov
cd ../sphere_cartopy
ffmpeg -r 3 -i %d.png -pix_fmt yuv420p h.mov