# This shell script is to execute the whole program, plot all graphs, and make videos 
#!/bin/bash
rm -rf build outputs
mkdir build
cd build/ && cmake ../ && make && ./csswm

cd ../src && python plot_h_wind.py && python plot_zeta.py

cd ../graphs/h/curvilinear
ffmpeg -r 3 -i %d.png -pix_fmt yuv420p curvilinear.mov -y
cd ../sphere
ffmpeg -r 3 -i %d.png -pix_fmt yuv420p sphere.mov -y
cd ../sphere_cartopy
ffmpeg -r 3 -i %d.png -pix_fmt yuv420p sphere_cartopy.mov -y

cd ../../zeta
ffmpeg -r 3 -i %d.png -pix_fmt yuv420p zeta.mov -y