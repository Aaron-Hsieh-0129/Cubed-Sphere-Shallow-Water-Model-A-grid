# This shell script is to execute the whole program, plot all graphs, and make videos 
#!/bin/bash
rm -rf build outputs
mkdir build
cd build/ && cmake ../ && make && ./csswm

conda activate python3.8
cd ../src && python plot_h_wind.py && python plot_zeta.py

cd ../graphs/h/curvilinear
ffmpeg -r 3 -i %d.png -pix_fmt yuv420p h+wind.mov -y
ffmpeg -r 3 -i %d.png h+wind.gif -y

cd ../sphere
ffmpeg -r 3 -i %d.png -pix_fmt yuv420p h+wind.mov -y
ffmpeg -r 3 -i %d.png h+wind.gif -y

cd ../sphere_cartopy
ffmpeg -r 3 -i %d.png -pix_fmt yuv420p h+wind.mov -y
ffmpeg -r 3 -i %d.png h+wind.gif -y

cd ../../zeta/curvilinear
ffmpeg -r 3 -i %d.png -pix_fmt yuv420p zeta.mov -y
ffmpeg -r 3 -i %d.png zeta.gif -y

cd ../sphere_cartopy
ffmpeg -r 3 -i %d.png -pix_fmt yuv420p zeta.mov -y
ffmpeg -r 3 -i %d.png zeta.gif -y
conda deactivate