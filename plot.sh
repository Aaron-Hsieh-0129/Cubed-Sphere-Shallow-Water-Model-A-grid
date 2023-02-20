cd src && python plot_h_wind.py && python plot_zeta.py

cd ../graphs/h/curvilinear
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