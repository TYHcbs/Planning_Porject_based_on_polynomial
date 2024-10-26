##Install
install system dependences
`sudo apt-get install cmake libopenblas-dev liblapack-dev libarpack-dev libarpack2-dev libsuperlu-dev`
install Armadillo
#if not unzip yet
`xz -d armadillo-9.870.2.tar.xz
tar -xvf armadillo-9.870.2.tar`

`cd armadillo-9.870.2`
# if "build" directory already exits
`rm -rf build`

`mkdir build
cd build
cmake ..
make
sudo make install`
