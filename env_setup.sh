#!/bin/bash

sudo apt update
sudo apt upgrade
sudo apt install cmake g++ m4 openssl libssl-dev unzip

sudo unzip envir.zip -d ~
sudo unzip data.zip -d .


# eigen3:
cd ~/envir/eigen-3.3.3/
sudo mkdir build
cd build
sudo cmake ..
sudo make install
sudo cp -r /usr/local/include/eigen3 /usr/include/

# gmp:
cd ~/envir/gmp-6.1.2/
sudo ./configure --prefix=/usr/local
sudo make
sudo make check
sudo make install

# relic:
cd ~/envir/relic-relic-toolkit-0.4.0
sudo mkdir build
cd build
sudo cmake ..
sudo make
sudo make install

# emp-tool, emp-ot, emp-sh2pc:
cd ~/envir/emp-tool-b0292c370c6a4c4acc3f4788b178814cc1898937
sudo mkdir build
cd build
sudo cmake ..
sudo make
sudo make install

cd ~/envir/emp-ot-72cdb09d29af731285646001c046b3b2fba2d351
sudo mkdir build
cd build
sudo cmake ..
sudo make
sudo make install

cd ~/envir/emp-sh2pc-d422f85a39b6a59eb8a8c4db2950f181b84b8912
sudo mkdir build
cd build
sudo cmake ..
sudo make
sudo make install
