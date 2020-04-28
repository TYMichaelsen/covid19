# Build guide

It is advisable to use GNU `make` command to build:

Just run `make` at the directory where this README.md file is located. On the first-time build, `make` will build the base image as well. On the subsequent builds, when base image (covid19_1.5.sif) already exists, `make` will go ahead and build other parts. 

In case you make some changes on the base defition file, `make clean` and then `make` will clean up and build from srcatch. This will take more time.  
