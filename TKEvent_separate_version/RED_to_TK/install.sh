#!/bin/bash  

	echo "                                          "
	echo "Please enter the full path to *Event3D* root folder:"
	read EVENT3DP
	echo " "

## /pbs/home/m/mmacko/Event3D/

	rm -rf build
	mkdir build

	cd build
	
	cmake -DEvent3D_PATH=$EVENT3DP ..

	make
