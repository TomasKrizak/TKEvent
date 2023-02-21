#!/bin/bash  

	#echo "                                          "
	#echo "Please enter the full path to *TKEvent* root folder:"
	#read TKEVENT_P
	#echo " "
	TKEVENT_P=../TKEvent/

	rm -rf build
	mkdir build

	cd build
	
	cmake -DTKEvent_PATH=$TKEVENT_P ..

	make
	
	mkdir ../../runs
