#!/bin/bash 

	rm -rf ./build ./src/dicts/ ./lib
	mkdir build lib

	cd include

echo "#############################################"
echo "#         GENERATE ROOT DICTIONARIES        #"
echo "#############################################"
echo " "

	echo "rootcling: ./include/Event_3D.hh         -> ./include/Event_3Ddict_rdict.pcm         + ./include/Event_3Ddict.cpp"
      	      rootcling -f ../lib/Event_3Ddict.cpp  Event_3D.hh+
	echo " "
	echo "Dictionaries generated!"
	echo " "

	cd ..

	cp -r ./lib/ ./src/dicts/
	rm -rf ./lib/*.cpp
	rm -rf ./src/dicts/*.pcm

echo "#############################################"
echo "#         MODULE COMPILATION (CMAKE)        #"
echo "#############################################"
echo " "

	cd build
	cmake .. || { echo $'\n****Cmake failed, installation of Event3D aborted!****' ; exit 1; } 

echo " "
echo "#############################################"
echo "#         MODULE COMPILATION (MAKE)         #"
echo "#############################################"
echo " "

	make    || { echo $'\n****Make failed, installation of Event3D aborted!****' ; exit 1; }

	cp libEvent3D.so ../lib/libEvent3D.so
	rm libEvent3D.so
	
