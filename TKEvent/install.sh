#!/bin/bash 

	rm -rf ./build ./src/dicts/ ./lib
	mkdir build lib

	cd include

echo "#############################################"
echo "#         GENERATE ROOT DICTIONARIES        #"
echo "#############################################"
echo " "

	echo "rootcling: ./include/TKEvent.h         -> ./include/TKEventdict.pcm         + ./include/TKEventdict.cpp"
      	      rootcling -f ../lib/TKEventdict.cpp  TKEvent.h+
	echo "rootcling: ./include/TKOMhit.h         -> ./include/TKOMhitdict.pcm         + ./include/TKOMhitdict.cpp"
      	      rootcling -f ../lib/TKOMhitdict.cpp  TKOMhit.h+
	echo "rootcling: ./include/TKtrhit.h         -> ./include/TKtrhitdict.pcm         + ./include/TKtrhitdict.cpp"
      	      rootcling -f ../lib/TKtrhitdict.cpp  TKtrhit.h+
      	echo "rootcling: ./include/TKtrack.h         -> ./include/TKtrackdict.pcm         + ./include/TKtrackdict.cpp"
      	      rootcling -f ../lib/TKtrackdict.cpp  TKtrack.h+
      	echo "rootcling: ./include/TKcluster.h       -> ./include/TKclusterdict.pcm       + ./include/TKclusterdict.cpp"
      	      rootcling -f ../lib/TKclusterdict.cpp  TKcluster.h+
      	echo "rootcling: ./include/TKtrajectory.h     -> ./include/TKtrajectorydict.pcm       + ./include/TKtrajectorydict.cpp"
      	      rootcling -f ../lib/TKtrajectorydict.cpp  TKtrajectory.h+
      	echo "rootcling: ./include/TKpoint.h     -> ./include/TKpointdict.pcm       + ./include/TKpointdict.cpp"
      	      rootcling -f ../lib/TKpointdict.cpp  TKpoint.h+
      	      
      	echo "rootcling: ./include/TKtimer.h     -> ./include/TKtimerdict.pcm       + ./include/TKtimerdict.cpp"
      	      rootcling -f ../lib/TKtimerdict.cpp  TKtimer.h+
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
	cmake .. || { echo $'\n****Cmake failed, installation of TKEvent aborted!****' ; exit 1; } 

echo " "
echo "#############################################"
echo "#         MODULE COMPILATION (MAKE)         #"
echo "#############################################"
echo " "

	make    || { echo $'\n****Make failed, installation of TKEvent aborted!****' ; exit 1; }

	cp libTKEvent.so ../lib/libTKEvent.so
	rm libTKEvent.so
	
