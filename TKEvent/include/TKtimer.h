#ifndef _TKTIMER_H_
#define _TKTIMER_H_

// Standard headers
#include <iostream>
#include <cmath>

// ROOT headers
#include "TObject.h"

class TKtimer: public TObject
{		
	private:
		clock_t start;	
		
	public:
		double global_time;
		
		double time_reconstruct_ML;
		double time_clustering_hough;
		double time_clustering_legendre;
		double time_build_trajectories;
		double time_create_trajectories;
		double time_extrapolate;
		
		int count_reconstruct_ML;
		int count_clustering_hough;
		int count_clustering_legendre;
		int count_build_trajectories;
		int count_create_trajectories;
		int count_extrapolate;

		TKtimer();
		~TKtimer();

		void print();	
							
		ClassDef(TKtimer,1);
};

#endif
