// TK headers
#include "TKtimer.h"

using namespace std;

ClassImp(TKtimer);

TKtimer::TKtimer()
{
	cout << "timer created" << endl;

	start = clock();
	global_time = 0.0;
	
	time_reconstruct_ML = 0.0;
	time_clustering_hough = 0.0;
	time_clustering_legendre = 0.0;
	time_build_trajectories = 0.0;
	time_extrapolate = 0.0;

	count_reconstruct_ML = 0;
	count_clustering_hough = 0;
	count_clustering_legendre = 0;
	count_build_trajectories = 0;
	count_extrapolate = 0;
}

TKtimer::~TKtimer()
{	
	this->print();
}


void TKtimer::print()
{
	std::cout << "timer info: " << std::endl;
	std::cout << "	global time: " << (clock() - start) / (double) CLOCKS_PER_SEC << " seconds." <<endl;
	
	std::cout << "	reconstruct ML: " << time_reconstruct_ML << "	count: " << count_reconstruct_ML << std::endl;
	std::cout << "	clustering hough: " << time_clustering_hough << "	count: " << count_clustering_hough << std::endl;
	std::cout << "	clustering legendre: " << time_clustering_legendre << "	count: " << count_clustering_legendre << std::endl;
	std::cout << "	build trajectories: " << time_build_trajectories << "	count: " << count_build_trajectories << std::endl;
	std::cout << "	extrapolating: " << time_extrapolate << "	count: " << count_extrapolate << std::endl;
	double residual_time = (clock() - start) / (double) CLOCKS_PER_SEC -
	 (time_reconstruct_ML + time_clustering_hough + time_clustering_legendre + 
	 time_build_trajectories + time_create_trajectories + time_extrapolate);
	std::cout << "	others: " << residual_time << std::endl; 
}



