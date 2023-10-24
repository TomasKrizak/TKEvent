// TK headers
#include "TKcluster.h"

ClassImp(TKcluster);


TKcluster::TKcluster()
{
	phi_min = 0.0;
	phi_max = M_PI;
	std::vector<TKtrhit*> cluster_tr_hits;
}

TKcluster::TKcluster(std::vector<TKtrhit*> tr_hits, double _phi_min, double _phi_max)
{
	side = tr_hits.at(0)->get_SRL('s');
	phi_min = _phi_min;
	phi_max = _phi_max;
	cluster_tr_hits = std::vector<TKtrhit*>();
	for(int i = 0; i < tr_hits.size(); i++)
	{
		cluster_tr_hits.push_back( tr_hits.at(i) );
	}
}

TKcluster::~TKcluster()
{
	cluster_tr_hits.clear();
}

void TKcluster::add_tr_hit(TKtrhit* tracker_hit)
{
	cluster_tr_hits.push_back( tracker_hit );
}

std::vector<TKtrhit*> TKcluster::get_tr_hits()
{
	return cluster_tr_hits;
}

void TKcluster::set_track(TKtrack *_track)
{
	track = _track;
}

TKtrack* TKcluster::get_track()
{
	return track;
}

void TKcluster::set_side(double _side)
{
	side = _side;
}

void TKcluster::set_phi_min(double _phi_min)
{
	phi_min = _phi_min;
}

void TKcluster::set_phi_max(double _phi_max)
{
	phi_max = _phi_max;
}

int    TKcluster::get_side()
{
	return side;
}

double TKcluster::get_phi_min()
{
	return phi_min;
}

double TKcluster::get_phi_max()
{
	return phi_max;
}

void TKcluster::print()
{
	std::cout << "	side "       << side
		  << ", size "       << cluster_tr_hits.size()
	     	  << ", phi_min: " << phi_min 
	     	  << ", phi_max: "   << phi_max << std::endl;
}
