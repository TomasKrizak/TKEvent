#ifndef _TKCLUSTER_H_
#define _TKCLUSTER_H_

// Standard headers
#include <iostream>
#include <math.h>

// ROOT headers
#include "TObject.h"

#include "TKtrhit.h"
#include "TKtrack.h"


class TKcluster: public TObject
{
	private:
		
		int    side;
		double phi_min;
		double phi_max;
	
		std::vector<TKtrhit*> cluster_tr_hits;
		TKtrack *track;

	public:
		
		TKcluster();
		TKcluster(std::vector<TKtrhit*> tr_hits, double _phi_min, double _phi_max);
		~TKcluster();
		
		void add_tr_hit(TKtrhit* tracker_hit);
		std::vector<TKtrhit*> get_tr_hits();
		
		void set_track(TKtrack *_track);
		TKtrack* get_track();
		
		void set_side(double _side);
		void set_phi_min(double _phi_min);
		void set_phi_max(double _phi_max);

		int    get_side();
		double get_phi_min();
		double get_phi_max();

		void print();
		
		ClassDef(TKcluster,1);
};

#endif
