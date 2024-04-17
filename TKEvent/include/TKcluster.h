#ifndef _TKCLUSTER_H_
#define _TKCLUSTER_H_

// Standard headers
#include <iostream>
#include <math.h>

// ROOT headers
#include "TCanvas.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TF1.h"
#include "TROOT.h"
#include "TGraph.h"
#include "TObject.h"

#include "TKtrhit.h"
#include "TKtrack.h"


class TKcluster: public TObject
{
	private:
		
		int    side;
		double phi_min;
		double phi_max;
		
		// 0 == no ambiguity
		// 1 == mirror image along line x = x0 
		// 2 == mirror image along line y = y0 
		// 3 == mirror image along line y = x + (y0-x0) 
		// 4 == mirror image along line y = -x + (y0-x0) 
		int ambiguity_type;
		
		std::vector<TKtrhit*> cluster_tr_hits;
		TKtrack* track;

	public:
		
		TKcluster();
		TKcluster(std::vector<TKtrhit*>& tr_hits, double _phi_min, double _phi_max);
		~TKcluster();
		
		void add_tr_hit(TKtrhit* tracker_hit);
		std::vector<TKtrhit*> get_tr_hits();
		
		void set_track(TKtrack* _track);
		TKtrack* get_track();
		
		void set_side(double _side);
		void set_phi_min(double _phi_min);
		void set_phi_max(double _phi_max);
		
		int    get_side();
		double get_phi_min();
		double get_phi_max();

		void detect_ambiguity_type();
		int get_ambiguity_type();
		
		void reconstruct_ambiguity();
		void reconstruct_MLM(bool save_sinograms, int run_number, int event_number);
		void reconstruct_MLM_3D(bool save_sinograms, int run_number, int event_number);
		void draw_ML_vertical(int run_number, int event_number);

		void print();
		
		ClassDef(TKcluster,1);
};

#endif
