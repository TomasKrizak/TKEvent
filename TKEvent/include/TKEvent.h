#ifndef _TKEVENT_H_
#define _TKEVENT_H_

// Standard headers
#include <vector>
#include <iostream>

// ROOT headers
#include "TCanvas.h"
#include "TColor.h"
#include "TEllipse.h"
#include "TH2D.h"
#include "TAttLine.h"
#include "TGLViewer.h"
#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TFile.h"
#include "TROOT.h"
#include "TPolyLine3D.h"
#include "TBox.h"
#include "TLatex.h"
#include "TLine.h"
#include "TObject.h"
#include <TStyle.h>
#include <TF1.h>

// TK headers
#include "TKOMhit.h"
#include "TKtrack.h"
#include "TKtrhit.h"
#include "TKcluster.h"


// note (29.10.2023): currently the main functions to use are the following:

// 	reconstruct_single(bool save_sinograms);
//	reconstruct_multi(bool save_sinograms);
//	reconstruct_ML(bool save_sinograms);
//	make_top_projection(int option);
//	build_event();

// rest of the functions are important mainly for tracking development or as 
// inner functions of important reconstruction functions listed above

// In case you are looking for the implementation - all functions concerning 
// maximum likelihood method are implemented in a file "tracking_tools.cpp"

class TKEvent: public TObject
{
	private:
		
		int run_number;
		int event_number;
		
		std::vector<TKOMhit*>   OM_hits;
		std::vector<TKtrhit*>   tr_hits;
		std::vector<TKtrack*>   tracks;
		std::vector<TKcluster*> clusters;
	
	public:
	
	// basic functionality section:
	
		TKEvent();
		TKEvent(int _run_number ,int _event_number);
		~TKEvent();

		std::vector<TKOMhit*>   get_OM_hits();
		std::vector<TKtrhit*>   get_tr_hits();
		std::vector<TKtrack*>   get_tracks();
		std::vector<TKcluster*> get_clusters();
		TKOMhit*   get_OM_hit(int _i);
		TKtrhit*   get_tr_hit(int _i);
		TKtrack*   get_track(int _i);
		TKcluster* get_cluster(int _i);
	
		int get_run_number();
		int get_event_number();
		int get_no_tracks();
		
		void print();
		void print_tracks();

		void add_OM_hit(int _OM_num,  bool _is_HT, int64_t _OM_TDC, int16_t _OM_pcell);
		void add_OM_hit(int _SWCR[4], bool _is_HT, int64_t _OM_TDC, int16_t _OM_pcell);
		
		void add_tracker_hit(int _cell_num, int64_t _tsp[7]);
		void add_tracker_hit(int _SRL[3],   int64_t _tsp[7]);
	
	// possibly future drift model section - currently only radius calculation
	
		// drift model: "Manchester" or "Betsy"
		// association_mode: "time" or "distance"
		// 	"distance": minimazes distance between OM and tracker hit
		// 	"time": minimazes time difference between OM and tracker hit
		void set_r(std::string drift_model, std::string association_mode);	
	
	// clustering functions section
	
		std::vector<TKtrhit*> filter_side(std::vector<TKtrhit*> _hits, int side);
		std::vector<TKtrhit*> filter_usable(std::vector<TKtrhit*> _hits);
		std::vector<TKtrhit*> filter_close_hits(std::vector<TKtrhit*> _hits, double phi, double r, double distance_limit);
		
		// basic clustering - finds a largest subgroup of given hits such that is geometrically possible to have a single common line
		TKcluster* find_cluster(std::vector<TKtrhit*> tr_hits);
		
	// reconstruction section	
	
		// basic reconstruction - no uncertainties, one candidate
		void reconstruct_track(bool save_sinograms);
		void reconstruct_track_from_hits(std::vector<TKtrhit*> hits, bool save_sinograms);
		
		// with uncertainties, one candidate - recommended function
		void reconstruct_single(bool save_sinograms);
				
		// with uncertainties, multiple candidates
		void reconstruct_multi(bool save_sinograms);

		// reconstruction based on maximum likelihood - currently best algorithm
		// a combination of basic clustering and maximum likelihood method
		// currently finds only solution per detector side
		void reconstruct_ML(bool save_sinograms);
		
	// vizualization section
		
		// options:
		// 	0 - no unused hits	    
		//	1 - yellow for unused hits in recontstruction
		//	2 - yellow for unused hits in recontstruction 
		//	    green for associated hits to track
		void make_top_projection(int option);
		void build_event();	

	// tools for drawing certain mathematical functions 
		void hough_transform(std::vector<TKtrhit*> hits, double phi_min, double phi_max, double R_min, double R_max, int ID);
		void draw_likelihood();
		void draw_likelihood_centred();
		void draw_sinusoids();

		ClassDef(TKEvent,1);	
};

#endif
