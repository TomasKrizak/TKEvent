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
#include "TPolyLine.h"
#include "TBox.h"
#include "TLatex.h"
#include "TLine.h"
#include "TPoint.h"
#include "TObject.h"
#include <TStyle.h>
#include <TF1.h>

// TK headers
#include "TKOMhit.h"
#include "TKtrack.h"
#include "TKtrhit.h"
#include "TKcluster.h"
#include "TKpoint.h"
#include "TKtrajectory.h"


// note (9.4.): currently the main functions to use are the following:

// 	reconstruct(bool save_sinograms); full reconstruction including trajectory 
//					   builder and verteces extrapolation
//	reconstruct_ML(bool save_sinograms); for quick track finding
//	make_top_projection(int hits_option, int tracking_option);
//	build_event(int tracking_option);

// In case you are looking for the implementation - all functions concerning 
// the tracking are implemented in a file "tracking_tools.cpp"

class TKEvent: public TObject
{
	private:
		
		int run_number;
		int event_number;
		
		std::vector<TKOMhit*>   OM_hits;
		std::vector<TKtrhit*>   tr_hits;
		std::vector<TKtrack*>   tracks;
		std::vector<TKcluster*> clusters;
		std::vector<TKtrajectory*> trajectories;
		
	public:
	
	// basic functionality section:
	
		TKEvent();
		TKEvent(int _run_number ,int _event_number);
		~TKEvent();

		std::vector<TKOMhit*>&      get_OM_hits();
		std::vector<TKtrhit*>&      get_tr_hits();
		std::vector<TKtrack*>       get_tracks();
		std::vector<TKcluster*>&    get_clusters();
		std::vector<TKtrajectory*>& get_trajectories();
		
		TKOMhit*      get_OM_hit(int _i);
		TKtrhit*      get_tr_hit(int _i);
		TKcluster*    get_cluster(int _i);
		TKtrajectory* get_trajectory(int _i);
	
		int get_run_number();
		int get_event_number();
		int get_no_tracks();
		int get_no_trajectories();
		
		void print();
		void print_tracks();
		void print_trajectories();

		void add_OM_hit(int _OM_num,  bool _is_HT, int64_t _OM_TDC, int16_t _OM_pcell);
		void add_OM_hit(int _SWCR[4], bool _is_HT, int64_t _OM_TDC, int16_t _OM_pcell);
		
		void add_tracker_hit(int _cell_num, int64_t _tsp[7]);
		void add_tracker_hit(int _SRL[3],   int64_t _tsp[7]);
	
	// drift model and plasma propagation section
	
		// associates tracker hits to OM hits and calculates hit radii
		// drift model: "Manchester" or "Betsy"
		// association_mode: "time" or "distance"
		// 	"distance": minimazes distance between OM and tracker hit
		// 	"time": minimazes time difference between OM and tracker hit
		void set_r(std::string drift_model = "Manchester", std::string association_mode = "distance");
		
		// calls sigma_R, sigma_Z, h set functions for all tracker hits
		void set_sigma_R(); // currently unnecessary
		void set_h();
		void set_sigma_Z(); // currently unnecessary
	
	// tracker hit collection filtering
	
		std::vector<TKtrhit*> filter_side(std::vector<TKtrhit*>& _hits, int side);
		std::vector<TKtrhit*> filter_usable(std::vector<TKtrhit*>& _hits);
		std::vector<TKtrhit*> filter_unassociated(std::vector<TKtrhit*>& _hits);
		std::vector<TKtrhit*> filter_unclustered(std::vector<TKtrhit*>& _hits);
		std::vector<TKtrhit*> filter_distant(std::vector<TKtrhit*>& _hits);
		std::vector<TKtrhit*> filter_close_hits(std::vector<TKtrhit*>& _hits, double phi, double r, double distance_limit);
		
	// clustering 
	
		// Hough transform based clutering - finds a largest subgroup of given hits
		// such that is geometrically possible to have a single common line
		TKcluster* find_cluster(std::vector<TKtrhit*>& tr_hits);
		// Legendre based clustering
		TKcluster* find_cluster_legendre(std::vector<TKtrhit*>& hits, bool save_sinograms = false);
	
	// full reconstruction functions
	
		// full reconstruction algorithm:
		//	1. different clusterings to safe failed events 
		//	2. maximum likelihood to obtain line tracks 
		//	3. ambiguity checker and solver
		//	4. trajectory builder from found segments
		//	5. trajectory extrapolator
		 
		void reconstruct(bool save_sinograms = false); // full reconstruction
		void reconstruct_simple(bool save_sinograms = false); // simpler quick algo for one track per side
		
	// line track reconstruction section	
	
		// basic reconstruction - no uncertainties, one candidate
		void reconstruct_track(bool save_sinogram = false);
		void reconstruct_track_from_hits(std::vector<TKtrhit*>& hits, bool save_sinograms = false);
		
		// with uncertainties, one candidate - recommended function
		void reconstruct_single(bool save_sinograms = false);
		void reconstruct_single_from_hits(std::vector<TKtrhit*>& hits, bool save_sinograms = false);
				
		// with uncertainties, multiple candidates
		void reconstruct_multi(bool save_sinograms = false);

		// reconstruction based on maximum likelihood - currently best algorithm
		// a combination of basic clustering and maximum likelihood method
		// currently finds only solution per detector side
		void reconstruct_ML(bool save_sinograms = false);
		void reconstruct_ML_3D(bool save_sinograms = false);
	
	// trajectory builder
	
		void calculate_tr_hit_points();
		void build_trajectories();
		void extrapolate_trajectories();
		
	// vizualization section
		
		// tracker hits options:
		// 	0 - no unused hits	
		//	    red	= used hits for reconstruction    
		//
		//	1 - red	= used hits for reconstruction
		//	    yellow 	= unused hits for recontstruction
		//
		//	2 - red	= used hits for reconstruction (unassociated)
		//	    yellow 	= unused hits for recontstruction 
		//	    green 	= associated hits to track
		//
		// 	3 - red	= used hits for reconstruction (unassociated + good vertical position)
		//	    yellow 	= unused hits for recontstruction
		//	    magenta 	= failed vertical position reconstruction but good drift radius (unassociated)
		//	    green 	= associated hits to track, good vertical position
		//	    teal	= associated hits to track, failed vertical position
				 
		// tracking options:
		//	0 - only tracks
		//	1 - tracks
		//	    reconstructed tracker hit avalanche origin points
		//	2 - trajectories
		//	    avalanche origin points
		//	3 - tracks
		//	    trajectories
		//	    avalanche origin points
				 
				 
		void make_top_projection(int hits_option = 3, int tracking_option = 3);
		void build_event(int tracking_option = 3);	

	// tools for drawing certain mathematical functions 
		void hough_transform(std::vector<TKtrhit*>& hits, double phi_min, double phi_max, double R_min, double R_max, int ID);
		void draw_likelihood();
		void draw_likelihood_centred();
		void draw_sinusoids();

		ClassDef(TKEvent,1);	
};

#endif
