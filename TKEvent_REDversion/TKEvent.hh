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

#include <datatools/clhep_units.h>

#include <vector>
#include <iostream>

// detector dimensions in mm
// origin in the center of detector

static double foil_spacex = 58.0; // probably wrong

const double tc_radius = 22.0;
const double tc_sizez = 3030.0;

const double mw_sizex = 194.0;
const double mw_sizey = 256.0;
const double mw_sizez = 256.0;

const double gv_sizex = 308.0;
const double gv_sizey = 310.0;
const double gv_sizez = 150.0;

const double xw_sizex = 200.0;
const double xw_sizey = 150.0;
const double xw_sizez = 208.5;

const double Bi_source_x = 1.4;
const double Bi_source_y = 14.2;
const double Bi_source_z = 21.1;
const double Bi_source_dist_y = 835.0; // possibly wrong or not accurate (originally 850)
const double Bi_source_dist_z = 475.0; // possibly wrong or not accurate (originally 425)

class TKEvent
{
	private:
		
		int run_number;
		int event_number;
		
		std::vector<int>     hit_om_num;
		std::vector<bool>    hit_is_HT;
		std::vector<int64_t> hit_om_tdc;
		std::vector<int16_t> hit_om_peak_cell;
		
		std::vector<int>    hit_cell_num;
		std::vector<double> hit_height;
		std::vector<double> hit_radius;
		
		std::vector<double> rec_tracks_side;
		std::vector<double> rec_tracks_a;
		std::vector<double> rec_tracks_b;
		std::vector<double> rec_tracks_c;
		std::vector<double> rec_tracks_d;
		
	public:
		TKEvent(int _run_number ,int _event_number);
		~TKEvent();
	
		int get_run_number();
		int get_event_number();
		int get_track_side(int index);
		double get_track_a(int index);
		double get_track_b(int index);
		double get_track_c(int index);
		double get_track_d(int index);
		int get_no_tracks();
		
		void add_event(snemo::datamodel::raw_event_data red);
		void add_calo_hit(const snemo::datamodel::calo_digitized_hit & red_calo_hit);
		void add_calo_hit(int _hit_om_num, bool _hit_is_HT, int64_t calo_tdc, int16_t peak_cell);
		void add_calo_hit(int cell_side, int cell_wall, int cell_column, int cell_row, bool is_HT, int64_t calo_tdc, int16_t peak_cell);
		
		void add_tracker_hit(const snemo::datamodel::tracker_digitized_hit & red_tracker_hit);
		void add_tracker_hit(int _hit_cell_num, double _hit_height, double _hit_radius);
		void add_tracker_hit(int hit_cell_side, int hit_cell_row, int hit_cell_layer, double hit_height, double hit_radius);
		
		void reconstruct_track();
		void make_top_projection();
		void build_event_3D();	
		
		double static get_tracker_hit_height(const snemo::datamodel::tracker_digitized_hit::gg_times & gg_timestamps);
		
		// drift model (you can switch between "Manu" and "Betsy")
		double static drift_time_to_radius(double time_, std::string type);
				
						
		// tracker: [side, row, layer] to Cell number
		int static SRL_to_Cellnum(int cell_side, int cell_row, int cell_layer);

		// tracker: Cell number to [side, row, layer]
		void static Cellnum_to_SRL(int cell_num, int SRL[]);

		// tracker: [side, row, layer] to [x,y] coordinates (centers of tracker cells)
		void static SRL_to_xy(int cell_side, int cell_row, int cell_layer, double position_xy[]);
				
		// tracker: Cell number to [x,y] coordinates (centers of tracker cells)
		void static Cellnum_to_xy(int cell_num, double position_xy[]);	


		// calo: [side, wall, column, row] to OM number
		int static SWCR_to_OMnum(int om_side, int om_wall, int om_column, int om_row);

		// calo: OM number to [side, wall, column, row]
		void static OMnum_to_SWCR(int om_num, int OM[]);

		// calo: [side, wall, column, row] to [x,y,z] coordinates (centers of calos)
		void static SWCR_to_xyz(int om_side, int om_wall, int om_column, int om_row, double position_xyz[]);
		// warning: dimensions are approximately taken from Falaise and might not be exactly right

		// calo: OM number to [x,y,z] coordinates (centers of calos)
		void static OMnum_to_xyz(int om_num, double position_xyz[]);
		// warning: dimensions are approximately taken from Falaise and might not be exactly right

};
