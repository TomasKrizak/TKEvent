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

// TK headers
#include "TKOMhit.h"
#include "TKtrhit.h"

// MIRO: PRIDAŤ TIE DIREKTIVY PREPROCESORU ČI AKO SA TO VOLÁ IFDEF A TAK - ABY DVAKRÁT NEINCLUDOVAL TO ISTÉ

class TKEvent: public TObject
{
	private:
		
		int run_number;
		int event_number;
		
		std::vector<TKOMhit*> OM_hits;
		std::vector<TKtrhit*> tr_hits;

// MIRO: Pridať triedu TKtrack	
// MIRO: Pridať sinogramy?
// MIRO: Pridať triedu TKBisource?	
		std::vector<double> rec_tracks_side;
		std::vector<double> rec_tracks_a;
		std::vector<double> rec_tracks_b;
		std::vector<double> rec_tracks_c;
		std::vector<double> rec_tracks_d;
		
	public:
		TKEvent();
		TKEvent(int _run_number ,int _event_number);
		~TKEvent();

		std::vector<TKOMhit*> get_OM_hits();
		std::vector<TKtrhit*> get_tr_hits();
		TKOMhit* 	      get_OM_hit(int _i);
		TKtrhit* 	      get_tr_hit(int _i);
	
		int get_run_number();
		int get_event_number();
		int get_track_side(int index);
		double get_track_a(int index);
		double get_track_b(int index);
		double get_track_c(int index);
		double get_track_d(int index);

		int  get_no_tracks();
		void print();
		void print_tracks();

		void add_OM_hit(int _OM_num,  bool _is_HT, int64_t _OM_TDC, int16_t _OM_pcell);
		void add_OM_hit(int _SWCR[4], bool _is_HT, int64_t _OM_TDC, int16_t _OM_pcell);
		
		void add_tracker_hit(int _cell_num, int64_t _tsp[7]);
		void add_tracker_hit(int _SRL[3],   int64_t _tsp[7]);
		
		void reconstruct_track();
		void make_top_projection();
		void build_event();	

		// drift model (you can switch between "Manchester" and "Betsy")
		void set_r(std::string _model_n);				

		ClassDef(TKEvent,1);	
};
