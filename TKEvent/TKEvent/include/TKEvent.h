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

// TK headers
#include "TKOMhit.h"
#include "TKtrhit.h"
#include "TKtrack.h"

// MIRO: PRIDAŤ TIE DIREKTIVY PREPROCESORU ČI AKO SA TO VOLÁ IFDEF A TAK - ABY DVAKRÁT NEINCLUDOVAL TO ISTÉ

class TKEvent: public TObject
{
	private:
		
		int run_number;
		int event_number;
		
		std::vector<TKOMhit*> OM_hits;
		std::vector<TKtrhit*> tr_hits;
		std::vector<TKtrack*> tracks;
	
// MIRO: Pridať triedu TKBisource?	
		
	public:
		TKEvent();
		TKEvent(int _run_number ,int _event_number);
		~TKEvent();

		std::vector<TKOMhit*> get_OM_hits();
		std::vector<TKtrhit*> get_tr_hits();
		std::vector<TKtrack*> get_tracks();
		TKOMhit* 	      get_OM_hit(int _i);
		TKtrhit* 	      get_tr_hit(int _i);
		TKtrack* 	      get_track(int _i);
	
		int get_run_number();
		int get_event_number();

		int  get_no_tracks();
		void print();
		void print_tracks();

		void add_OM_hit(int _OM_num,  bool _is_HT, int64_t _OM_TDC, int16_t _OM_pcell);
		void add_OM_hit(int _SWCR[4], bool _is_HT, int64_t _OM_TDC, int16_t _OM_pcell);
		
		void add_tracker_hit(int _cell_num, int64_t _tsp[7]);
		void add_tracker_hit(int _SRL[3],   int64_t _tsp[7]);
		
		void reconstruct_track(bool save_sinograms);
		void make_top_projection();
		void build_event();	

		// drift model (you can switch between "Manchester" and "Betsy")
		void set_r(std::string _model_n);				

		ClassDef(TKEvent,1);	
};

#endif
