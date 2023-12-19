#ifndef _TKTRHIT_H_
#define _TKTRHIT_H_

// Standard headers
#include <iostream>

// ROOT headers
#include "TObject.h"

#include "TKOMhit.h"
class TKtrack;

class TKtrhit: public TObject
{
	private:
// MIRO: Pridať is_corner, is_side, is_inside, is_broken - určí drift model
		
		int    cell_num;
		int    SRL[3];	// 0 = Side, 1 = Row, 2 = Layer
		double xy [2];	// (x,y) coordinate of the tracker cell xy[0] = x, xy[1] = y

		int64_t tsp[7]; // 0-4 = timestamp 0-4, 5 = cathode bottom, 6 = cathode top
		TKOMhit *associated_OMhit;
		TKtrack *associated_track;

		// vertical position of the hit
		double h;
		double sigma_Z; // default value: 17.0 mm (Gaussian model)

		// r = drift radius 		
		double r;
		double sigma_R; // default value: 2.0 mm (Gaussian model)
		
		void set_SRL_xy();

	public:
		
		TKtrhit();
		TKtrhit(int _cell_num);
		TKtrhit(int _SRL[3]);
		TKtrhit(int _cell_num, int64_t _tsp[7]);
		TKtrhit(int _SRL[3],   int64_t _tsp[7]);
		~TKtrhit();

		void set_cell_num(int _cell_num);
		void set_tsp	  (int64_t _tsp[7]); // sets timestamps

		void set_h(double _h);
		void set_h();	// space for better model implementation
		void set_sigma_Z(double _sigma_Z);
		void set_sigma_Z(); // space for better model implementation
		
		void set_r(double _r);
		void set_sigma_R(double _sigma_R);
		void set_sigma_R(); // space for better model implementation
		
		void set_associated_OMhit(TKOMhit *_associated_OMhit);
		void set_associated_track(TKtrack *_associated_track);

		int     get_cell_num  ();
		int     get_SRL       (char _SRL_n);
		double  get_xy        (char _xy_n);
		int64_t get_tsp       (char _tsp_n); // returns timestamp
		double  get_r         ();
		double  get_sigma_R   ();
		double  get_h         ();
		double  get_sigma_Z   ();
		TKOMhit* get_associated_OMhit();
		TKtrack* get_associated_track();
		
		void print();
		
		ClassDef(TKtrhit,1);
};

#endif
