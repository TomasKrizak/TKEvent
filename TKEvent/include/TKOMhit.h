#ifndef _TKOMHIT_H_
#define _TKOMHIT_H_

// Standard headers
#include <iostream>

// ROOT headers
#include "TObject.h"

class TKOMhit: public TObject
{
	private:
// MIRO: Pozrieť či sa dá nejako už pridať metóda  TKOMhit::calibrate() od filipa,  takže bude treba pridať aj double E; - energia, prípadne aj deltaE
// MIRO: Pridať niečo ako OM_type - aby sám vedel či je main wall, XW alebo GV
// MIRO: Pridať operator metódy aby sa dali porovnávať, napríklad operator> ktorý testuje či sú vedľa seba a operator>> či sú o dva ďalej - alebo zvoliť nejaký lepší
		int     OM_num;
		int 	SWCR[4];	// 0 = Side, 1 = Wall, 2 = Column, 3 = Row
		double  xyz[3];	// (x,y,z) coordinate of the center of OM xyz[0] = x, xyz[1] = y, xyz[2] = z
		
		bool    HT;		// high treshold flag
		int32_t charge;
		int16_t amplitude;
		int16_t baseline;
		int64_t OM_TDC;
		int16_t OM_pcell;
		
		void 	set_SWCR_xyz();

	public:
		
		TKOMhit();
		TKOMhit(int _OM_num);
		TKOMhit(int _SWCR[4]);
		TKOMhit(int _OM_num,  bool _HT, int64_t _OM_TDC, int16_t _OM_pcell);
		TKOMhit(int _SWCR[4], bool _HT, int64_t _OM_TDC, int16_t _OM_pcell);
		~TKOMhit();
		
		void set_OM_num   (int     _OM_num);
		
		void set_HT       (bool    _HT);
		void set_charge   (int32_t _charge);
		void set_amplitude(int16_t _amplitude);
		void set_baseline (int16_t _baseline);
		void set_OM_TDC   (int64_t _OM_TDC);
		void set_OM_pcell (int16_t _OM_pcell);

		int     get_OM_num   ();
		int     get_SWCR     (char _SWCR_n);
		double  get_xyz      (char _xyz_n);
		
		bool    is_HT        ();
		int32_t get_charge   ();
		int16_t get_amplitude();
		int16_t get_baseline ();
		int64_t get_OM_TDC   ();
		int64_t get_OM_pcell ();

		void print();
		
		ClassDef(TKOMhit,1);
};

#endif
