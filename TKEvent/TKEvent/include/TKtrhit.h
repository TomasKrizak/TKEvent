// Standard headers
#include <iostream>

// ROOT headers
#include "TObject.h"

class TKtrhit: public TObject
{
	private:
// MIRO: Pridať is_corner, is_side, is_inside, is_broken - určí driftmodel
// MIRO: Neskôr pridať aj double deltaR a deltaZ neistoty keď s tým bude vedieť drift model pracovať		
		int    cell_num;
		int    SRL[3];	// 0 = Side, 1 = Row, 2 = Layer
		double xy [2];	// (x,y) coordinate of the tracker cell xy[0] = x, xy[1] = y

		int64_t tsp[7]; // 0-4 = timestamp 0-4, 5 = cathode bottom, 6 = cathode top

		double h;
		double r;

		void set_SRL_xy();

	public:
		
		TKtrhit();
		TKtrhit(int _cell_num);
		TKtrhit(int _SRL[3]);
		TKtrhit(int _cell_num, int64_t _tsp[7]);
		TKtrhit(int _SRL[3],   int64_t _tsp[7]);
		~TKtrhit();

		void set_cell_num(int _cell_num);
		void set_tsp	 (int64_t _tsp[7]); // sets timestamps

		void set_h();
		void set_r(double _r);

		int     get_cell_num  ();
		int     get_SRL       (char _SRL_n);
		double  get_xy        (char _xy_n);
		int64_t get_tsp       (char _tsp_n); // returns timestamp
		double  get_r         ();
		double  get_h         ();

		void print();
		
		ClassDef(TKtrhit,1);
};
