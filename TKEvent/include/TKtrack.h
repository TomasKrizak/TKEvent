#ifndef _TKTRACK_H_
#define _TKTRACK_H_

// Standard headers
#include <iostream>

// ROOT headers
#include "TObject.h"

#include "TKOMhit.h"
#include "TKtrhit.h"

//a line in a form:
//	y = ax + b
//	z = cx + d

class TKtrack: public TObject
{
	private:
		
		int    side;
		double a;
		double b;
		double c;
		double d;
		double likelihood; // metrics for determining quality of track (only in horizontal plane - does not evaluate Z at the moment)
		std::vector<TKtrhit*> associated_tr_hits; // association_distance can be changed in reconstruction functions (3 sigma by default = 6mm)

	public:
		
		TKtrack();
		TKtrack(int _side, double _a, double _b, double _c, double _d);
		~TKtrack();
		
		void add_associated_tr_hit(TKtrhit* trakcer_hit);
		std::vector<TKtrhit*> get_associated_tr_hits();
		
		void set_side(double _side);
		void set_a   (double _a);
		void set_b   (double _b);
		void set_c   (double _c);
		void set_d   (double _d);
		void set_likelihood(double _likelihood);

		int    get_side();
		double get_a   ();
		double get_b   ();
		double get_c   ();
		double get_d   ();
		double get_likelihood();

		void print();
		
		ClassDef(TKtrack,1);
};

#endif
