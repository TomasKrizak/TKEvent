#ifndef _TKTRACK_H_
#define _TKTRACK_H_

// Standard headers
#include <iostream>

// ROOT headers
#include "TObject.h"

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

	public:
		
		TKtrack();
		TKtrack(int _side, double _a, double _b, double _c, double _d);
		~TKtrack();
		
		void set_side(double _side);
		void set_a   (double _a);
		void set_b   (double _b);
		void set_c   (double _c);
		void set_d   (double _d);

		int    get_side();
		double get_a   ();
		double get_b   ();
		double get_c   ();
		double get_d   ();

		void print();
		
		ClassDef(TKtrack,1);
};

#endif
