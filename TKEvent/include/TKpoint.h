#ifndef _TKPOINT_H_
#define _TKPOINT_H_

// Standard headers
#include <iostream>
#include <cmath>

// ROOT headers
#include "TObject.h"

class TKpoint: public TObject
{
	private:
		
		double x;
		double y;
		double z;
		int type; 

	public:
		
		TKpoint();
		TKpoint(double _x, double _y, double _z);
		~TKpoint();
		
		void set_x(double _x);
		void set_y(double _y);
		void set_z(double _z);
		
		double get_x();
		double get_y();
		double get_z();

		void print();
		
		ClassDef(TKpoint,1);
		
		
};

double distance_2D(TKpoint &point1, TKpoint &point2);
double distance_3D(TKpoint &point1, TKpoint &point2);

#endif
