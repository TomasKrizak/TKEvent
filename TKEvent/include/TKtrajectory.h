#ifndef _TKTRAJECTORY_H_
#define _TKTRAJECTORY_H_

// Standard headers
#include <iostream>
#include <vector>
#include <cmath>

// ROOT headers
#include "TObject.h"

#include "TKOMhit.h"
#include "TKtrhit.h"
#include "TKtrack.h"
#include "TKpoint.h"

class TKtrajectory: public TObject
{
	private:
		
		int side;
		bool composite_trajectory;
		std::vector<TKpoint*> track_points;
		std::vector<TKtrack*> segments;
		
	public:
		
		TKtrajectory();
		TKtrajectory(TKtrack* _segment);
		TKtrajectory(std::vector<TKtrack*> _segments);
		~TKtrajectory();
		
		void set_side(double _side);		
		void add_track_point(TKpoint* track_point);

		int get_side();
		std::vector<TKtrack*> get_segments();
		std::vector<TKpoint*> get_track_points();
		
		void extrapolate();

		void print();
		
		ClassDef(TKtrajectory,1);
};

#endif
