// TK headers
#include "TKtrajectory.h"
using namespace std;

ClassImp(TKtrajectory);

TKtrajectory::TKtrajectory()
{
	side = -1;
	composite_trajectory = false;
	track_points = std::vector<TKpoint*>();
	segments = std::vector<TKtrack*>();
}

TKtrajectory::TKtrajectory(TKtrack* segment)
{
	//cout << "creating new simple trajectory" << endl;
	side = segment->get_side();
	composite_trajectory = false;
	
	track_points = std::vector<TKpoint*>();
	segments = std::vector<TKtrack*>();
	
	segments.push_back( segment );
	vector<TKpoint*> tr_hit_points = segment->get_associated_tr_hit_points();
	// TODO: Run 974 | Event 1376: only 1 associated hit
	if(tr_hit_points.size() < 2) cout << "WARNING: only 1 associated tracker hit... Cannot create a trajectory" << endl;
	
	double y_min = 2500.0;
	double y_max = -2500.0;
	
	int index_max;
	int index_min;
	// TODO: Run 974 | Event 1158: a = 0 -> does not work
	for(int i = 0; i < tr_hit_points.size(); i++)
	{
		double y = tr_hit_points[i]->get_y();
		if(y >= y_max)
		{
			y_max = y;
			index_max = i; 
		}
		if(y < y_min)
		{
			y_min = y; 
			index_min = i;
		}
	}
	track_points.push_back(tr_hit_points[index_max]);
	track_points.push_back(tr_hit_points[index_min]);
}

TKtrajectory::TKtrajectory(std::vector<TKtrack*> _segments)
{
	//cout << "creating new composite trajectory" << endl;
	segments = std::vector<TKtrack*>();
	for(int i = 0; i < _segments.size(); i++)
	{
		segments.push_back(_segments[i]);
	}
	side = segments[0]->get_side();
	if(segments.size() > 1)
	{
		composite_trajectory = true;
	}
	else
	{
		composite_trajectory = false;
	}
	
	track_points = std::vector<TKpoint*>();	
	
	for(int i = 0; i < segments.size()-1; i++)
	{
		double a1 = segments[i]->get_a();
		double a2 = segments[i+1]->get_a();
		double b1 = segments[i]->get_b();
		double b2 = segments[i+1]->get_b();
		
		double x = (b2-b1)/(a1-a2);
		double y = a1*x + b1;			
		double z1 = segments[i]->get_c()*x + segments[i]->get_d();
		double z2 = segments[i+1]->get_c()*x + segments[i+1]->get_d();
		double z = (z1+z2)/2.0;
		track_points.push_back(new TKpoint(x,y,z));
	}

	vector<TKpoint*> tr_hit_points = segments[0]->get_associated_tr_hit_points();
	double y_kink = track_points[0]->get_y();
	double x_kink = track_points[0]->get_x();
	double y_min = 2500.0;
	double y_max = -2500.0;
	
	int side_positive_counter = 0;
	int side_negative_counter = 0;
	int index_max;
	int index_min;
	for(int i = 0; i < tr_hit_points.size(); i++)
	{
		double y = tr_hit_points[i]->get_y();
		if( y > y_kink ) side_positive_counter++;
		if( y < y_kink ) side_negative_counter++;
		
		
		if(y > y_max)
		{
			y_max = y;
			index_max = i; 
		}
		if(y < y_min)
		{
			y_min = y; 
			index_min = i;
		}
	}
	if( side_positive_counter > side_negative_counter )
	{
		track_points.insert(track_points.begin(), tr_hit_points[index_max]);
	} 
	else
	{
		track_points.insert(track_points.begin(), tr_hit_points[index_min]);
	}
	
	tr_hit_points = segments.back()->get_associated_tr_hit_points();
	y_kink = track_points.back()->get_y();
	x_kink = track_points.back()->get_x();
	y_min = 2500.0;
	y_max = -2500.0;
	
	side_positive_counter = 0;
	side_negative_counter = 0;
	index_max;
	index_min;
	for(int i = 0; i < tr_hit_points.size(); i++)
	{
		double y = tr_hit_points[i]->get_y();
		if( y > y_kink ) side_positive_counter++;
		if( y < y_kink ) side_negative_counter++;
		
		
		if(y > y_max)
		{
			y_max = y;
			index_max = i; 
		}
		if(y < y_min)
		{
			y_min = y; 
			index_min = i;
		}
	}
	if( side_positive_counter > side_negative_counter )
	{
		track_points.push_back(tr_hit_points[index_max]);
	} 
	else
	{
		track_points.push_back(tr_hit_points[index_min]);
	}
}

TKtrajectory::~TKtrajectory()
{
	segments.clear();
	for(int i = 0; i < track_points.size(); i++)
	{
		delete track_points[i];
	}
	track_points.clear();

}

void TKtrajectory::set_side(double _side)
{
	side = _side;
}

int TKtrajectory::get_side()
{
	return side;
} 

std::vector<TKtrack*> TKtrajectory::get_segments()
{
	return segments;
}

void TKtrajectory::add_track_point(TKpoint* track_point)
{
	track_points.push_back(track_point);
}

std::vector<TKpoint*> TKtrajectory::get_track_points()
{
	return track_points;
}

void TKtrajectory::print()
{
	cout <<"Trajectory: " << endl;
	cout << "	side: " << side << endl;
	cout << "	" << segments.size() << " segments: " << endl;
	for(int i = 0; i < segments.size(); i++)
	{
		cout << "	";
		segments[i]->print();
	}
	cout << "	" << track_points.size() << " track points: " << endl;
	for(int i = 0; i < track_points.size(); i++)
	{
		std::cout << "	" << i+1 << ". ";
		track_points[i]->print(); 
	}
}

