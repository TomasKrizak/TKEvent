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
	side = segment->get_side();
	composite_trajectory = false;
	
	track_points = std::vector<TKpoint*>();
	segments = std::vector<TKtrack*>();
	
	segments.push_back( segment );
	vector<TKpoint*> tr_hit_points = segment->get_associated_tr_hit_points();
	// TODO: Run 974 | Event 1376: only 1 associated hit
	if(tr_hit_points.size() < 2)
	{
		cout << "WARNING: only 1 associated tracker hit... Cannot create a trajectory" << endl;
	}
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

void TKtrajectory::extrapolate()
{
	const double max_allowed_extrapolation = 150.0;

	const double mainwall_x_position = 435.0;
	const double Xwall_y_position = 2505.5;
	const double Gveto_z_position = 1550.0;

	if(track_points.size() < 2) return;
	
	double x0(0.0), y0(0.0), z0(0.0);
	double x1 = track_points[0]->get_x();
	double y1 = track_points[0]->get_y();
	double z1 = track_points[0]->get_z();
	double x2 = track_points[1]->get_x();
	double y2 = track_points[1]->get_y();
	double z2 = track_points[1]->get_z();

	bool found = false;
	if(x1 != x2)
	{
		x1 < x2 ? x0 = mainwall_x_position * (side-1) : x0 = side*mainwall_x_position; 
		double Cx = (x1-x0)/(x2-x1);
		y0 = y1 - Cx*(y2-y1);
		z0 = z1 - Cx*(z2-z1);
		if( Xwall_y_position > abs(y0) && Gveto_z_position > abs(z0) )
		{
			found = true;
		}
	}
	if(!found && y1 != y2)
	{
		y1 < y2 ? y0 = -Xwall_y_position : y0 = Xwall_y_position; 
		double Cy = (y1-y0)/(y2-y1);
		x0 = x1 - Cy*(x2-x1);
		z0 = z1 - Cy*(z2-z1);
		if( mainwall_x_position/2.0 > abs(x0 - (double(side)-0.5)*mainwall_x_position) && Gveto_z_position > abs(z0) )
		{
			found = true;
		}
	}
	if(!found && z1 != z2)
	{
		z1 < z2 ? z0 = -Gveto_z_position : z0 = Gveto_z_position; 
		double Cz = (z1-z0)/(z2-z1);
		x0 = x1 - Cz*(x2-x1);
		y0 = y1 - Cz*(y2-y1);
		if( mainwall_x_position/2.0 > abs(x0 - (double(side)-0.5)*mainwall_x_position) && Xwall_y_position > abs(y0) )
		{
			found = true;
		}
	}
	
	TKpoint* start_point = new TKpoint(x0, y0, z0);
	if( distance_2D( *start_point, *track_points[0]) < max_allowed_extrapolation )
	{
		track_points.insert(track_points.begin(), start_point);
	}
	else
	{
		delete start_point;
	}
	x0 = 0.0;
	y0 = 0.0; 
	z0 = 0.0;
	x1 = track_points.end()[-1]->get_x();
	y1 = track_points.end()[-1]->get_y();
	z1 = track_points.end()[-1]->get_z();
	x2 = track_points.end()[-2]->get_x();
	y2 = track_points.end()[-2]->get_y();
	z2 = track_points.end()[-2]->get_z();

	found = false;
	if(x1 != x2)
	{
		x1 < x2 ? x0 = mainwall_x_position * (side-1) : x0 = side*mainwall_x_position; 
		double Cx = (x1-x0)/(x2-x1);
		y0 = y1 - Cx*(y2-y1);
		z0 = z1 - Cx*(z2-z1);
		if( Xwall_y_position > abs(y0) && Gveto_z_position > abs(z0) )
		{
			found = true;
		}
	}
	if(!found && y1 != y2)
	{
		y1 < y2 ? y0 = -Xwall_y_position : y0 = Xwall_y_position; 
		double Cy = (y1-y0)/(y2-y1);
		x0 = x1 - Cy*(x2-x1);
		z0 = z1 - Cy*(z2-z1);
		if( mainwall_x_position/2.0 > abs(x0 - (double(side)-0.5)*mainwall_x_position) && Gveto_z_position > abs(z0) )
		{
			found = true;
		}
	}
	if(!found && z1 != z2)
	{
		z1 < z2 ? z0 = -Gveto_z_position : z0 = Gveto_z_position; 
		double Cz = (z1-z0)/(z2-z1);
		x0 = x1 - Cz*(x2-x1);
		y0 = y1 - Cz*(y2-y1);
		if( mainwall_x_position/2.0 > abs(x0 - (double(side)-0.5)*mainwall_x_position) && Xwall_y_position > abs(y0) )
		{
			found = true;
		}
	}
	
	TKpoint* end_point = new TKpoint(x0, y0, z0);
	if( distance_2D( *end_point, *track_points.end()[-1]) < max_allowed_extrapolation )
	{
		track_points.push_back(end_point);
	}
	else
	{
		delete end_point;
	}
	
	
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

