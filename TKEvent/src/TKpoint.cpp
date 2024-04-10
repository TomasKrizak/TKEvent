// TK headers
#include "TKpoint.h"

using namespace std;

ClassImp(TKpoint);


TKpoint::TKpoint() :
	x(0.0), y(0.0), z(0.0)
{
}

TKpoint::TKpoint(double _x, double _y, double _z) : 
	x(_x), y(_y), z(_z)
{
}

TKpoint::~TKpoint()
{

}

void TKpoint::set_x(double _x)
{
	x = _x;
}

void TKpoint::set_y(double _y)
{
	y = _y;
}

void TKpoint::set_z(double _z)
{
	z = _z;
}

double TKpoint::get_x()
{
	return x;
}

double TKpoint::get_y()
{
	return y;
}

double TKpoint::get_z()
{
	return z;
}

void TKpoint::print()
{
	std::cout << "point (" << x << ", " << y << ", " << z << ")" << std::endl;
}

double distance_2D(TKpoint &point1, TKpoint &point2)
{
	double temp = pow(point1.get_x() - point2.get_x(), 2) + pow(point1.get_y() - point2.get_y(), 2);
	return sqrt(temp);
}

double distance_3D(TKpoint &point1, TKpoint &point2)
{
	double temp = pow(point1.get_x() - point2.get_x(), 2) + pow(point1.get_y() - point2.get_y(), 2) + pow(point1.get_z() - point2.get_z(), 2);
	return sqrt(temp);
}
