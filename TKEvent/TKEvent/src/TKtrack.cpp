// TK headers
#include "TKtrack.h"

ClassImp(TKtrack);

TKtrack::TKtrack()
{
}

TKtrack::TKtrack(int _side, double _a, double _b, double _c, double _d)
{
	side = _side;
	a = _a;
	b = _b;
	c = _c;
	d = _d;
}

TKtrack::~TKtrack()
{
}

void TKtrack::set_side(double _side)
{
	side = _side;
}

void TKtrack::set_a(double _a)
{
	a = _a;
}

void TKtrack::set_b(double _b)
{
	b = _b;
}

void TKtrack::set_c(double _c)
{
	c = _c;
}

void TKtrack::set_d(double _d)
{
	d = _d;
}

int TKtrack::get_side()
{
	return side;
} 

double TKtrack::get_a()
{
	return a;
} 

double TKtrack::get_b()
{
	return b;
} 

double TKtrack::get_c()
{
	return c;
} 

double TKtrack::get_d()
{
	return d;
} 

void TKtrack::print()
{
	std::cout << "Track | side: " << side 
	     	  << ", a = " << a 
	     	  << ", b = " << b 
	     	  << ", c = " << c 
	     	  << ", d = " << d << std::endl;
}
