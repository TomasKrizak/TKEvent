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
	associated_tr_hits = std::vector<TKtrhit*>();
}

TKtrack::~TKtrack()
{
	associated_tr_hits.clear();
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

void TKtrack::set_likelihood(double _likelihood)
{
	likelihood = _likelihood;
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

double TKtrack::get_likelihood()
{
	return likelihood;
} 

void TKtrack::print()
{
	std::cout << "Track | side: " << side 
	     	  << ", a = " << a 
	     	  << ", b = " << b 
	     	  << ", c = " << c 
	     	  << ", d = " << d << std::endl; 
	//std::cout << "	likelihood: " << likelihood << std::endl;
	std::cout << "	number of associated tracker hits: " << associated_tr_hits.size() << std::endl << std::endl;
}

void TKtrack::add_associated_tr_hit(TKtrhit* tracker_hit)
{
	associated_tr_hits.push_back(tracker_hit);
}

std::vector<TKtrhit*> TKtrack::get_associated_tr_hits()
{
	return associated_tr_hits;
}





