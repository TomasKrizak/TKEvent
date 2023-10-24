// TK headers
#include "TKtrack.h"
using namespace std;

ClassImp(TKtrack);

TKtrack::TKtrack()
{
	likelihood   = 0.0;
	likelihood_R = 0.0;
	likelihood_Z = 0.0;
	
	quality   = 0.0;
	quality_R = 0.0;
	quality_Z = 0.0;
}

TKtrack::TKtrack(int _side, double _phi, double _r)
{
	side  = _side;
	phi   = _phi;
	r     = _r;
	associated_tr_hits = std::vector<TKtrhit*>();
	
	a = tan(phi);
	b = -r / cos(phi);
	
	likelihood   = 0.0;
	likelihood_R = 0.0;
	likelihood_Z = 0.0;
	
	quality   = 0.0;
	quality_R = 0.0;
	quality_Z = 0.0;
}

TKtrack::TKtrack(int _side, double _a, double _b, double _c, double _d)
{
	side = _side;
	a = _a;
	b = _b;
	c = _c;
	d = _d;
	associated_tr_hits = std::vector<TKtrhit*>();
	
	phi = atan(a);
	r = -b/sqrt(a*a+1.0);
	theta = atan(c/sqrt(a*a+1.0));
	h = d - a*b*c/(a*a+1.0);
	
	likelihood   = 0.0;
	likelihood_R = 0.0;
	likelihood_Z = 0.0;
	
	quality   = 0.0;
	quality_R = 0.0;
	quality_Z = 0.0;
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

void TKtrack::set_phi(double _phi)
{
	phi = _phi;
}

void TKtrack::set_r(double _r)
{
	r = _r;
}

void TKtrack::set_theta(double _theta)
{
	theta = _theta;
}

void TKtrack::set_h(double _h)
{
	h = _h;
}

void TKtrack::set_chi_squared(double _chi_squared)
{
	chi_squared = _chi_squared;
}

void TKtrack::set_chi_squared_R(double _chi_squared_R)
{
	chi_squared_R = _chi_squared_R;
}

void TKtrack::set_chi_squared_Z(double _chi_squared_Z)
{
	chi_squared_Z = _chi_squared_Z;
}

void TKtrack::set_quality(double _quality)
{
	quality = _quality;
}

void TKtrack::set_quality_R(double _quality_R)
{
	quality_R = _quality_R;
}

void TKtrack::set_quality_Z(double _quality_Z)
{
	quality_Z = _quality_Z;
}

void TKtrack::set_likelihood(double _likelihood)
{
	likelihood = _likelihood;
}

void TKtrack::set_likelihood_R(double _likelihood_R)
{
	likelihood_R = _likelihood_R;
}

void TKtrack::set_likelihood_Z(double _likelihood_Z)
{
	likelihood_Z = _likelihood_Z;
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

double TKtrack::get_phi()
{
	return phi;
} 

double TKtrack::get_r()
{
	return r;
} 

double TKtrack::get_theta()
{
	return theta;
} 

double TKtrack::get_h()
{
	return h;
} 

double TKtrack::get_chi_squared()
{
	return chi_squared;
}

double TKtrack::get_chi_squared_R()
{
	return chi_squared_R;
}

double TKtrack::get_chi_squared_Z()
{
	return chi_squared_Z;
}

double TKtrack::get_quality()
{
	return quality;
} 

double TKtrack::get_quality_R()
{
	return quality_R;
} 

double TKtrack::get_quality_Z()
{
	return quality_Z;
} 

double TKtrack::get_likelihood()
{
	return likelihood;
} 

double TKtrack::get_likelihood_R()
{
	return likelihood_R;
} 

double TKtrack::get_likelihood_Z()
{
	return likelihood_Z;
} 

void TKtrack::reconstruct_vertical_least_square()
{
	std::vector<double> hits_p;
	double projection_distance;
	double distance_from_wire;
	double denominator = sqrt((a*a) + 1);	
	for(int i = 0; i < associated_tr_hits.size(); i++)
	{
		double x = associated_tr_hits[i]->get_xy('x');
		double y = associated_tr_hits[i]->get_xy('y');
		distance_from_wire = abs(y - a*x - b) / denominator;
		projection_distance = pow(x, 2.0) + pow(y - b, 2.0) - pow(distance_from_wire, 2.0);
		projection_distance = pow(projection_distance, 0.5);
		hits_p.push_back(projection_distance);
	} 
	
	int sum_n = 0;
	double sum_Z = 0.0;
	double sum_P = 0.0;
	double sum_ZP = 0.0;
	double sum_PP = 0.0;
	for(int i = 0; i < associated_tr_hits.size(); i++)
	{
		if(associated_tr_hits[i]->get_h() != 0.0)
		{
			sum_Z += associated_tr_hits[i]->get_h();
			sum_P += hits_p[i];
			sum_ZP += associated_tr_hits[i]->get_h()*hits_p[i];
			sum_PP += hits_p[i]*hits_p[i];
			sum_n++;		
		}
	}
	
	double c_temp;
	denominator = sum_n*sum_PP - sum_P*sum_P; 
	if(sum_n > 0 && denominator != 0.0)
	{
		c_temp = (sum_n*sum_ZP - sum_P*sum_Z) / denominator;
		
		c = c_temp*(2.0*side-1.0)*pow(1+a*a, 0.5);
		d = (sum_PP*sum_Z - sum_P*sum_ZP) / denominator;
	}	
	else
	{
		c = 0.0;
		d = 0.0;
	}
	
	theta = atan(c/sqrt(a*a+1.0));
	h = d - a*b*c/(a*a+1.0);
	
	this->update_likelihood();
}

void TKtrack::reconstruct_vertical_MLM()
{
	double S = 0.0;
	double Sx = 0.0;
	double Sy = 0.0;
	double Sz = 0.0;
	double Sxy = 0.0;
	double Sxz = 0.0;
	double Syz = 0.0;
	double Sxx = 0.0;
	double Syy = 0.0;
	double x,y,z;
	
	for(int i = 0; i < associated_tr_hits.size(); i++)
	{
		z = associated_tr_hits[i]->get_h();
		if(z != 0.0)
		{	
			x = associated_tr_hits[i]->get_xy('x');
			y = associated_tr_hits[i]->get_xy('y');
			
			S++;
			Sx = Sx + x;
			Sy = Sy + y;
			Sz = Sz + z;
			Sxy = Sxy + x*y;
			Sxz = Sxz + x*z;
			Syz = Syz + y*z;
			Sxx = Sxx + x*x;
			Syy = Syy + y*y;
		}
	}
	
	double denominator = (Sy*Sy-S*Syy)*sin(phi)*sin(phi) + 2.0*(Sx*Sy-S*Sxy)*sin(phi)*cos(phi) + (Sx*Sx-S*Sxx)*cos(phi)*cos(phi);
	if( denominator != 0.0)
	{
		double tan_theta = ((Sy*Sz-S*Syz)*sin(phi) + (Sx*Sz-S*Sxz)*cos(phi)) / denominator;
		double h_temp = Sz / S - (Sx*cos(phi) + Sy*sin(phi)) * tan_theta / S;
		
		h = h_temp;
		theta = atan(tan_theta);

		c = tan_theta/cos(phi);
		d = h - r*tan(phi)*tan_theta;
		
		this->update_likelihood();	
	}
	else 
	{
		h = 0.0;
		theta = 0.0;

		c = 0.0;
		d = 0.0;
	}
}

void TKtrack::update_likelihood()
{
	const double sigma_z = 17.0; 
		
	double x,y,z;
	
	double sin_phi = sin(phi);
	double cos_phi = cos(phi);
	double tan_theta = tan(theta); 
	
	double sum = 0.0;
	double temp;

	int used_hits = 0.0;
	for(int i = 0; i < associated_tr_hits.size(); i++)
	{
		z = associated_tr_hits[i]->get_h();
		if(z != 0.0)
		{	
			used_hits++;
			x = associated_tr_hits[i]->get_xy('x');
			y = associated_tr_hits[i]->get_xy('y');
			if(r > 0)
			{
				temp = h + tan_theta*(y*sin_phi + x*cos_phi) - z; 
			}			
			else
			{
				temp = h - tan_theta*(y*sin_phi + x*cos_phi) - z; 
			}
			
			sum = sum + pow(temp, 2.0);
		}		
	}
	if(used_hits > 0)
	{
	
		chi_squared_Z = sum/(sigma_z*sigma_z*double(used_hits));
		chi_squared = chi_squared_Z * chi_squared_R; 
	
		quality_Z = exp( -0.5*chi_squared_Z );
		quality = quality_R * quality_Z;
		
		double con = pow( 2.0*M_PI*sigma_z*sigma_z, -2.0*double(used_hits)); 

		likelihood_Z = con *  exp( -0.5*chi_squared_Z*double(used_hits) );
		likelihood = likelihood_R * likelihood_Z; 	
	}
}

void TKtrack::print()
{
	std::cout << "Track | side: " << side 
	     	  << ", a = " << a 
	     	  << ", b = " << b 
	     	  << ", c = " << c 
		  << ", d = " << d << std::endl; 
	std::cout << "		 phi = " << phi 
	     	  << ", r = " << r 
	     	  << ", theta = " << theta 
	     	  << ", h = " << h << std::endl; 
	std::cout << "	quality: " << quality << std::endl;
	std::cout << "	quality R: " << quality_R << std::endl;
	std::cout << "	quality Z: " << quality_Z << std::endl;
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


