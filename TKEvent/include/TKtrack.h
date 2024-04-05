#ifndef _TKTRACK_H_
#define _TKTRACK_H_

// Standard headers
#include <iostream>
#include <vector>
#include <cmath>

// ROOT headers
#include "TObject.h"

#include "TKOMhit.h"
#include "TKtrhit.h"
#include "TKpoint.h"

// line in the form:
//	y = ax + b
//	z = cx + d

// parametrically:
//	x(s) = s
//	y(s) = a*s + b
//	z(s) = c*s + d

// or in the form:
// 	x(t) = cos(phi)*cos(theta)*t + r*sin(phi)
// 	y(t) = sin(phi)*cos(theta)*t - r*cos(phi)
// 	z(t) = sin(theta)*t + h

// set of transformations:
//	a = tan(phi)
//	b = -r/cos(phi)
//	c = tan(theta)/cos(phi)
//	d = h - r*tan(phi)*tan(theta)

//	phi = atan(a) 
//	r = -b/sqrt(a*a+1.0)
//	theta = atan(c/sqrt(a*a+1.0))
//	h = d - a*b*c/(a*a+1.0)

// associated tracker hit point: (tracker hit - anode wire (xi,yi), vertical position zi)
// 	x(t) = (xi*cos(phi)*+yi*sin(phi))*cos(phi) + r*sin(phi)
// 	y(t) = (xi*cos(phi)*+yi*sin(phi))*sin(phi) - r*cos(phi)
// 	z(t) = (xi*cos(phi)*+yi*sin(phi))*tan(theta) + h - zi

class TKtrack: public TObject
{
	private:
		
		int    side;
		double a;
		double b;
		double c;
		double d;
		
		double phi;
		double r;
		double theta;
		double h;
		
		// metrics for determining quality of track 
		
		double chi_squared;
		double chi_squared_R;
		double chi_squared_Z;
		
		// it is the n-th root of likelihood, scaled so that the theoretical maximum value is 1,
		//	where n is number of hits in the cluster 
		// value from [0,1] where 1 is absolutely perfect track
		// quality = quality_R * quality_Z 
		double quality;
		double quality_R;
		double quality_Z;
		
		// most reasonable statistical evaluation
		// L = L_R * L_Z
		double likelihood; 
		// horizontal part of likelihood
		double likelihood_R; 
		// vertical part of likelihood
		double likelihood_Z;
		
		// in case if ambiguities the mirror images are linked
		// ambiguity_type == 0 means no ambiguity
		int ambiguity_type;
		TKtrack* mirror_image;
		
		// association_distance can be changed in reconstruction functions (3 sigma by default = 6mm)
		std::vector<TKtrhit*> associated_tr_hits; 
		std::vector<TKpoint*> associated_tr_hit_points; 

	public:
		
		TKtrack();
		TKtrack(int _side, double _phi, double _r);
		TKtrack(int _side, double _a, double _b, double _c, double _d);
		~TKtrack();
		
		void add_associated_tr_hit(TKtrhit* tracker_hit);
		std::vector<TKtrhit*> get_associated_tr_hits();
		
		void add_associated_tr_hit_point(TKpoint* tracker_hit_point);
		std::vector<TKpoint*> get_associated_tr_hit_points();
		
		void calculate_tr_hit_points();
			
		void set_side(double _side);
		
		void set_a   (double _a);
		void set_b   (double _b);
		void set_c   (double _c);
		void set_d   (double _d);
		
		void set_phi   (double _phi);
		void set_r     (double _r);
		void set_theta (double _theta);
		void set_h     (double _h);
		
		void set_chi_squared(double _chi_squared);
		void set_chi_squared_R(double _chi_squared_R);
		void set_chi_squared_Z(double _chi_squared_Z);
		
		void set_quality(double _quality);
		void set_quality_R(double _quality_R);
		void set_quality_Z(double _quality_Z);
		
		void set_likelihood(double _likelihood);
		void set_likelihood_R(double _likelihood_R);
		void set_likelihood_Z(double _likelihood_Z);
		
		void set_ambiguity_type(int _ambiguity_type);
		
		// links a mirror solution to the original track and copies common information: 
		// 	associated hits, quality_R, likelihood_R, chi_squared_R
		//	note: link is not mutual - original track is a "master track"
		void link_mirror_image(TKtrack* _mirror_image);
		
		int get_side();
		
		double get_a();
		double get_b();
		double get_c();
		double get_d();
		
		double get_phi  ();
		double get_r    ();
		double get_theta();
		double get_h    ();
		
		double get_chi_squared();
		double get_chi_squared_R();
		double get_chi_squared_Z();
		
		double get_quality();
		double get_quality_R();
		double get_quality_Z();
		
		double get_likelihood();
		double get_likelihood_R();
		double get_likelihood_Z();

		int get_ambiguity_type();
		TKtrack* get_mirror_image();

		// calculates likelihood_Z and likelihood based on associated hits
		void update_likelihood();
		void reconstruct_vertical_least_square();
		void reconstruct_vertical_MLM();

		void print();
		
		ClassDef(TKtrack,1);
};

#endif
