// TK headers
#include "TKcluster.h"

ClassImp(TKcluster);


TKcluster::TKcluster()
{
	ambiguity_type = 0;
	phi_min = 0.0;
	phi_max = M_PI;
	std::vector<TKtrhit*> cluster_tr_hits;
	track = nullptr;
}

TKcluster::TKcluster(std::vector<TKtrhit*> tr_hits, double _phi_min, double _phi_max)
{
	ambiguity_type = 0;
	side = tr_hits.at(0)->get_SRL('s');
	phi_min = _phi_min;
	phi_max = _phi_max;
	cluster_tr_hits = std::vector<TKtrhit*>();
	for(int i = 0; i < tr_hits.size(); i++)
	{
		cluster_tr_hits.push_back( tr_hits.at(i) );
	}
	track = nullptr;
}

TKcluster::~TKcluster()
{
	cluster_tr_hits.clear();
}

void TKcluster::add_tr_hit(TKtrhit* tracker_hit)
{
	cluster_tr_hits.push_back( tracker_hit );
}

std::vector<TKtrhit*> TKcluster::get_tr_hits()
{
	return cluster_tr_hits;
}

void TKcluster::set_track(TKtrack *_track)
{
	track = _track;
}

TKtrack* TKcluster::get_track()
{
	return track;
}

void TKcluster::set_side(double _side)
{
	side = _side;
}

void TKcluster::set_phi_min(double _phi_min)
{
	phi_min = _phi_min;
}

void TKcluster::set_phi_max(double _phi_max)
{
	phi_max = _phi_max;
}

int    TKcluster::get_side()
{
	return side;
}

double TKcluster::get_phi_min()
{
	return phi_min;
}

double TKcluster::get_phi_max()
{
	return phi_max;
}

int TKcluster::get_ambiguity_type()
{
	return ambiguity_type;
}

void TKcluster::detect_ambiguity_type()
{
	int ambiguity = 0;
	double x0 = cluster_tr_hits[0]->get_xy('x');
	double y0 = cluster_tr_hits[0]->get_xy('y');
	bool ambiguous;
	
	// type 1 == mirror image along line x = x0 
	ambiguous = true;
	for(int i = 1; i < cluster_tr_hits.size(); i++)
	{	
		if( x0 != cluster_tr_hits[i]->get_xy('x') )
		{
			ambiguous = false;
			break;
		}
	}
	if( ambiguous == true )
	{
		ambiguity = 1;
	}
	
	// type 2 == mirror image along line y = y0 
	ambiguous = true;
	for(int i = 1; i < cluster_tr_hits.size(); i++)
	{	
		if( y0 != cluster_tr_hits[i]->get_xy('y') )
		{
			ambiguous = false;
			break;
		}
	}
	if( ambiguous == true )
	{
		ambiguity = 2;
	}
	
	// type 3 == mirror image along line y = x + (y0-x0) 
	ambiguous = true;
	for(int i = 1; i < cluster_tr_hits.size(); i++)
	{	
		if( y0-cluster_tr_hits[i]->get_xy('y') != x0-cluster_tr_hits[i]->get_xy('x') )
		{
			ambiguous = false;
			break;
		}
	}
	if( ambiguous == true )
	{
		ambiguity = 3;
	}
	
	// type 4 == mirror image along line y = -x + (y0-x0) 
	ambiguous = true;
	for(int i = 1; i < cluster_tr_hits.size(); i++)
	{	
		if( y0-cluster_tr_hits[i]->get_xy('y') != cluster_tr_hits[i]->get_xy('x')-x0 )
		{
			ambiguous = false;
			break;
		}
	}
	if( ambiguous == true )
	{
		ambiguity = 4;
	}

	ambiguity_type = ambiguity;
}


void TKcluster::reconstruct_ambiguity()
{
	if( ambiguity_type != 0 )
	{
		double a0 = track->get_a();
		double b0 = track->get_b();
		double x0 = cluster_tr_hits[0]->get_xy('x');
		double y0 = cluster_tr_hits[0]->get_xy('y');
		double a, b;
		
		switch(ambiguity_type)
		{
			case 1:
				a = -a0;
				b = 2.0*a0*x0 + b0;
				break;
			case 2:
				a = -a0;
				b = 2.0*y0 - b0;
				break;
			case 3:
				a = 1.0/a0;
				b = y0-x0 - b0/a0 + (y0-x0)/a0;
				break;
			case 4:
				a = 1.0/a0;
				b = y0+x0 + b0/a0 - (y0+x0)/a0;
				break;
		}
		
		TKtrack* mirror_image = new TKtrack(side, a, b, 0.0, 0.0);
		track->set_ambiguity_type( ambiguity_type );
		track->link_mirror_image( mirror_image );
		mirror_image->reconstruct_vertical_MLM();
	
	}
}


void TKcluster::reconstruct_MLM(bool save_sinograms, int run_number, int event_number)
{
	gROOT->SetBatch(kTRUE);
	
	const int resolution = 70;
	const int iterations = 2;
	const double sigma = 2.0;

	double S_r = 0.0;
	double S_rX = 0.0;
	double S_rY = 0.0;
	for(int i = 0; i < cluster_tr_hits.size();i++)
	{
		S_r = S_r + 1.0/(sigma*sigma);
		S_rX = S_rX + cluster_tr_hits[i]->get_xy('x')/(sigma*sigma);
		S_rY = S_rY + cluster_tr_hits[i]->get_xy('y')/(sigma*sigma);
	}
	
	double phi1 = phi_min;
	double phi2 = phi_max;
	double R1 = -35.0;
	double R2 = 35.0;
	
	double peak_Theta;
	double peak_R;
	
	double minimum;
	for(int iter = 0; iter < iterations; iter++)
	{	
		minimum = 1e100;
		double delta_phi = phi2 - phi1;
		double delta_R = R2 - R1;
		
		double offset_phi = delta_phi/(2.0*resolution);
		double offset_R = delta_R/(2.0*resolution);
		
		TH2F *chi_squared = new TH2F("chi_squared", "chi_squared; phi[rad]; r[mm]", resolution, phi1 + offset_phi, phi2 + offset_phi, resolution, R1 + offset_R, R2 + offset_R);	
				
		//double norm_const = pow(2.0*M_PI*sigma*sigma, -1.0*hits.size()/2.0); 
		double weight;
		double r, theta;
		double sin_theta, cos_theta;
		double R0;
		double rho;
		double shift;
		double diff_R;
		
		for(int j = 1; j <= resolution; j++)
		{
			theta = phi1 + (j * delta_phi / resolution);
			sin_theta = sin(theta);
			cos_theta = cos(theta);
			
			// shifting origin to "center of mass" of hits			
			shift = (S_rX/S_r)*sin_theta - (S_rY/S_r)*cos_theta;
			for(int k = 1; k <= resolution; k++)
			{
				r = R1 + (k * delta_R / resolution);
				r = r + shift;
				double diff_R_sum = 0.0;
				for(int hit = 0; hit < cluster_tr_hits.size(); hit++)
				{
					rho = cluster_tr_hits[hit]->get_xy('x')*sin_theta - cluster_tr_hits[hit]->get_xy('y')*cos_theta;
					R0 = cluster_tr_hits[hit]->get_r();
					diff_R = pow( abs(r-rho) - R0, 2 );
					diff_R_sum = diff_R_sum + diff_R;
				}
				
				//weight = -diff_R_sum / (2.0*sigma*sigma); 
				//weight = exp( -diff_R_sum / (2.0*sigma*sigma*double(cluster_tr_hits.size())) );
				//weight = weight * norm_const;
				chi_squared->SetBinContent(j, k, diff_R_sum );
			}	
		}
		
		
		for(int i = 1; i <= resolution; i++)
		{
			for(int j = 1; j <= resolution; j++)
			{
				if(minimum > chi_squared->GetBinContent(i,j))
				{
					minimum = chi_squared->GetBinContent(i,j);
					peak_Theta = i;
					peak_R = j;
				}
			}
		}
		
		if(minimum == 1e100) 
		{
			delete chi_squared;			
			break;
		}
		
		if( save_sinograms == true ) 
		{
			TCanvas* c2 = new TCanvas("chi_squared","chi_squared", 1000, 1000);
			c2->SetRightMargin(0.15);
			chi_squared->SetStats(0);
			chi_squared->SetContour(100);
			chi_squared->Draw("COLZ");
			//gStyle->SetPalette(62);
			
			// Hough transform of anode wires
			if( 1 )
			{
				for(int hit = 0; hit < cluster_tr_hits.size(); hit++)
				{
					double xi = cluster_tr_hits[hit]->get_xy('x');
					double yi = cluster_tr_hits[hit]->get_xy('y');
					auto function = new TF1("function", "[0]*sin(x)-[1]*cos(x)", phi1, phi2);
					function->SetParameter(0, xi-(S_rX/S_r));
					function->SetParameter(1, yi-(S_rY/S_r));
					function->SetLineWidth(2);
					function->Draw("Same");
				}
			}
			
			c2->SaveAs(Form("Events_visu/chi_squared-run-%d_event-%d_side-%d_zoom-%d.png", run_number, event_number, side, iter));
			c2->Close();
		}
		delete chi_squared;
		
		peak_Theta = phi1 + delta_phi * (peak_Theta-0.5) / double(resolution);
		peak_R = R1 + delta_R * (peak_R-0.5) / double(resolution);
		
		phi1 = peak_Theta - 0.05*delta_phi;
		phi2 = peak_Theta + 0.05*delta_phi;
		R1 = peak_R - 0.05*delta_R;
		R2 = peak_R + 0.05*delta_R;

		peak_R = peak_R + (S_rX/S_r)*sin(peak_Theta) - (S_rY/S_r)*cos(peak_Theta);
	}
			
	if(minimum != 1e100)
	{	
		/*
		if(peak_Theta < 0.0)
		{
			peak_Theta = peak_Theta + M_PI;
			peak_R = -peak_R;
		}*/
		
		track = new TKtrack(side, peak_Theta, peak_R);
				
		double a = track->get_a();
		double b = track->get_b();

		double association_distance = 3.0*sigma;
		double denominator = sqrt((a*a) + 1);	
		double distance_from_wire;
		for(int i = 0; i < cluster_tr_hits.size(); i++)
		{
			distance_from_wire = abs(cluster_tr_hits[i]->get_xy('y') - a*cluster_tr_hits[i]->get_xy('x') - b) / denominator;
			if( abs(distance_from_wire - cluster_tr_hits[i]->get_r()) < association_distance )
			{
				track->add_associated_tr_hit(cluster_tr_hits[i]);
				cluster_tr_hits[i]->set_associated_track(track);
			}
		}		
		
		double chi_squared = minimum/(double(cluster_tr_hits.size())*sigma*sigma);
		double norm_const = pow( 2.0*M_PI*sigma*sigma, -1.0*double(cluster_tr_hits.size())/2.0 ); 
		
		track->set_chi_squared_R( chi_squared );
		track->set_quality_R( exp(-0.5*chi_squared) );		
		track->set_likelihood_R( norm_const * exp(-0.5*chi_squared*double(cluster_tr_hits.size())) );
			
		//track->reconstruct_vertical_least_square();
		track->reconstruct_vertical_MLM();
	}
}

void TKcluster::print()
{
	std::cout << "Cluster | side "       << side
		  << ", size "       << cluster_tr_hits.size()
	     	  << ", phi_min: " << phi_min 
	     	  << ", phi_max: "   << phi_max 
	     	  << ", ambiguity type: " << ambiguity_type << std::endl;
}
