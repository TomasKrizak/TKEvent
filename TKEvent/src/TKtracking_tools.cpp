// TK headers
#include "TKEvent.h"

using namespace std;

void TKEvent::draw_sinusoids()
{
	gROOT->SetBatch(kTRUE);
	
	for(int side = 0; side < 2; side++)
	{
		vector<TKtrhit*> hits = filter_usable( filter_side(tr_hits, side) );
		
		if( hits.size() < 1 ) continue;
		
		TCanvas* canvas = new TCanvas("sinusoids");
		for(int hit = 0; hit < hits.size(); hit++)
		{
			double xi = hits[hit]->get_xy('x');
			double yi = hits[hit]->get_xy('y');
			auto function = new TF1("function", "[0]*sin(x)-[1]*cos(x)", 0.0, 2.0*M_PI);
			function->SetParameter(0, xi);
			function->SetParameter(1, yi);
			function->SetLineWidth(1);
			
			if(hit != 0)
			{
				function->Draw("Same");
			}
			else
			{
				function->Draw();
			}
		}
		canvas->SaveAs(Form("Events_visu/sinusoids-run-%d_event-%d_side-%d.png", run_number, event_number, side));
		canvas->Close();		
	}
}

void TKEvent::draw_likelihood()
{
	gROOT->SetBatch(kTRUE);
	
	const int resolution = 1250;
	const int iterations = 3;
	const double sigma = 2.0;

	for(int side = 0; side < 2; side++)
	{
		vector<TKtrhit*> hits = filter_usable( filter_side(tr_hits, side) );

		if( hits.size() < 3 ) continue;
	
		double phi1 = 0.0;
		double phi2 = 2.0*M_PI;
		double R1 = 0.0;
		double R2 = 2500.0;
		
		double peak_Theta;
		double peak_R;
		double delta_phi, delta_R;
		
		for(int iter = 0; iter < iterations; iter++)
		{	
			double offset_phi = (phi2 - phi1)/(2.0*resolution);
			double offset_R = (R2 - R1)/(2.0*resolution);
			TH2F *sinogram = new TH2F("sinogram", "sinogram; phi[rad]; r[mm]", resolution, phi1 + offset_phi, phi2 + offset_phi, resolution, R1 + offset_R, R2 + offset_R);		
				
			double norm_const = 1.0 / (sqrt(2.0*M_PI)*sigma); 
				
			for(int hit = 0; hit < hits.size(); hit++)
			{
				double R0 = hits[hit]->get_r();
				
				double weight;
				double r, theta;
				for(int j = 0; j <= resolution; j++)
				{
					theta = phi1 + (j * (phi2 - phi1) / resolution);
					double rho = hits[hit]->get_xy('x')*sin(theta) - hits[hit]->get_xy('y')*cos(theta);
					
					for(int k = 0; k <= resolution; k++)
					{
						r = R1 + (k * (R2 - R1) / resolution);
						
						weight = exp( -pow( abs(r-rho) - hits[hit]->get_r() , 2 ) / (2.0*sigma*sigma) );
						weight = weight * norm_const;
						
						if(hit != 0)
						{ 	
							weight = weight * sinogram->GetBinContent(j, k);
						}
						sinogram->SetBinContent(j, k, weight );
					}
				}	
			}
							
			double maximum = 0.0;
			for(int i = 1; i <= resolution; i++)
			{
				for(int j = 1; j <= resolution; j++)
				{
					if(maximum < sinogram->GetBinContent(i,j))
					{
						maximum = sinogram->GetBinContent(i,j);
						peak_Theta = i;
						peak_R = j;
					}
				}
			}
			
			
			if(maximum == 0.0) 
			{
				delete sinogram;			
				break;
			}
			
			if( true ) 
			{
				TCanvas* c2 = new TCanvas("sinogram","sinogram", 2000, 1600);
				c2->SetRightMargin(0.15);
				sinogram->SetStats(0);
				sinogram->SetContour(100);
				sinogram->Draw("COLZ");
				//gStyle->SetPalette(62);
				
				
				// Hough transform of anode wires
				if( 0 )
				{
					for(int hit = 0; hit < hits.size(); hit++)
					{
						double xi = hits[hit]->get_xy('x');
						double yi = hits[hit]->get_xy('y');
						auto function = new TF1("function", "[0]*sin(x)-[1]*cos(x)", phi1, phi2);
						function->SetParameter(0, xi);
						function->SetParameter(1, yi);
						function->SetLineWidth(1);
						function->Draw("Same");
					}
				}
				

				c2->SaveAs(Form("Events_visu/sinogram-run-%d_event-%d_side-%d_zoom-%d.png", run_number, event_number, side, iter));
				c2->Close();
			}
						
			delete sinogram;
			
			delta_phi = phi2 - phi1;
			delta_R = R2 - R1;
			
			peak_Theta = phi1 + delta_phi * (peak_Theta-0.5) / double(resolution);
			peak_R = R1 + delta_R * (peak_R-0.5) / double(resolution);
			
						
			phi1 = peak_Theta - 0.05*delta_phi;
			phi2 = peak_Theta + 0.05*delta_phi;
			R1 = peak_R - 0.05*delta_R;
			R2 = peak_R + 0.05*delta_R;

			cout << "Event: " << event_number << "	r: " << peak_R << "	phi: " << peak_Theta << endl;
			
			double a = tan(peak_Theta);
			double b = -peak_R / cos(peak_Theta);
			
			double joined_probability = 1.0;
			double probability;
			double norm = 1.0 / (sqrt(2.0*M_PI)*sigma);
			for(int i = 0; i < hits.size();i++)
			{	
				probability = exp(-pow(abs( peak_R + hits[i]->get_xy('x')*cos(peak_Theta) - hits[i]->get_xy('y')*sin(peak_Theta) ) - hits[i]->get_r() , 2)/(2.0*sigma*sigma));
				joined_probability = joined_probability * probability; 
			}		
			
			joined_probability = pow(norm, hits.size()) * joined_probability;
			
			cout << "likelihood R: " << joined_probability << endl;
		}
	}
}

void TKEvent::hough_transform(std::vector<TKtrhit*> hits, double phi_min, double phi_max, double R_min, double R_max, int ID)
{
	gROOT->SetBatch(kTRUE);

	const int resolution = 1000;
	double r, theta;
	
	const double sigma = 2.0;

	double S_r = 0.0;
	double S_rX = 0.0;
	double S_rY = 0.0;
	for(int i = 0; i < hits.size();i++)
	{
		S_r = S_r + 1.0/(sigma*sigma);
		S_rX = S_rX + hits[i]->get_xy('x')/(sigma*sigma);
		S_rY = S_rY + hits[i]->get_xy('y')/(sigma*sigma);
	}
					
	double peak_Theta;
	double peak_R;
	double delta_phi, delta_R;

	double offset1 = (phi_max - phi_min)/(2.0*resolution);
	double offset2 = (R_max - R_min)/(2.0*resolution);
	TCanvas *canvas = new TCanvas("canvas", "canvas", 1000, 1000);
	TH2F *hough = new TH2F("hough", "hough; theta; r", resolution, phi_min+offset1, phi_max+offset1, resolution, R_min+offset2, R_max+offset2);
	
	double boundary_down_1;
	double boundary_up_1;
	double boundary_down_2;		
	double boundary_up_2;

	double shift;
	
	for(int i = 0; i < hits.size(); i++)
	{	
		double x = hits.at(i)->get_xy('x');
		double y = hits.at(i)->get_xy('y');
		for(int k = 0; k < resolution; k++)
		{
			theta = phi_min + k * (phi_max - phi_min) / double(resolution);			
			shift = (S_rX/S_r)*sin(theta) - (S_rY/S_r)*cos(theta);
				
			boundary_down_1 = (x-22.0) * sin(theta) - (y+22.0) * cos(theta);
			boundary_up_1   = (x+22.0) * sin(theta) - (y-22.0) * cos(theta);
			boundary_down_2 = (x-22.0) * sin(theta) - (y-22.0) * cos(theta);
			boundary_up_2   = (x+22.0) * sin(theta) - (y+22.0) * cos(theta);
			
			
			for(int l = 0; l < resolution; l++)
			{	
				r = R_min + l * (R_max - R_min) / double(resolution);
				r = r + shift;
				if(theta <= M_PI/2.0)
				{
					if(boundary_down_1 <= r && r <= boundary_up_1  )
					{
						hough->SetBinContent(k, l,  hough->GetBinContent(k, l)+1);				
					}
				}
				else
				{
					if(boundary_down_2 <= r && r <= boundary_up_2  )
					{
						hough->SetBinContent(k, l,  hough->GetBinContent(k, l)+1);
					}				
				}	
			}
		}	
	}
	
	int side = hits.at(0)->get_SRL('S');
	hough->SetStats(0);
	//gStyle->SetPalette(62);
	hough->Draw("COLZ");
	// Hough transform of anode wires
	if( true )
	{
		for(int hit = 0; hit < hits.size(); hit++)
		{
			double xi = hits[hit]->get_xy('x');
			double yi = hits[hit]->get_xy('y');
			auto function = new TF1("function", "[0]*sin(x)-[1]*cos(x)", phi_min, phi_max);
			function->SetParameter(0, xi-(S_rX/S_r));
			function->SetParameter(1, yi-(S_rY/S_r));
			function->SetLineWidth(2);
			function->Draw("Same");
		}
	}
	canvas->SaveAs(Form("Events_visu/Hough_transform-run-%d_event-%d_side-%d_%d.png", run_number, event_number, side, ID));
	canvas->Close();
	delete hough;
}

void TKEvent::draw_likelihood_centred()
{
	gROOT->SetBatch(kTRUE);
	
	const int resolution = 1000;
	const int iterations = 3;
	const double sigma = 2.0;

	for(int side = 0; side < 2; side++)
	{
		vector<TKtrhit*> hits = filter_usable( filter_side(tr_hits, side) );
		if( hits.size() < 3 ) continue;

		double S_r = 0.0;
		double S_rX = 0.0;
		double S_rY = 0.0;
		for(int i = 0; i < hits.size();i++)
		{
			S_r = S_r + 1.0/(sigma*sigma);
			S_rX = S_rX + hits[i]->get_xy('x')/(sigma*sigma);
			S_rY = S_rY + hits[i]->get_xy('y')/(sigma*sigma);
		}
		
		double phi1 = 0.0;
		double phi2 = M_PI;
		double R1 = -150.0;
		double R2 = 150.0;
		
		double peak_Theta;
		double peak_R;
		
		for(int iter = 0; iter < iterations; iter++)
		{	
			double delta_phi = phi2 - phi1;
			double delta_R = R2 - R1;
			
			double offset_phi = delta_phi/(2.0*resolution);
			double offset_R = delta_R/(2.0*resolution);
			
			TH2F *sinogram_centred = new TH2F("sinogram_centred", "sinogram_centred; phi[rad]; r[mm]", resolution, phi1 + offset_phi, phi2 + offset_phi, resolution, R1 + offset_R, R2 + offset_R);	
			double norm_const = 1.0 / (sqrt(2.0*M_PI)*sigma); 
			
			double r, theta;
			for(int j = 0; j <= resolution; j++)
			{
				theta = phi1 + (j * delta_phi / resolution);
					
				for(int k = 0; k <= resolution; k++)
				{
					r = R1 + (k * delta_R / resolution);
					// shifting origin to "center of mass" of hits
					r = r + (S_rX/S_r)*sin(theta) - (S_rY/S_r)*cos(theta);
					
					double weight = 1.0;
					double temp;
					for(int hit = 0; hit < hits.size(); hit++)
					{
						double rho = hits[hit]->get_xy('x')*sin(theta) - hits[hit]->get_xy('y')*cos(theta);
						double R0 = hits[hit]->get_r();
					
						temp = exp( -pow( abs(r-rho) - hits[hit]->get_r() , 2 ) / (2.0*sigma*sigma) );
						temp = temp * norm_const;
						
						weight = weight * temp;
						
					}
					sinogram_centred->SetBinContent(j, k, log(weight) );
				}	
			}
							
			double maximum = -10e100;
			for(int i = 1; i <= resolution; i++)
			{
				for(int j = 1; j <= resolution; j++)
				{
					if(maximum < sinogram_centred->GetBinContent(i,j))
					{
						maximum = sinogram_centred->GetBinContent(i,j);
						peak_Theta = i;
						peak_R = j;
					}
				}
			}
			if(maximum == -10e100) 
			{
				delete sinogram_centred;			
				break;
			}
			
			
			if( true ) 
			{
				TCanvas* c2 = new TCanvas("sinogram_centred","sinogram_centred", 2000, 1600);
				c2->SetRightMargin(0.15);
				sinogram_centred->SetStats(0);
				sinogram_centred->SetContour(100);
				sinogram_centred->Draw("COLZ");
				//gStyle->SetPalette(62);
				
				// Hough transform of anode wires
				if( 1 )
				{
					for(int hit = 0; hit < hits.size(); hit++)
					{
						double xi = hits[hit]->get_xy('x');
						double yi = hits[hit]->get_xy('y');
						auto function = new TF1("function", "[0]*sin(x)-[1]*cos(x)", phi1, phi2);
						function->SetParameter(0, xi-(S_rX/S_r));
						function->SetParameter(1, yi-(S_rY/S_r));
						function->SetLineWidth(2);
						function->Draw("Same");
					}
				}
				
				
				c2->SaveAs(Form("Events_visu/sinogram_centred-run-%d_event-%d_side-%d_zoom-%d.png", run_number, event_number, side, iter));
				c2->Close();
			}
			delete sinogram_centred;
			
			peak_Theta = phi1 + delta_phi * (peak_Theta-0.5) / double(resolution);
			peak_R = R1 + delta_R * (peak_R-0.5) / double(resolution);
			
			phi1 = peak_Theta - 0.06*delta_phi;
			phi2 = peak_Theta + 0.06*delta_phi;
			R1 = peak_R - 0.06*delta_R;
			R2 = peak_R + 0.06*delta_R;

			peak_R = peak_R + (S_rX/S_r)*sin(peak_Theta) - (S_rY/S_r)*cos(peak_Theta);
			cout << "Event: " << event_number << "	r: " << peak_R << "	phi: " << peak_Theta << endl;
			
			double joined_probability = 1.0;
			double probability;
			double norm = 1.0 / (sqrt(2.0*M_PI)*sigma);
			for(int i = 0; i < hits.size();i++)
			{	
				probability = exp(-pow(abs( peak_R + hits[i]->get_xy('x')*cos(peak_Theta) - hits[i]->get_xy('y')*sin(peak_Theta) ) - hits[i]->get_r() , 2)/(2.0*sigma*sigma));
				joined_probability = joined_probability * probability; 
			}		
			
			joined_probability = pow(norm, hits.size()) * joined_probability;
			
			cout << "likelihood R: " << joined_probability << endl;
			
			if(iter == iterations-1)
			{
				TKtrack* track = new TKtrack(side, peak_Theta, peak_R);
				track->set_likelihood_R(joined_probability);
				tracks.push_back(track);
			}
		}
	}
}


void TKEvent::reconstruct_ML(bool save_sinograms)
{
	for(int side = 0; side < 2; side++)
	{
		vector<TKtrhit*> hits = filter_usable( filter_side(tr_hits, side) );
		if( hits.size() < 3 ) continue;
		
		TKcluster* cluster = find_cluster(hits);
		if( cluster != nullptr )
		{
			clusters.push_back( cluster );
			cluster->detect_ambiguity_type();
			cluster->reconstruct_MLM( save_sinograms, run_number, event_number );
			cluster->reconstruct_ambiguity();
			//hough_transform(cluster->get_tr_hits(), cluster->get_phi_min(), cluster->get_phi_max(), -50.0, 50.0, 0);
		}
	}
}

TKcluster* TKEvent::find_cluster(std::vector<TKtrhit*> hits)
{
	//number of different values of phi among which the cluster is being searched for
	const int bins_phi = 100;

	for(int i = 0; i < hits.size(); i++)
	{
		if( hits[0]->get_SRL('s') != hits[i]->get_SRL('s'))
		{
			cout << "cluster not found - hits from both sides included" << endl;
			return nullptr;
		}
	}

	if( hits.size() < 3 ) 
	{
		cout << "cluster not found - too few usable hits" << endl;
		return nullptr;
	}
	
	int max_count[bins_phi] = {0};
	double argmax_R[bins_phi];
	int global_max = 0;
	for(int i = 0; i < bins_phi; i++)
	{
		vector<double> boundaries;
		vector<int> hit_count;
		double theta = i * M_PI / double(bins_phi);		
		double boundary_down, boundary_up;
		
		boundaries.push_back(-10000);
		hit_count.push_back(0);
		boundaries.push_back(10000);
		
		for(int j = 0; j < hits.size(); j++)
		{
			double x = hits.at(j)->get_xy('x');
			double y = hits.at(j)->get_xy('y');
			if(i < bins_phi/2)
			{
				boundary_down = (x-22.0) * sin(theta) - (y+22.0) * cos(theta);
				boundary_up   = (x+22.0) * sin(theta) - (y-22.0) * cos(theta);
			}
			else
			{		
				boundary_down = (x-22.0) * sin(theta) - (y-22.0) * cos(theta);
				boundary_up   = (x+22.0) * sin(theta) - (y+22.0) * cos(theta);	
			}
			
			int k = 0;
			while( boundary_down > boundaries.at(k) )
			{
				k++;
			} 
			if( boundary_down < boundaries.at(k) )
			{
				boundaries.insert(boundaries.begin()+k, boundary_down);
				hit_count.insert(hit_count.begin()+k, hit_count.at(k-1));
			}
			while( boundary_up > boundaries.at(k) )			
			{
				hit_count.at(k)++;
				k++;
			}
			if( boundary_up < boundaries.at(k) )
			{
				boundaries.insert(boundaries.begin()+k, boundary_up);
				hit_count.insert(hit_count.begin()+k, hit_count.at(k-1)-1);
			}

		}
		
		for(int j = 0; j < hit_count.size(); j++)
		{
			if(hit_count.at(j) > max_count[i])
			{
				max_count[i] = hit_count.at(j);
				argmax_R[i] = (boundaries.at(j) + boundaries.at(j+1))/2.0;
			}
		}

		if(max_count[i] > global_max)
			global_max = max_count[i];
	}

	int phi_bin_min = 0;
	int phi_bin_max = 0;	
	double phi_min;
	double phi_max;
	
	if(max_count[0] == global_max)
	{
		while(max_count[phi_bin_max] == global_max)
		{
			phi_bin_max++;
		}
		phi_bin_min = bins_phi-1;
		while(max_count[phi_bin_min] == global_max)
		{
			phi_bin_min--;
		}

		phi_max = phi_bin_max * M_PI / double(bins_phi);
		phi_min = phi_bin_min * M_PI / double(bins_phi) - M_PI;
	}
	else
	{
		while(max_count[phi_bin_min] < global_max)
		{
			phi_bin_min++;
		}		
		phi_bin_max = phi_bin_min;
		while(max_count[phi_bin_max] == global_max)
		{
			phi_bin_max++;
		}

		phi_min = (phi_bin_min-1.0) * M_PI / double(bins_phi);
		phi_max = (phi_bin_max) * M_PI / double(bins_phi);
	}

	double R_0 = argmax_R[phi_bin_max-1];
	double theta_0 = double(phi_bin_max-1) * M_PI / double(bins_phi);
	
	vector<TKtrhit*> cluster_hits;
	for(int j = 0; j < hits.size(); j++)
	{
		double boundary_down, boundary_up;
		double x = hits.at(j)->get_xy('x');
		double y = hits.at(j)->get_xy('y');
		if(phi_bin_max-1 < bins_phi/2)
		{
			boundary_down = (x-22.0) * sin(theta_0) - (y+22.0) * cos(theta_0);
			boundary_up   = (x+22.0) * sin(theta_0) - (y-22.0) * cos(theta_0);
		}
		else
		{		
			boundary_down = (x-22.0) * sin(theta_0) - (y-22.0) * cos(theta_0);
			boundary_up   = (x+22.0) * sin(theta_0) - (y+22.0) * cos(theta_0);	
		}
		if(boundary_down <= R_0 && R_0 <= boundary_up)
		{
			cluster_hits.push_back(hits.at(j));
		}
	}
	TKcluster* cluster = new TKcluster(cluster_hits, phi_min, phi_max);
	return cluster;
}

