// TK headers
#include "TKEvent.h"

using namespace std;

bool debug_mode = false;

void TKEvent::reconstruct_track_from_hits(std::vector<TKtrhit*>& hits, bool save_sinograms)
{
	// ancient algorithm - does not provide complete likelihood and (r,phi,theta,h) coordinates

	// finds only one track per side of the tracker seperately from given set of tracker hits
	// no uncertainties involved here
	
	// each iteration zooms into a region with 10% size in both directions
	const int resolution = 250;
	const int iterations = 2;

	for(int side = 0; side < 2; side++)
	{	
		vector<double> hits_x;
		vector<double> hits_y;
		vector<double> hits_z;
		vector<double> hits_r;
		vector<double> hits_sigma;	
		vector<int>    index;
		
		for(int i = 0; i < hits.size(); i++)
		{
			if( side == hits[i]->get_SRL('s'))
			{
				// not using broken or (wrongly associated) tracker hits with too large radius
				if( hits[i]->get_r() != -1.0 && hits[i]->get_r() < 35.0 )
				{
					hits_x.push_back( hits[i]->get_xy('x') );
					hits_y.push_back( hits[i]->get_xy('y') );
					hits_z.push_back( hits[i]->get_h() );
					hits_r.push_back( hits[i]->get_r() );
					hits_sigma.push_back( hits[i]->get_sigma_R() );
					index.push_back( i );
				}
			}
		}

		if( hits_x.size() < 3 ) continue;

		// fitting top projections of tracks using Legender transform
		double r, theta;
		double phi1 = 0.0;
		double phi2 = M_PI;
		double R1 = -2500;
		double R2 = 2500;
		
		double peak_Theta;
		double peak_R;
		double delta_phi, delta_R;
		
		for(int q = 0; q < iterations; q++)
		{	
			double offset = (phi2 - phi1)/(2.0*resolution);
			TH2F *sinograms = new TH2F("sinograms", "sinograms; theta; r", resolution, phi1+offset, phi2+offset, resolution, R1, R2);		
			for(int i = 0; i < hits_x.size(); i++)
			{
				for(int k = 0; k <= resolution; k++)
				{
					theta = phi1 + (k * (phi2 - phi1) / resolution);
					r = hits_x[i]*sin(theta) - hits_y[i]*cos(theta);
					sinograms->Fill( theta, r - hits_r[i] );
					sinograms->Fill( theta, r + hits_r[i] );
				}	
			}
								
			double maximum = 0.0;
			for(int i = 1; i <= resolution; i++)
			{
				for(int j = 1; j <= resolution; j++)
				{
					if(maximum < sinograms->GetBinContent(i,j))
					{
						maximum = sinograms->GetBinContent(i,j);
						peak_Theta = i;
						peak_R = j;
					}
				}
			}
			
			delta_phi = phi2 - phi1;
			delta_R = R2 - R1;
			
			peak_Theta = phi1 + delta_phi * peak_Theta / double(resolution);
			peak_R = R1 + delta_R * (peak_R-0.5) / double(resolution);
			
			phi1 = peak_Theta - 0.05*delta_phi;
			phi2 = peak_Theta + 0.05*delta_phi;
			R1 = peak_R - 0.05*delta_R;
			R2 = peak_R + 0.05*delta_R;

			if( save_sinograms == true ) 
			{
				TCanvas* c2 = new TCanvas("sinograms");
				sinograms->SetContour(100);
				sinograms->SetStats(0);
				sinograms->Draw("COLZ");
				c2->SaveAs(Form("Events_visu/sinograms-run-%d_event-%d_side-%d_zoom-%d.png", run_number, event_number, side, q));
				c2->Close();
			}
			delete sinograms;
		}
		
		// creating track
		TKtrack* track = new TKtrack(side, peak_Theta, peak_R);	
		
		double a = track->get_a();
		double b = track->get_b();
		
		// associating hits to track
		vector<bool> hits_associated;
		double denominator = sqrt((a*a) + 1);	
		for(int i = 0; i < hits_x.size(); i++)
		{
			double distance_from_wire = abs(hits_y[i] - a*hits_x[i] - b) / denominator;
			double association_distance = 3.0 * hits_sigma[i];
			if( abs(distance_from_wire - hits_r[i]) < association_distance )
			{
				hits_associated.push_back(true);
				track->add_associated_tr_hit(hits[index[i]]);
				hits[index[i]]->set_associated_track(track);
			}
			else
			{
				hits_associated.push_back(false);
			}
		
		}					
	
		// calculating likelihood_R
		double joined_probability = 1.0;
		double probability;
		
		double norm = 1.0 / sqrt(2.0*M_PI);
		for(int i = 0; i < hits.size();i++)
		{	
			
			probability = exp(-pow(abs( r + hits_x[i]*cos(theta) - hits_y[i]*sin(theta) ) - hits_r[i] , 2)/(2.0*hits_sigma[i]*hits_sigma[i]));
			joined_probability = joined_probability * probability / hits_sigma[i]; 
		}		
		
		joined_probability = pow(norm, hits.size()) * joined_probability;
		track->set_likelihood_R(joined_probability);
		
		// fitting z coordinates of hits
		vector<double> hits_p;
		for(int i = 0; i < hits_x.size(); i++)
		{
			double distance_from_wire = abs(hits_y[i] - a*hits_x[i] - b) / denominator;
			double projection_distance = pow(hits_x[i], 2.0) + pow((hits_y[i] - b), 2.0) - pow(distance_from_wire, 2.0);
			projection_distance = pow(projection_distance, 0.5);
			
			hits_p.push_back(projection_distance);
		} 
		
		int sum_n = 0;
		double sum_Z = 0.0;
		double sum_P = 0.0;
		double sum_ZP = 0.0;
		double sum_PP = 0.0;
		for(int i = 0; i < hits_x.size(); i++)
		{
			if(hits_z[i] != 0.0 && hits_associated[i] == true)
			{
				sum_Z += hits_z[i];
				sum_P += hits_p[i];
				sum_ZP += hits_z[i]*hits_p[i];
				sum_PP += hits_p[i]*hits_p[i];
				sum_n++;		
			}
		}
		
		double c = 0.0;
		double d = 0.0;
		denominator = sum_n*sum_PP - sum_P*sum_P; 
		if(sum_n > 0 && denominator != 0.0)
		{
			c = (sum_n*sum_ZP - sum_P*sum_Z) / denominator;
			c = c*(2.0*side-1.0)*pow(1+a*a, 0.5);
			d = (sum_PP*sum_Z - sum_P*sum_ZP) / denominator;
		}
		
		track->set_c(c);
		track->set_d(d);
		
		track->set_theta( atan(c/sqrt(a*a+1.0)) ); 
		track->set_h( d - a*b*c/(a*a+1.0) );	
		
		track->update_likelihood();
				
		tracks.push_back(track);
		
		hits_x.erase(hits_x.begin(), hits_x.end());
		hits_y.erase(hits_y.begin(), hits_y.end());
		hits_z.erase(hits_z.begin(), hits_z.end());
		hits_r.erase(hits_r.begin(), hits_r.end());
	}
}

void TKEvent::reconstruct_track(bool save_sinograms)
{
	reconstruct_track_from_hits(tr_hits, save_sinograms);
}

void TKEvent::reconstruct_multi(bool save_sinograms)
{
	gROOT->SetBatch(kTRUE);

	// resolution of both dimensions of every 2D histogram	
	int resolution = 50; 
	
	// number of zooming iterations into histogram
	const int iterations = 3;
	
	// how high does a peak has to be in comparison to the highest one in order to be considered as a candidate as well
	const double treshold = 0.75;
	
	// into how many segments (in both directions) is histogram divided during each iteration (except the final one)
	// new peak is calculated on each segment 
	const int start_no_segments = 5; 
	
	// both sides are dealt with separately
	for(int side = 0; side < 2; side++)
	{		
		// filtering good tracker hits for reconstruction (no missing time stamps...)
		vector<TKtrhit*> hits_from_side = filter_side(tr_hits, side); 
		vector<TKtrhit*> hits = filter_usable(hits_from_side);

		if( hits.size() < 3 ) continue;
		
		// peaks_Theta, peak_R and peaks_value stores information about peak candidates
		// each iteration takes those candidates zooms around them looks for more candidates in that region
		// new set of candidates (more precise ones) comes out from each iteration
		vector<double> peaks_Theta = {M_PI/2.0};
		vector<double> peaks_R = {0.0};
		vector<double> peaks_value = {0.0};
		
		// size of region (delta_phi x delta_R) on witch each peak is calculated more precisely
		double delta_phi = M_PI;
		double delta_R = 5000.0;
		
		int segments = start_no_segments;
		for(int iter = 0; iter < iterations; iter++)
		{
			// the performance if better with only one segment in final iteration (algorithm looks only for one peak) 
			if( iter == (iterations - 1) )
			{
				segments = 1;
			}
			
			// stores candidates during iteration
			vector<double> peaks_Theta_temp;
			vector<double> peaks_R_temp;
			vector<double> peaks_value_temp;
			
			for(int i = 0; i < peaks_R.size(); i++)
			{		
				// arrays stores information about candidates from each segment
				double segment_peaks[segments][segments][3];
				
				// max_peak - stores value of the highest peak in all segments
				double max_peak = 0.0;
				
				// loop over segments
				for(int seg_r = 0; seg_r < segments; seg_r++)
				{
					for(int seg_theta = 0; seg_theta < segments; seg_theta++)
					{
						// boundaries of segments
						double r_min = peaks_R[i] - (delta_R/2.0) + (seg_r * delta_R/double(segments));
						double r_max = peaks_R[i] - (delta_R/2.0) + ((seg_r+1) * delta_R/double(segments));
						double theta_min = peaks_Theta[i] - (delta_phi/2.0) + (seg_theta * delta_phi/double(segments));
						double theta_max = peaks_Theta[i] - (delta_phi/2.0) + ((seg_theta+1) * delta_phi/double(segments));
						
						// sinograms are calculated for each bin of theta range - (offset = 1/2 of bin widht) 
						double offset = (delta_phi/double(segments))/(2.0*resolution);
						TH2F *sinograms = new TH2F("sinograms", "sinograms; theta; r", resolution, theta_min + offset, theta_max + offset, resolution, r_min, r_max);								
						
						double theta;
						double r;			
						for(int hit = 0; hit < hits.size(); hit++)
						{
							double sigma = hits[hit]->get_sigma_R();
							for(int k = 0; k <= resolution; k++)
							{
								theta = theta_min + ( k * (delta_phi/double(segments)) / double(resolution) );
								// r - legendre transform of a center of a circle
								r = hits[hit]->get_xy('x')*sin(theta) - hits[hit]->get_xy('y')*cos(theta);
								
								double weight;
								for(int half = 0; half < 2; half++)
								{	
									// mu - legendre transform of half circle (+r/-r)
									double mu = (r + (2.0*half - 1.0)*hits[hit]->get_r());	
									
									// gauss is calculated only for -3 to 3 sigma region to cut time							
									double r1 = mu - 3.0*sigma;
									double r2 = mu + 3.0*sigma;
									
									// bin numbers coresponding to r1 and r2 values
									int bin1 = (double(resolution) * (r1-r_min) / (r_max-r_min)) + 1;
									int bin2 = (double(resolution) * (r2-r_min) / (r_max-r_min)) + 1;
									
									// real values of r coresponding to each bin 
									double r_j1 = r_min + (r_max - r_min) * double(bin1) / double(resolution);
									double r_j2;
									for(int binj = bin1; binj < bin2 + 1; binj++)
									{
										r_j2 = r_j1 + (r_max - r_min) / double(resolution);
										
										// integrated probability
										//weight = (delta_phi / double(resolution * segments)) * ( erf( (r_j2 - mu)/sqrt(2.0*sigma) ) - erf( (r_j1 - mu)/sqrt(2.0*sigma) ) ) / (4.0*sigma); 
										
										// average probability density in a bin given by gauss distribution with mean in mu 
										// (uniformly distributed with respect to theta)
										weight = ( erf( (r_j2 - mu)/(sqrt(2.0)*sigma) ) - erf( (r_j1 - mu)/(sqrt(2.0)*sigma) ) ) / (2.0 * delta_R / double(resolution * segments));
										
										// result is 2D histogram of several sinusoid functions f(theta) in convolution with gauss with respect to R
										sinograms->Fill( theta, (r_j2 + r_j1)/2.0, weight );
										r_j1 = r_j2;			
									}								
								}			
							}	
						}										
						// sinogram for given segment is complete - next step is looking for peaks 
						
						double maximum = 0.0;
						double peak_Theta, peak_R;
						for(int j = 1; j <= resolution; j++)
						{
							for(int k = 1; k <= resolution; k++)
							{
								if(maximum < sinograms->GetBinContent(j,k))
								{
									// information (value, bin_theta, bin_R) about peak in given segment
									maximum = sinograms->GetBinContent(j,k);
									peak_Theta = j;
									peak_R = k;
								}
							}
						}
						
						if( save_sinograms == true ) 
						{
							TCanvas* c2 = new TCanvas("sinograms");
							sinograms->SetStats(0);
							sinograms->SetContour(100);
							sinograms->Draw("COLZ");
							c2->SaveAs(Form("Events_visu/sinograms-run-%d_event-%d_side-%d_iter-%d_R-%d_Th-%d.png", run_number, event_number, side, iter, seg_r, seg_theta));
							c2->Close();
						}
						delete sinograms;
						
						// converting bin numbers to real values
						peak_Theta = peaks_Theta[i] - (delta_phi/2.0) + (seg_theta * delta_phi/double(segments)) + ( peak_Theta * (delta_phi/double(segments)) / double(resolution) );
						peak_R = peaks_R[i] - (delta_R/2.0) + (seg_r * delta_R/double(segments)) + (delta_R/double(segments)) * (peak_R-0.5) / double(resolution);
						
						segment_peaks[seg_r][seg_theta][0] = peak_R;
						segment_peaks[seg_r][seg_theta][1] = peak_Theta;
						segment_peaks[seg_r][seg_theta][2] = maximum;
						if(maximum > max_peak)
						{
							// finding out the value of the highest peak in all segments
							max_peak = maximum;
						}
					}
				}
				
				for(int seg_r = 0; seg_r < segments; seg_r++)
				{
					for(int seg_theta = 0; seg_theta < segments; seg_theta++)
					{
						// filtering peak candidates obtained from region around given previous peak 
						if(segment_peaks[seg_r][seg_theta][2] >= treshold * max_peak /*&& segment_peaks[seg_r][seg_theta][2] > */)
						{
							peaks_R_temp.push_back(segment_peaks[seg_r][seg_theta][0]);
							peaks_Theta_temp.push_back(segment_peaks[seg_r][seg_theta][1]);
							peaks_value_temp.push_back(segment_peaks[seg_r][seg_theta][2]);
						}
					}
				}			
			}
			
			peaks_R.clear();
			peaks_Theta.clear();
			peaks_value.clear();
			
			// maximum value of all peaks 
			double maximum = 0.0;
			for(int i = 0; i < peaks_value_temp.size(); i++)
			{
				if( peaks_value_temp[i] > maximum )
				{
					maximum = peaks_value_temp[i];
				}
			}
			
			for(int i = 0; i < peaks_value_temp.size(); i++)
			{
				if( peaks_value_temp[i] > treshold * maximum )
				{	
					// if two peaks are so close that they would be duplicit in next iteration one is removed
					// otherwise we obtain multiple identical peaks
					bool is_new = true;
					for(int j = 0; j < peaks_R.size(); j++)
					{
						if(abs(peaks_R_temp[i] - peaks_R[j]) < delta_R/(2.0*segments))
						{
							if(abs(peaks_Theta_temp[i] - peaks_Theta[j]) < delta_phi/(2.0*segments))
							{
								is_new = false;
							}
						}
					}
					if(is_new)
					{
						peaks_R.push_back(peaks_R_temp[i]);
						peaks_Theta.push_back(peaks_Theta_temp[i]);
						peaks_value.push_back(peaks_value_temp[i]);						
					}
				}
			}
			
			// temp vectors are cleared for next iteration
			peaks_R_temp.clear();
			peaks_Theta_temp.clear();
			peaks_value_temp.clear();
			
			// regions of delta and theta are scaled down
			delta_phi = delta_phi/double(segments);
			delta_R = delta_R/double(segments);

			
			// last iteration - higher resolution to obtain more accuracy without getting more solutions
			// (not ideal solution) 
			if( iter == (iterations - 2) )
			{
				resolution = resolution * 4;
			}
		}
		// end of finding peaks
		// next step - creating tracks
	
		for(int j = 0; j < peaks_R.size(); j++)
		{
			
			double theta = peaks_Theta[j];
			double r     = peaks_R[j];
			
			TKtrack* track = new TKtrack(side, theta, r);	
			
			double a = track->get_a();
			double b = track->get_b();
			
			// associating hits to track
			double denominator = sqrt((a*a) + 1);	
			for(int i = 0; i < hits.size(); i++)
			{
				double distance_from_wire = abs(hits[i]->get_xy('y') - a*hits[i]->get_xy('x') - b) / denominator;
				double association_distance = 3.0*hits[i]->get_sigma_R();
				if( abs(distance_from_wire - hits[i]->get_r()) < association_distance )
				{
					track->add_associated_tr_hit(hits[i]);
					hits[i]->set_associated_track(track);
				}
			}			
			
			// calculating likelihood_R
			double joined_probability = 1.0;
			double probability;
			double norm = 1.0 / sqrt(2.0*M_PI);
			for(int i = 0; i < hits.size();i++)
			{	
				double sigma = hits[i]->get_sigma_R();
				double association_distance = 3.0*sigma;
				if(abs( r + hits[i]->get_xy('x')*cos(theta) - hits[i]->get_xy('y')*sin(theta) ) < association_distance )
				{
					probability = norm*exp(-pow(abs( r + hits[i]->get_xy('x')*cos(theta) - hits[i]->get_xy('y')*sin(theta) ) - hits[i]->get_r() , 2)/(2.0*sigma*sigma));
					joined_probability = joined_probability * probability / sigma; 
				}
			}		
			
			track->set_likelihood_R(joined_probability);
			track->reconstruct_vertical_least_square();
			tracks.push_back(track);				
		}
	}	
}


void TKEvent::reconstruct_single(bool save_sinograms)
{
	reconstruct_single_from_hits(tr_hits, save_sinograms);
}

void TKEvent::reconstruct_single_from_hits(std::vector<TKtrhit*>& _hits, bool save_sinograms)
{
	gROOT->SetBatch(kTRUE);

	// resolution of both dimensions of every 2D histogram	
	int resolution = 250; 
	
	// number of zooming iterations into histogram
	const int iterations = 2;
	
	// both sides are dealt with separately
	for(int side = 0; side < 2; side++)
	{	
		// filtering good tracker hits for reconstruction (no missing time stamps...)
		vector<TKtrhit*> hits_from_side = filter_side(_hits, side); 
		vector<TKtrhit*> hits = filter_usable(hits_from_side);
		if( hits.size() < 3 ) continue;

		// peaks_Theta, peak_R and peaks_value stores information about peak candidate
		// each iteration takes those candidate and zooms around them
		double peak_Theta = M_PI/2.0;
		double peak_R = 0.0;
		double peak_value = 0.0;
		
		// size of region (delta_phi x delta_R) on witch each peak is calculated more precisely
		double delta_phi = M_PI;
		double delta_R = 5000.0;
		
		for(int iter = 0; iter < iterations; iter++)
		{
			double r_min = peak_R - (delta_R/2.0);
			double r_max = peak_R + (delta_R/2.0);
			double theta_min = peak_Theta - (delta_phi/2.0);
			double theta_max = peak_Theta + (delta_phi/2.0);
			
			// sinograms are calculated for each bin of theta range - (offset = 1/2 of bin widht) 
			double offset = (delta_phi)/(2.0*resolution);
			TH2F *sinograms = new TH2F("sinograms", "sinograms; theta; r", resolution, theta_min + offset, theta_max + offset, resolution, r_min, r_max);								
			
			double theta;
			double r;			
			for(int hit = 0; hit < hits.size(); hit++)
			{
				double sigma = hits[hit]->get_sigma_R();
				for(int k = 0; k <= resolution; k++)
				{
					theta = theta_min + ( k * delta_phi / double(resolution) );
					// r - legendre transform of a center of a circle (Hough transform) - in old coordinates
					//r = ( -hits[hit]->get_xy('x')*cos(theta) ) - ( -hits[hit]->get_xy('y')*sin(theta) );
					// in new coordinates
					r = hits[hit]->get_xy('x')*sin(theta) - hits[hit]->get_xy('y')*cos(theta);
					
					double weight;
					for(int half = 0; half < 2; half++)
					{	
						// mu - legendre transform of half circle (+r/-r)
						double mu = (r + (2.0*half - 1.0)*hits[hit]->get_r());	
						
						// gauss is calculated only for -3 to 3 sigma region to cut time							
						double r1 = mu - 3.0*sigma;
						double r2 = mu + 3.0*sigma;
						
						// bin numbers coresponding to r1 and r2 values
						int bin1 = (double(resolution) * (r1-r_min) / (r_max-r_min)) + 1;
						int bin2 = (double(resolution) * (r2-r_min) / (r_max-r_min)) + 1;
						
						// real values of r coresponding to each bin 
						double r_j1 = r_min + (r_max - r_min) * double(bin1) / double(resolution);
						double r_j2;
						for(int binj = bin1; binj < bin2 + 1; binj++)
						{
							r_j2 = r_j1 + (r_max - r_min) / double(resolution);
							
							// average probability density in a bin given by gauss distribution with mean in mu 
							// (uniformly distributed with respect to theta)
							weight = ( erf( (r_j2 - mu)/(sqrt(2.0)*sigma) ) - erf( (r_j1 - mu)/(sqrt(2.0)*sigma) ) ) / (2.0 * delta_R / double(resolution));
							
							// result is 2D histogram of several sinusoid functions f(theta) in convolution with gauss with respect to R
							sinograms->Fill( theta, (r_j2 + r_j1)/2.0, weight );
							r_j1 = r_j2;			
						}								
					}			
				}	
			}										

			double maximum = 0.0;
			int peak_Theta_bin = 0;
			int peak_R_bin = 0;
			for(int j = 1; j <= resolution; j++)
			{
				for(int k = 1; k <= resolution; k++)
				{
					if(maximum < sinograms->GetBinContent(j,k))
					{
						// information (value, bin_theta, bin_R) about peak in given segment
						maximum = sinograms->GetBinContent(j,k);
						peak_Theta_bin = j;
						peak_R_bin = k;
					}
				}
			}
			
			if( save_sinograms == true ) 
			{
				TCanvas* c2 = new TCanvas("sinograms", "sinograms", 2000, 1600);
				sinograms->SetStats(0);
				sinograms->SetContour(100);
				sinograms->Draw("COLZ");
				c2->SaveAs(Form("Events_visu/sinograms-run-%d_event-%d_side-%d_iter-%d.png", run_number, event_number, side, iter));
				c2->Close();
			}
			delete sinograms;
			
			// converting bin numbers to real values
			peak_Theta = peak_Theta - (delta_phi/2.0) + (double(peak_Theta_bin) * (delta_phi) / double(resolution) );
			peak_R = peak_R - (delta_R/2.0) + (delta_R) * (double(peak_R_bin) - 0.5) / double(resolution);
			
			delta_phi = delta_phi*0.1;
			delta_R = delta_R*0.1;
		}
	
		// creating track
		double phi = peak_Theta;
		double r = peak_R;
		TKtrack* track = new TKtrack(side, phi, r);	
		
		double a = track->get_a();
		double b = track->get_b();
	
		// associating hits to track
		// calculating chi_squared_R, quality_R and likelihood_R
		double sum_R = 0.0;
		
		for(int i = 0; i < hits.size(); i++)
		{
			double distance = abs( hits[i]->get_xy('x')*sin(phi) - hits[i]->get_xy('y')*cos(phi) - r );
			distance = abs( distance - hits[i]->get_r() );
			double sigma = hits[i]->get_sigma_R();
			double association_distance = 3.0 * sigma;
			if( distance < association_distance )
			{
				track->add_associated_tr_hit(hits[i]);
				hits[i]->set_associated_track(track);
				
				sum_R = sum_R + distance/(sigma*sigma);
			}
		}	
		
		double no_hits = track->get_associated_tr_hits().size();
		double chi_squared = sum_R / no_hits;
		double quality = exp( -0.5*chi_squared );
		
		double norm = pow( 2.0*M_PI, -no_hits/2.0 ); 
		for(int i = 0; i < track->get_associated_tr_hits().size(); i++)
		{
			norm = norm / track->get_associated_tr_hits()[i]->get_sigma_R();
		}
		double likelihood = norm * exp( -0.5*chi_squared*no_hits );
		
		track->set_chi_squared_R( chi_squared );	
		track->set_quality_R( quality );	
		track->set_likelihood_R( likelihood );	

		track->reconstruct_vertical_least_square();
		//track->reconstruct_vertical_MLM();
		tracks.push_back(track);				
	}	
}

void TKEvent::draw_sinusoids()
{
	gROOT->SetBatch(kTRUE);
	
	for(int side = 0; side < 2; side++)
	{
		vector<TKtrhit*> hits_from_side = filter_side(tr_hits, side); 
		vector<TKtrhit*> hits = filter_usable(hits_from_side);
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

	for(int side = 0; side < 2; side++)
	{
		vector<TKtrhit*> hits_from_side = filter_side(tr_hits, side); 
		vector<TKtrhit*> hits = filter_usable(hits_from_side);

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
				
			for(int hit = 0; hit < hits.size(); hit++)
			{
				double R0 = hits[hit]->get_r();
				double sigma_R = hits[hit]->get_sigma_R();
				double norm_const = 1.0 / (sqrt(2.0*M_PI)*sigma_R); 
				
				double weight;
				double r, theta;
				for(int j = 0; j <= resolution; j++)
				{
					theta = phi1 + (j * (phi2 - phi1) / resolution);
					double rho = hits[hit]->get_xy('x')*sin(theta) - hits[hit]->get_xy('y')*cos(theta);
					
					for(int k = 0; k <= resolution; k++)
					{
						r = R1 + (k * (R2 - R1) / resolution);
						
						weight = exp( -pow( abs(r-rho) - hits[hit]->get_r() , 2 ) / (2.0*sigma_R*sigma_R) );
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
			double norm_const = pow(2.0*M_PI, -0.5*hits.size());
			for(int i = 0; i < hits.size(); i++)
			{	
				double sigma_R = hits[i]->get_sigma_R();
				norm_const = norm_const / sigma_R;
				probability = exp(-pow(abs( peak_R + hits[i]->get_xy('x')*cos(peak_Theta) - hits[i]->get_xy('y')*sin(peak_Theta) ) - hits[i]->get_r() , 2)/(2.0*sigma_R*sigma_R));
				joined_probability = joined_probability * probability; 
			}		
			
			joined_probability = norm_const * joined_probability;
			
			cout << "likelihood R: " << joined_probability << endl;
		}
	}
}

void TKEvent::hough_transform(std::vector<TKtrhit*>& hits, double phi_min, double phi_max, double R_min, double R_max, int ID)
{
	gROOT->SetBatch(kTRUE);

	const int resolution = 1000;
	double r, theta;

	double S_r = 0.0;
	double S_rX = 0.0;
	double S_rY = 0.0;
	double const_R;
	for(int i = 0; i < hits.size();i++)
	{
		const_R = pow(hits[i]->get_sigma_R(), -2.0);
		S_r = S_r + const_R;
		S_rX = S_rX + hits[i]->get_xy('x')*const_R;
		S_rY = S_rY + hits[i]->get_xy('y')*const_R;
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

	for(int side = 0; side < 2; side++)
	{
		vector<TKtrhit*> hits_from_side = filter_side(tr_hits, side); 
		vector<TKtrhit*> hits = filter_usable(hits_from_side);
		if( hits.size() < 3 ) continue;

		double S_r = 0.0;
		double S_rX = 0.0;
		double S_rY = 0.0;
		double const_R;
		for(int i = 0; i < hits.size();i++)
		{
			const_R = hits[i]->get_sigma_R();
			S_r = S_r + const_R;
			S_rX = S_rX + hits[i]->get_xy('x')*const_R;
			S_rY = S_rY + hits[i]->get_xy('y')*const_R;
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
						double sigma_R = hits[hit]->get_sigma_R();
						double norm_const = 1.0 / (sqrt(2.0*M_PI)*sigma_R); 
						
						temp = exp( -pow( abs(r-rho) - hits[hit]->get_r() , 2 ) / (2.0*sigma_R*sigma_R) );
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
			double norm_const = pow(2.0*M_PI, -0.5*hits.size());
			for(int i = 0; i < hits.size(); i++)
			{	
				double sigma_R = hits[i]->get_sigma_R();
				norm_const = norm_const / sigma_R;
				probability = exp(-pow(abs( peak_R + hits[i]->get_xy('x')*cos(peak_Theta) - hits[i]->get_xy('y')*sin(peak_Theta) ) - hits[i]->get_r() , 2)/(2.0*sigma_R*sigma_R));
				joined_probability = joined_probability * probability; 
			}	
			
			joined_probability = norm_const * joined_probability;
			
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

void TKEvent::reconstruct(bool save_sinograms)
{
	//const bool debug_mode = false;
	const double chi_square_treshold = 5;
	for(int side = 0; side < 2; side++)
	{
		//if(debug_mode) cout << "side: " << side << endl;
		vector<TKtrhit*> hits = filter_side(tr_hits, side); 
		hits = filter_usable( hits );
		hits = filter_distant( hits );
		int no_hits_before = 1e10;
		while( hits.size() > 2 && no_hits_before > hits.size() )
		{
			//if(debug_mode) cout << "iteration start: " << endl;
			no_hits_before = hits.size();
			bool failed = true;
			TKcluster* cluster = find_cluster(hits);			
			if( cluster != nullptr )
			{
				//if(debug_mode) cout << "find_cluster succes" << endl; 
				cluster->reconstruct_MLM( save_sinograms, run_number, event_number );						
				if( cluster->get_track()->get_chi_squared_R() < chi_square_treshold && cluster->get_track()->get_associated_tr_hits().size() > 3 )
				{
					//if(debug_mode) cout << "chi_squared good" << endl;
					failed = false;
					cluster->detect_ambiguity_type();
					cluster->reconstruct_ambiguity(); 
					cluster->get_track()->calculate_tr_hit_points();
					if(cluster->get_track()->get_mirror_image() != nullptr)
					{	
						cluster->get_track()->get_mirror_image()->calculate_tr_hit_points();
					}
					clusters.push_back( cluster );
				}
				else
				{
					//if(debug_mode) cout << "chi_squared bad" << endl; 
					delete cluster;
				}
			}
			if( failed )
			{
				TKcluster* cluster = find_cluster_legendre(hits, save_sinograms);
				if( cluster != nullptr )
				{
					//if(debug_mode) cout << "find_cluster_legendre succes" << endl; 
					cluster->reconstruct_MLM( save_sinograms, run_number, event_number );
					if( cluster->get_track()->get_chi_squared_R() < chi_square_treshold && cluster->get_track()->get_associated_tr_hits().size() > 3 )
					{
						//if(debug_mode) cout << "chi_squared good" << endl; 
						cluster->detect_ambiguity_type();
						cluster->reconstruct_ambiguity();
						cluster->get_track()->calculate_tr_hit_points();
						if(cluster->get_track()->get_mirror_image() != nullptr)
						{	
							cluster->get_track()->get_mirror_image()->calculate_tr_hit_points();
						}
						clusters.push_back( cluster );
					}	
					else
					{
						//if(debug_mode) cout << "chi_squared bad" << endl; 
						delete cluster;
					}
				}
			}
			hits = filter_unassociated( hits );
			hits = filter_distant( hits );
			//hits = filter_unclustered( hits );
		}
	}
	this->build_trajectories();
	this->extrapolate_trajectories();
}

void TKEvent::reconstruct_simple(bool save_sinograms)
{
	//const bool debug_mode = true;
	const double chi_square_treshold = 5;
	for(int side = 0; side < 2; side++)
	{
		//if(debug_mode) cout << "side: " << side << endl;
		vector<TKtrhit*> hits = filter_side(tr_hits, side); 
		hits = filter_usable( hits );
		hits = filter_distant( hits );
		if( hits.size() < 3 ) continue; 
		bool failed = true;
		TKcluster* cluster = find_cluster(hits);			
		if( cluster != nullptr )
		{
			//if(debug_mode) cout << "find_cluster succes" << endl; 
			cluster->reconstruct_MLM( save_sinograms, run_number, event_number );						
			if( cluster->get_track()->get_chi_squared_R() < chi_square_treshold && cluster->get_track()->get_associated_tr_hits().size() > 2 )
			{
				failed = false;
				//if(debug_mode) cout << "chi_squared good" << endl;
				cluster->detect_ambiguity_type();
				cluster->reconstruct_ambiguity(); 
				cluster->get_track()->calculate_tr_hit_points();
				if(cluster->get_track()->get_mirror_image() != nullptr)
				{	
					cluster->get_track()->get_mirror_image()->calculate_tr_hit_points();
				}
				clusters.push_back( cluster );
			}
			else
			{
				//if(debug_mode) cout << "chi_squared bad" << endl; 
				delete cluster;
			}
		}
		if( failed )
		{
			TKcluster* cluster = find_cluster_legendre(hits, save_sinograms);
			if( cluster != nullptr )
			{
				//if(debug_mode) cout << "find_cluster_legendre succes" << endl; 
				cluster->reconstruct_MLM( save_sinograms, run_number, event_number );
				if( cluster->get_track()->get_chi_squared_R() < chi_square_treshold && cluster->get_track()->get_associated_tr_hits().size() > 2 )
				{
					//if(debug_mode) cout << "chi_squared good" << endl; 
					cluster->detect_ambiguity_type();
					cluster->reconstruct_ambiguity();
					cluster->get_track()->calculate_tr_hit_points();
					if(cluster->get_track()->get_mirror_image() != nullptr)
					{	
						cluster->get_track()->get_mirror_image()->calculate_tr_hit_points();
					}
					clusters.push_back( cluster );
				}	
				else
				{
					//if(debug_mode) cout << "chi_squared bad" << endl; 
					delete cluster;
				}
			}
		}
	}
	this->build_trajectories();
	this->extrapolate_trajectories();
}

void TKEvent::reconstruct_ML(bool save_sinograms)
{
	for(int side = 0; side < 2; side++)
	{
		vector<TKtrhit*> hits = filter_side(tr_hits, side); 
		hits = filter_usable( hits );
		hits = filter_distant( hits );
		if( hits.size() < 3 ) continue;
		
		TKcluster* cluster = find_cluster(hits);
		if( cluster != nullptr )
		{
			clusters.push_back( cluster );
			cluster->reconstruct_MLM( save_sinograms, run_number, event_number );
			cluster->detect_ambiguity_type();
			cluster->reconstruct_ambiguity();
		}
	}
}


void TKEvent::reconstruct_ML_3D(bool save_sinograms)
{
	for(int side = 0; side < 2; side++)
	{
		vector<TKtrhit*> hits = filter_side(tr_hits, side); 
		hits = filter_usable( hits );
		hits = filter_distant( hits );
		if( hits.size() < 3 ) continue;
		
		TKcluster* cluster = find_cluster(hits);
		if( cluster != nullptr )
		{
			clusters.push_back( cluster );
			cluster->reconstruct_MLM_3D( save_sinograms, run_number, event_number );
			cluster->detect_ambiguity_type();
			cluster->reconstruct_ambiguity();
		}
	}
}

TKcluster* TKEvent::find_cluster(std::vector<TKtrhit*>& hits)
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
	cluster_hits = filter_distant(cluster_hits);//TODO is this good? 
	if(cluster_hits.size() < 3)
	{
		return nullptr;
	}
	TKcluster* cluster = new TKcluster(cluster_hits, phi_min, phi_max);
	return cluster;
}

TKcluster* TKEvent::find_cluster_legendre(std::vector<TKtrhit*>& hits, bool save_sinograms)
{
	gROOT->SetBatch(kTRUE);

	if( hits.size() < 3 ) 
	{
		cout << "cluster not found - too few usable hits" << endl;
		return nullptr;
	}

	const double distance_limit = 6.0;
	const double limit_angle = 5.0*M_PI/180.0;

	// resolution of both dimensions of every 2D histogram	
	int resolution = 250; 
	
	// number of zooming iterations into histogram
	const int iterations = 2;

	// peaks_phi, peak_R and peaks_value stores information about peak candidate
	// each iteration takes those candidate and zooms around them
	double peak_phi = M_PI/2.0;
	double peak_R = 0.0;
	double peak_value = 0.0;
	
	// size of region (delta_phi x delta_R) on witch each peak is calculated more precisely
	double delta_phi = M_PI;
	double delta_R = 5000.0;
	
	for(int iter = 0; iter < iterations; iter++)
	{
		double r_min = peak_R - (delta_R/2.0);
		double r_max = peak_R + (delta_R/2.0);
		double phi_min = peak_phi - (delta_phi/2.0);
		double phi_max = peak_phi + (delta_phi/2.0);
		
		// sinograms are calculated for each bin of phi range - (offset = 1/2 of bin widht) 
		double offset = (delta_phi)/(2.0*resolution);
		TH2F *sinograms = new TH2F("sinograms", "sinograms; phi; r", resolution, phi_min + offset, phi_max + offset, resolution, r_min, r_max);								
		
		double phi;
		double r;			
		for(int hit = 0; hit < hits.size(); hit++)
		{
			double sigma = hits[hit]->get_sigma_R();
			for(int k = 0; k <= resolution; k++)
			{
				phi = phi_min + ( k * delta_phi / double(resolution) );
				// r - legendre transform of a center of a circle (Hough transform)
				r = hits[hit]->get_xy('x')*sin(phi) - hits[hit]->get_xy('y')*cos(phi);
				
				double weight;
				for(int half = 0; half < 2; half++)
				{	
					// mu - legendre transform of half circle (+r/-r)
					double mu = (r + (2.0*half - 1.0)*hits[hit]->get_r());	
					
					// gauss is calculated only for -3 to 3 sigma region to cut time							
					double r1 = mu - 3.0*sigma;
					double r2 = mu + 3.0*sigma;
					
					// bin numbers coresponding to r1 and r2 values
					int bin1 = (double(resolution) * (r1-r_min) / (r_max-r_min)) + 1;
					int bin2 = (double(resolution) * (r2-r_min) / (r_max-r_min)) + 1;
					
					// real values of r coresponding to each bin 
					double r_j1 = r_min + (r_max - r_min) * double(bin1) / double(resolution);
					double r_j2;
					for(int binj = bin1; binj < bin2 + 1; binj++)
					{
						r_j2 = r_j1 + (r_max - r_min) / double(resolution);
						
						// average probability density in a bin given by gauss distribution with mean in mu 
						// (uniformly distributed with respect to phi)
						weight = ( erf( (r_j2 - mu)/(sqrt(2.0)*sigma) ) - erf( (r_j1 - mu)/(sqrt(2.0)*sigma) ) ) / (2.0 * delta_R / double(resolution));
						
						// result is 2D histogram of several sinusoid functions f(phi) in convolution with gauss with respect to R
						sinograms->Fill( phi, (r_j2 + r_j1)/2.0, weight );
						r_j1 = r_j2;			
					}								
				}			
			}	
		}										

		// Get bin number of maximum value
		int maxBin = sinograms->GetMaximumBin();

		// Get X and Y values corresponding to the maximum bin
		int bin_phi, bin_R, bin_Z;
		sinograms->GetBinXYZ(maxBin, bin_phi, bin_R, bin_Z);
		peak_phi = sinograms->GetXaxis()->GetBinCenter(bin_phi);
		peak_R = sinograms->GetYaxis()->GetBinCenter(bin_R);

		delta_phi = delta_phi*0.1;
		delta_R = delta_R*0.1;


		if( save_sinograms ) 
		{
			TCanvas* c2 = new TCanvas("sinograms", "sinograms", 2000, 1600);
			sinograms->SetStats(0);
			sinograms->SetContour(100);
			sinograms->Draw("COLZ");
			c2->SaveAs(Form("Events_visu/clustering-run-%d_event-%d_iter-%d.png", run_number, event_number, iter));
			c2->Close();
			delete c2;
		}
		delete sinograms;
	}
	
	vector<TKtrhit*> cluster_candidate = filter_close_hits(hits, peak_phi, peak_R, distance_limit);
	
	cluster_candidate = filter_distant(cluster_candidate); //TODO is this good?
	if(cluster_candidate.size() < 3)
	{
		return nullptr;
	}
	else
	{
		TKcluster* cluster = new TKcluster(cluster_candidate, peak_phi - limit_angle, peak_phi + limit_angle);
		return cluster;
	}
}

void TKEvent::build_trajectories()
{
	//if(debug_mode) cout << "building trajectory" << endl;
	vector<TKtrack*> all_tracks = this->get_tracks();
	
	
	
	for(int side = 0; side < 2; side++)
	{
		//if(debug_mode) cout << "building side: " << side << endl;
		vector<TKtrack*> all_tracks_from_side;
		for(int i = 0; i < all_tracks.size(); i++)
		{
			if(all_tracks[i]->get_side() == side)
			{
				if(all_tracks[i]->get_associated_tr_hit_points().size() > 1)
				all_tracks_from_side.push_back(all_tracks[i]);
			}
		}		
		const int no_tracks = all_tracks_from_side.size();
		//if(debug_mode) cout << no_tracks << " tracks available" << endl;
		bool connections[no_tracks][no_tracks];
		for(int i = 0; i < no_tracks; i++)
		{
			for(int j = 0; j < no_tracks; j++)
			{
				connections[i][j] = 0;
			}
		}
		for(int i = 0; i < no_tracks; i++)
		{
			TKtrack* track1 = all_tracks_from_side[i];
			for(int j = i+1; j < no_tracks; j++)
			{
				TKtrack* track2 = all_tracks_from_side[j];
				if(track1->get_mirror_image() == track2) continue;
				
				double a1 = track1->get_a();
				double a2 = track2->get_a();
				double b1 = track1->get_b();
				double b2 = track2->get_b();
				double x = (b2-b1)/(a1-a2);
				double y = a1*x + b1;
				if(y < -2486.0 || y > 2486.0) continue;
				if(x < -425.0 || x > 425.0) continue;
				if(side == 0 && x > -29.0) continue;
				if(side == 1 && x < 29.0) continue;				
				double z1 = track1->get_c()*x + track1->get_d();
				double z2 = track2->get_c()*x + track2->get_d();
				if(z1 == 0 || z2 == 0) continue;
				if(abs(z2 - z1) > 40.0) continue;
				
				connections[i][j] = 1;
				connections[j][i] = 1;			
			}
		}
		int connection_counter[no_tracks];	
		for(int i = 0; i < no_tracks; i++)
		{
			connection_counter[i] = 0;
			for(int j = 0; j < no_tracks; j++)
			{
				if(connections[j][i] == 1)
				{
					connection_counter[i]++;
				}
			}
		}
		bool trajectorized[no_tracks];
		//if(debug_mode) cout << "number of connections for each segment:" << endl;
		for(int i = 0; i < no_tracks; i++)
		{
			trajectorized[i] = false;
			if(debug_mode) cout << connection_counter[i] << endl; 
		}
		/*
		if(debug_mode)
		{
			cout << "connection matrix:" << endl;
			for(int i = 0; i < no_tracks; i++)
			{
				for(int j = 0; j < no_tracks; j++)
				{
					cout << connections[i][j] << " ";	
				}
				cout << endl;
			}
			cout << endl;
		}*/
		
		for(int i = 0; i < no_tracks; i++)
		{
			//if(debug_mode) cout << "connecting track number: " << i << endl; 
			if(connection_counter[i] == 0)
			{
				//if(all_tracks_from_side[i]->get_associated_tr_hits().size() > 1) // in ideal case should be unnecessary
				trajectories.push_back(new TKtrajectory(all_tracks_from_side[i]));
				trajectorized[i] = true;
				//if(debug_mode) cout << i << " trajectorized" << endl;
			}
			else if(connection_counter[i] == 1 && trajectorized[i] == false)
			{
				vector<TKtrack*> composite_track;
				composite_track.push_back(all_tracks_from_side[i]);
				trajectorized[i] = true;
				//if(debug_mode) cout << i << " trajectorized" << endl;
				int next_index;
				bool is_valid = false;
				for(int j = 0; j < no_tracks; j++)
				{
					if(connections[i][j] == 1 && trajectorized[j] == false)
					{
						next_index = j;
						is_valid = true;
						composite_track.push_back(all_tracks_from_side[next_index]);
						trajectorized[next_index] = true;
						//if(debug_mode) cout << next_index << " trajectorized" << endl;
						break;
					}
				}
				
				while(connection_counter[next_index] >= 2 && is_valid)
				{
					is_valid = false;
					for(int j = 0; j < no_tracks; j++)
					{
						if(connections[next_index][j] == 1 && trajectorized[j] == false)
						{
							next_index = j;
							is_valid = true;
							composite_track.push_back(all_tracks_from_side[next_index]);
							trajectorized[next_index] = true;
							//if(debug_mode) cout << next_index << " trajectorized" << endl;
							break;
						}
					}
				}
				if(composite_track.size() > 1)
				{
					trajectories.push_back(new TKtrajectory(composite_track));
				}
				else if(composite_track.size() == 1)
				{
					//if(composite_track[0]->get_associated_tr_hits().size() > 1) // in ideal case should be unnecessary
					trajectories.push_back(new TKtrajectory(composite_track[0]));
				}
			}
		}
	}
	/*
	if(debug_mode)
	{
		cout << "trajectory:" << endl;
		for(int i = 0; i < trajectories.size(); i++)
		{
			trajectories[i]->print();
		}
	}*/
	
	return;
}

void TKEvent::extrapolate_trajectories()
{
	for(auto& trajectory : trajectories)
	{
		trajectory->extrapolate();
	}
}

void TKEvent::calculate_tr_hit_points()
{
	vector<TKtrack*> all_tracks = this->get_tracks();
	for(auto& track : all_tracks)
	{
		track->calculate_tr_hit_points();
	}
}

