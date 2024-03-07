// TK headers
#include "TKEvent.h"

// dimensions in mm
// origin in the center of detector
const double Bi_source_x = 1.4;
const double Bi_source_y = 14.2;
const double Bi_source_z = 21.1;
const double Bi_source_dist_y = 835.0; // possibly wrong or not accurate (originally 850)
const double Bi_source_dist_z = 475.0; // possibly wrong or not accurate (originally 425)

// dimensions in mm
// origin in the center of detector
static double foil_spacex = 58.0; // possibly inaccurate
const double tc_radius = 22.0;
const double tc_sizez = 3030.0;

// OM dimensions in mm
// origin in the center of detector
// warning: dimensions are approximately taken from Falaise and might not be exactly right
const double mw_sizex = 194.0;
const double mw_sizey = 256.0;
const double mw_sizez = 256.0;

const double gv_sizex = 308.0;
const double gv_sizey = 310.0;
const double gv_sizez = 150.0;

const double xw_sizex = 200.0;
const double xw_sizey = 150.0;
const double xw_sizez = 208.5;

using namespace std;
ClassImp(TKEvent);

TKEvent::TKEvent()
{
	OM_hits   = std::vector<TKOMhit*>();
	tr_hits   = std::vector<TKtrhit*>();
	tracks    = std::vector<TKtrack*>();
	clusters  = std::vector<TKcluster*>();
}

TKEvent::TKEvent(int _run_number ,int _event_number)
{
	event_number = _event_number;
	run_number = _run_number;

	OM_hits   = std::vector<TKOMhit*>();
	tr_hits   = std::vector<TKtrhit*>();
	tracks    = std::vector<TKtrack*>();
	clusters  = std::vector<TKcluster*>();
}

TKEvent::~TKEvent()
{
	for(int i = 0;i < OM_hits.size();i++)
	{
		delete OM_hits[i];
	}
	for(int i = 0;i < tr_hits.size();i++)
	{
		delete tr_hits[i];
	}
	for(int i = 0;i < tracks.size();i++)
	{
		delete tracks[i];
	}
	for(int i = 0;i < clusters.size();i++)
	{
		delete clusters[i];
	}
	
	OM_hits.clear();
	tr_hits.clear();
	tracks.clear();
	clusters.clear();
}

std::vector<TKOMhit*> TKEvent::get_OM_hits()
{
	return OM_hits;
}

std::vector<TKtrhit*> TKEvent::get_tr_hits()
{
	return tr_hits;
}

std::vector<TKtrack*> TKEvent::get_tracks()
{
	vector<TKtrack*> all_tracks;
	for(int i = 0; i < tracks.size(); i++)
	{
		all_tracks.push_back( tracks[i] );
	}
	for(int i = 0; i < clusters.size(); i++)
	{
		if( clusters[i]->get_track() != nullptr )
		{
			all_tracks.push_back( clusters[i]->get_track() );
			if( clusters[i]->get_track()->get_mirror_image() != nullptr )
			{
				all_tracks.push_back( clusters[i]->get_track()->get_mirror_image() );
			}
		}
	}
	return all_tracks;
}

std::vector<TKcluster*> TKEvent::get_clusters()
{
	return clusters;
}

TKOMhit* TKEvent::get_OM_hit(int _i)
{
	return OM_hits[_i];
}

TKtrhit* TKEvent::get_tr_hit(int _i)
{
	return tr_hits[_i];
}

TKtrack* TKEvent::get_track(int _i)
{
	vector<TKtrack*> all_tracks;
	for(int i = 0; i < tracks.size(); i++)
	{
		all_tracks.push_back( tracks[i] );
	}
	for(int i = 0; i < clusters.size(); i++)
	{
		if( clusters[i]->get_track() != nullptr )
		{
			all_tracks.push_back( clusters[i]->get_track() );
			if( clusters[i]->get_track()->get_mirror_image() != nullptr )
			{
				all_tracks.push_back( clusters[i]->get_track()->get_mirror_image() );
			}
		}
	}
	return all_tracks[_i];
}

TKcluster* TKEvent::get_cluster(int _i)
{
	return clusters[_i];
}

int TKEvent::get_run_number()
{
	return run_number;
}

int TKEvent::get_event_number()
{	
	return event_number;
}

int TKEvent::get_no_tracks()
{
	vector<TKtrack*> all_tracks;
	for(int i = 0; i < tracks.size(); i++)
	{
		all_tracks.push_back( tracks[i] );
	}
	for(int i = 0; i < clusters.size(); i++)
	{
		if( clusters[i]->get_track() != nullptr )
		{
			all_tracks.push_back( clusters[i]->get_track() );
			if( clusters[i]->get_track()->get_mirror_image() != nullptr )
			{
				all_tracks.push_back( clusters[i]->get_track()->get_mirror_image() );
			}
		}
	}
	return all_tracks.size();
}
		
void TKEvent::print()
{
	std::cout << std::endl;
	std::cout << "RUN " << run_number << " | EVENT " << event_number << std::endl << std::endl;
	std::cout << "Collection of the OM hits: " << std::endl;

	for (int i = 0; i < OM_hits.size(); i++)
	{
		OM_hits[i]->print();
	}

	std::cout << std::endl;
	std::cout << "Collection of the tracker hits: " << std::endl;
	
	for (int i = 0; i < tr_hits.size(); i++)
	{
		tr_hits[i]->print();
	}

	std::cout << std::endl;
	std::cout << "Collection of the tracks: " << std::endl;
	
	for (int i = 0; i < this->get_tracks().size(); i++)
	{
		this->get_tracks()[i]->print();
	}
	std::cout << std::endl;
}

void TKEvent::print_tracks()
{
	std::cout << std::endl;
	std::cout << "RUN " << run_number << " | EVENT " << event_number << std::endl << std::endl;
	for (int i = 0; i < this->get_tracks().size(); i++)
	{
		this->get_tracks()[i]->print();
	}
	std::cout << std::endl;
	
}

void TKEvent::add_OM_hit(int _OM_num, bool _is_HT, int64_t _OM_TDC, int16_t _OM_pcell)
{
	OM_hits.push_back(new TKOMhit(_OM_num, _is_HT, _OM_TDC, _OM_pcell));
}
		
void TKEvent::add_OM_hit(int _SWCR[4], bool _is_HT, int64_t _OM_TDC, int16_t _OM_pcell)
{
	OM_hits.push_back(new TKOMhit(_SWCR, _is_HT, _OM_TDC, _OM_pcell));
}
	
void TKEvent::add_tracker_hit(int _cell_num, int64_t _tsp[7])
{
	tr_hits.push_back(new TKtrhit(_cell_num, _tsp));
}
		
void TKEvent::add_tracker_hit(int _SRL[3], int64_t _tsp[7])
{
	tr_hits.push_back(new TKtrhit(_SRL, _tsp));
}

std::vector<TKtrhit*> TKEvent::filter_side(std::vector<TKtrhit*> _hits, int side)
{
	vector<TKtrhit*> hits;	
	for(int i = 0; i < _hits.size(); i++)
	{
		if( side == _hits[i]->get_SRL('s'))
		{
			hits.push_back( _hits[i] );
		}
	}
	return hits;
}

std::vector<TKtrhit*> TKEvent::filter_usable(std::vector<TKtrhit*> _hits)
{
	vector<TKtrhit*> hits;	
	for(int i = 0; i < _hits.size(); i++)
	{
		
		// not using broken or too big (incorrectly associated) tracker hits
		if( _hits[i]->get_r() != -1.0 && _hits[i]->get_r() < 35.0 && _hits[i]->get_r() > 2.0 )
		{
			hits.push_back( _hits[i] );
		}
	}
	return hits;
}

std::vector<TKtrhit*> TKEvent::filter_close_hits(std::vector<TKtrhit*> _hits, double phi, double r, double distance_limit)
{
	vector<TKtrhit*> hits;	
	for(int i = 0; i < _hits.size(); i++)
	{
		double R = _hits[i]->get_r();
		double x = _hits[i]->get_xy('x');
		double y = _hits[i]->get_xy('y');
		double distance = abs(r - x*sin(phi) + y*cos(phi)) - R;
		if( abs(distance) <= distance_limit )
		{
			hits.push_back( _hits[i] );
		}
	}
	return hits;
}

void TKEvent::set_r(std::string drift_model, std::string association_mode)
{	
	for(int tr_hit = 0; tr_hit < tr_hits.size(); tr_hit++)
	{
		double r;
		if(tr_hits[tr_hit]->get_tsp('0') != -1)
		{
			double min_time = 1e10;
			// associates tracker hits to OM with minimal time difference
			if( association_mode == "time" ) 
			{
				//double calo_hit = (calo_tdc * 6.25) - 400.0 + (400.0 * peak_cell / 1024.0);
				for(int om_hit = 0; om_hit < OM_hits.size(); om_hit++)
				{
					int64_t TDC_diff = 2*tr_hits[tr_hit]->get_tsp('0') - OM_hits[om_hit]->get_OM_TDC() + 44;
					if(       TDC_diff  < 800 && 
						  TDC_diff  > 0   &&
					   double(TDC_diff) < min_time)
					{
						min_time = 6.25 * TDC_diff;
						tr_hits[tr_hit]->set_associated_OMhit(OM_hits[om_hit]);
					}
				}
			}

			// associates tracker hits to OM with minimal distance
			else if( association_mode == "distance" )
			{
				double min_distance = 1e10;
				for(int om_hit = 0; om_hit < OM_hits.size(); om_hit++)
				{
					double delta_x = tr_hits[tr_hit]->get_xy('x') - OM_hits[om_hit]->get_xyz('x');
					double delta_y = tr_hits[tr_hit]->get_xy('y') - OM_hits[om_hit]->get_xyz('y');
					double delta_z = tr_hits[tr_hit]->get_h()     - OM_hits[om_hit]->get_xyz('z');
					double distance = sqrt( delta_x*delta_x + delta_y*delta_y + delta_z*delta_z );
					if( distance < min_distance ) 
					{
						int64_t TDC_diff = 2*tr_hits[tr_hit]->get_tsp('0') - OM_hits[om_hit]->get_OM_TDC() + 44;
						if( TDC_diff <= 0 ) 
						{
							continue;
						}
						min_time = 6.25 * TDC_diff;
						tr_hits[tr_hit]->set_associated_OMhit(OM_hits[om_hit]);
						min_distance = distance;
					}
					
				}
			}
			else cout << "invalid association model: choose \"distance\" or \"time\"." << endl;

			if( drift_model == "Manchester" )
			{
				const double A1 = 0.570947153108633;
				const double B1 = 0.580148313540993;
				const double C1 = 1.6567483468611;
				const double A2 = 1.86938462695651;
				const double B2 = 0.949912427483918;

				const double t_usec = min_time / 1000.0;
				const double ut = 10. * t_usec;
		
				r = A1 * ut / (std::pow(ut, B1) + C1);
				if (r > A2 * ut / (std::pow(ut, B2))) 
				{
					r = A2 * ut / (std::pow(ut, B2));
				}	

				r *= 10.0;	
			}
			// modification of Betsy's model
			else if( drift_model == "Betsy" )
			{
				const double a1 = 0.828;
				const double b1 = -0.907;
				const double a2 = 0.402;
				const double b2 = -1.955;
		
				min_time = min_time / 1000.0;
				if(min_time < 3.0845 && min_time > 0.0)
				{
					r = pow(min_time/a1, 1.0/(1.0-b1)) * 10.0;
				}
				else
				{
					r = pow(min_time/a2, 1.0/(1.0-b2)) * 10.0;
				}
			}
			else cout << "invalid drift time model: choose \"Betsy\" or \"Manchester\"." << endl;
		}
		else
		{
			r = -1.0;
		}
		
		tr_hits[tr_hit]->set_r(r);
		tr_hits[tr_hit]->set_sigma_R();
	}
}

void TKEvent::set_sigma_R()
{
	for(int tr_hit = 0; tr_hit < tr_hits.size(); tr_hit++)
	{
		tr_hits[tr_hit]->set_sigma_R();
	}
}

void TKEvent::set_h()
{
	for(int tr_hit = 0; tr_hit < tr_hits.size(); tr_hit++)
	{
		tr_hits[tr_hit]->set_h();
	}
}

void TKEvent::set_sigma_Z()
{
	for(int tr_hit = 0; tr_hit < tr_hits.size(); tr_hit++)
	{
		tr_hits[tr_hit]->set_sigma_Z();
	}
}
		
void TKEvent::make_top_projection(int option)
{
	gROOT->SetBatch(true);
	TCanvas *canvas = new TCanvas("canvas","", 5800, 1600);	
	canvas->Range(-2900.0, -700.0, 2900.0, 900.0);
	TLatex* title = new TLatex(-2600.0, 700.0, Form("Run %d | Event %d", run_number, event_number));
	title->SetTextSize(0.07);
	title->Draw();
	

	// Drawing mainwall and mainwall calo hits
	for(int om_side = 0; om_side < 2; om_side++) 
	{
		for(int om_column = 0; om_column < 20; om_column++) 
		{
			int swcr[4] = {om_side, -1, om_column, 0};
			TKOMhit* ohit = new TKOMhit(swcr); 
			
			TBox *calo = new TBox(ohit->get_xyz('y') - mw_sizey/2.0, 
					    -(ohit->get_xyz('x') + mw_sizex/2.0), 
					      ohit->get_xyz('y') + mw_sizey/2.0, 
					    -(ohit->get_xyz('x') - mw_sizex/2.0));
			
			TLatex *tex = new TLatex(ohit->get_xyz('y') - mw_sizey/4.0,-ohit->get_xyz('x')
				,Form("%d.%d.%d.-", ohit->get_SWCR('s'), ohit->get_SWCR('w'), ohit->get_SWCR('c')));
			tex->SetTextSize(0.035);
			
			bool is_hit = false;
			for(int hit = 0; hit < OM_hits.size(); hit++)
			{
// MIRO: Zaviesť TKOMhit::operator==(...) a porovnávať totožnosť OM pomocou neho!
				if(om_side   == OM_hits[hit]->get_SWCR('s') && 
				        -1   == OM_hits[hit]->get_SWCR('w') && 
				   om_column == OM_hits[hit]->get_SWCR('c') ) 
				{
				   is_hit = true;
				}
			}
			if(is_hit)
			{
				calo->SetFillColor(kRed);	
			}
			else
			{
				//calo->SetLineColor(kGray);	
			}
			calo->SetLineColor(kBlack);	
			calo->SetLineWidth(2);
			/*
			if(om_side == 0 && om_column == 0) calo->Draw();
			else calo->Draw("same");
			*/
			calo->Draw("same");
			if(is_hit)
			{
				//tex->Draw("same");
			}
		}
	}
	
	// Drawing Xwall and Xwall calo hits
	for(int om_side = 0; om_side < 2; om_side++) 
	{
		for(int om_wall = 0; om_wall < 2; om_wall++) 
		{
			for(int om_column = 0; om_column < 2; om_column++) 
			{
				int swcr[4] = {om_side, om_wall, om_column, 0};
				TKOMhit* ohit = new TKOMhit(swcr); 
								
				TBox *calo = new TBox(ohit->get_xyz('y') - xw_sizey/2.0, 
				          	    -(ohit->get_xyz('x') + xw_sizex/2.0), 
				          	      ohit->get_xyz('y') + xw_sizey/2.0, 
				          	    -(ohit->get_xyz('x') - xw_sizex/2.0));
				
				bool is_hit = false;
				for(int hit = 0; hit < OM_hits.size(); hit++)
				{
					if(om_side   == OM_hits[hit]->get_SWCR('s') && 
				           om_wall   == OM_hits[hit]->get_SWCR('w') && 
				   	   om_column == OM_hits[hit]->get_SWCR('c') )  
				   	{
				   		is_hit = true;
				   	}
				}
				if(is_hit)
				{
					calo->SetFillColor(kRed);	
				}
				else
				{
					//calo->SetLineColor(kGray);	
				}
				calo->SetLineColor(kBlack);	
				calo->SetLineWidth(2);
				calo->Draw("same");	
			}	
		}
	}
	
	// Drawing tracker and tracker hits
	for(int cell_num = 0; cell_num < 2034; cell_num++)
	{
		TKtrhit* thit = new TKtrhit(cell_num); 
		
		double radius = tc_radius;
		double sigma;
		bool is_hit = false;
		bool is_broken = false;
		bool is_associated = false;
		bool has_height = false;
		
		for(int hit = 0; hit < tr_hits.size(); hit++)
		{
			if(cell_num == tr_hits[hit]->get_cell_num())
			{
				is_hit = true;
				if( tr_hits[hit]->get_r() > 35.0 || tr_hits[hit]->get_r() == -1.0 )
				{
					is_broken = true;
				}
				else
				{
					radius = tr_hits[hit]->get_r();
					sigma = tr_hits[hit]->get_sigma_R();
				}	
				
				if( tr_hits[hit]->get_associated_track() != nullptr )
				{
					is_associated = true;
				}
				
				if( tr_hits[hit]->get_h() != 0.0 )
				{
					has_height = true;
				}		
				break;
			}
		}
		
		TEllipse* tracker_cell;
		if(is_hit)
		{
			if(is_broken)
			{
				tracker_cell = new TEllipse(thit->get_xy('y'), -thit->get_xy('x'), radius, radius);
				tracker_cell->SetLineWidth(1);
				switch(option)
				{
					case 0:
						break;
					case 1:
						tracker_cell->SetFillColor(kOrange);
						tracker_cell->Draw("same");	
						break;
					case 2:
						tracker_cell->SetFillColor(kOrange);
						tracker_cell->Draw("same");			
						break;
					case 3:
						tracker_cell->SetFillColor(kOrange);
						tracker_cell->Draw("same");			
						break;		
				}
			}
			else
			{
				tracker_cell = new TEllipse(thit->get_xy('y'), -thit->get_xy('x'), radius + sigma, radius + sigma);
				tracker_cell->SetLineWidth(0);	
				
				if(is_associated)
				{
					switch(option)
					{
						case 0:
							tracker_cell->SetFillColor(kRed);
							break;
						case 1:
							tracker_cell->SetFillColor(kRed);
							break;
						case 2:
							tracker_cell->SetFillColor(kGreen);			
							break;
						case 3:
							if(has_height)
							{
								tracker_cell->SetFillColor(kGreen);			
							}
							else
							{							
								tracker_cell->SetFillColor(kTeal);
							}
							break;		
					}
				} 
				else
				{
					switch(option)
					{
						case 3:
						if(has_height)
						{
							tracker_cell->SetFillColor(kRed);
						}
						else
						{							
							tracker_cell->SetFillColor(kMagenta);
						}
						break;
						
						default:
							tracker_cell->SetFillColor(kRed);
						break;
					}
				}
				
				tracker_cell->Draw("same");
				if( radius - sigma > 0.0 )
				{
					TEllipse* tracker_cell_in = new TEllipse(thit->get_xy('y'), -thit->get_xy('x'), radius - sigma, radius - sigma);	
					tracker_cell_in->SetLineWidth(0);
					tracker_cell_in->Draw("same");
				}
			}
		}
		else
		{
			tracker_cell = new TEllipse(thit->get_xy('y'), -thit->get_xy('x'), radius, radius);
			tracker_cell->SetLineWidth(1);
			tracker_cell->Draw("same");
		}
		
	}
	
	// Drawing Bi sources
	for(int column = 0; column < 6; column++)
	{
		TBox *Bi_source = new TBox((column-2.5)*Bi_source_dist_y - Bi_source_y/2.0, -(-Bi_source_x/2.0), (column-2.5)*Bi_source_dist_y + Bi_source_y/2.0, -(Bi_source_x/2.0));
			
		Bi_source->SetFillColor(kBlue);	
	
		Bi_source->SetLineWidth(2);
		Bi_source->Draw("same");
	}
	
	vector<TKtrack*> all_tracks = this->get_tracks();
	// Drawing tracks
	for (int i = 0; i < this->get_no_tracks(); i++)
	{
		
		TLine* track;
		double x, y;
		
		x = 435.0;
		if( all_tracks[i]->get_side() == 0) 
		{
			x = -x;
		}		
		y = all_tracks[i]->get_a()*x + all_tracks[i]->get_b();
		
		if( y > 2505.5 )
		{
			x = (2505.5-all_tracks[i]->get_b())/all_tracks[i]->get_a();
			y = all_tracks[i]->get_a()*x + all_tracks[i]->get_b();
			
		}
		else if( y < -2505.5 )
		{
			x = (-2505.5-all_tracks[i]->get_b())/all_tracks[i]->get_a();
			y = all_tracks[i]->get_a()*x + all_tracks[i]->get_b();
		}
		
		track = new TLine(all_tracks[i]->get_b(), 0.0, y, -x);			
		track->SetLineColor(kBlue);
		track->SetLineWidth(3);
		track->Draw("same");
	}
	
	canvas->SaveAs(Form("./Events_visu/Run-%d_event-%d_2D.png", run_number, event_number));
	delete canvas;
}
		
void TKEvent::build_event()
{	
	gROOT->SetBatch(true);
	TCanvas *canv = new TCanvas();
	
	TGeoManager *geom = new TGeoManager();
	TGeoMaterial *matVacuum = new TGeoMaterial("matVacuum", 0, 0, 0);    
	TGeoMedium *Vacuum = new TGeoMedium("Vacuum",1, matVacuum);
	TGeoVolume *top = gGeoManager->MakeBox("top", Vacuum, 1500, 1750, 2500);
	geom->SetTopVolume(top);

	int object_counter = 0;

	// Drawing calorimeter
	for(int omnum = 0; omnum < 712; omnum++) 
	{
		// Skipping calo hits
		bool is_hit = false;
		for(int hit = 0; hit < OM_hits.size(); hit++)
		{
			if(omnum == OM_hits[hit]->get_OM_num()) is_hit = true;
		}
		if(is_hit) continue;
		
		TGeoVolume *calo;
		
		TKOMhit* ohit = new TKOMhit(omnum); 
		
		if(omnum < 520) // Mainwall
		{
			calo = gGeoManager->MakeBox(Form("OM:%d.%d.%d.%d", ohit->get_SWCR('s'), ohit->get_SWCR('w'), ohit->get_SWCR('c'), ohit->get_SWCR('r')), Vacuum, mw_sizex/2.0, mw_sizey/2.0, mw_sizez/2.0);
		}
		else if(omnum < 648) // Xwall
		{
			calo = gGeoManager->MakeBox(Form("OM:%d.%d.%d.%d", ohit->get_SWCR('s'), ohit->get_SWCR('w'), ohit->get_SWCR('c'), ohit->get_SWCR('r')), Vacuum, xw_sizex/2.0, xw_sizey/2.0, xw_sizez/2.0);
		}
		else // Gveto
		{
			calo = gGeoManager->MakeBox(Form("OM:%d.%d.%d.%d", ohit->get_SWCR('s'), ohit->get_SWCR('w'), ohit->get_SWCR('c'), ohit->get_SWCR('r')), Vacuum, gv_sizex/2.0, gv_sizey/2.0, gv_sizez/2.0);
		}	
		
		TGeoHMatrix *trans = new TGeoHMatrix("Trans");
		
		trans->SetDx(ohit->get_xyz('x'));
		trans->SetDy(ohit->get_xyz('y'));
		trans->SetDz(ohit->get_xyz('z'));
		
		calo->SetLineColor(kGray);	
		
		top->AddNode(calo, object_counter, trans);
		object_counter++;
	}

	// Drawing Bi calibration sources
	for(int row = 0; row < 7; row++)
	{    
	    for(int column = 0; column < 6; column++)
	    {
	    	TGeoVolume *Bi_source = gGeoManager->MakeBox("Bi_source", Vacuum, Bi_source_x/2.0, Bi_source_y/2.0, Bi_source_z/2.0);
	    	
	    	TGeoHMatrix *trans = new TGeoHMatrix("Trans");
		trans->SetDy( (column - 2.5) * Bi_source_dist_y );
		trans->SetDz( (row - 3.0) * Bi_source_dist_z );
		
		Bi_source->SetLineColor(kCyan);
		Bi_source->SetLineWidth(3);
		
		top->AddNode(Bi_source, object_counter, trans);
		object_counter++;
	    }
	}

	// Adding calo hits
	for(int hit = 0; hit < OM_hits.size(); hit++)
	{		
		TGeoVolume *box;
		if(OM_hits[hit]->get_OM_num() < 520) // Mainwall
		{
			box = gGeoManager->MakeBox(Form("OM:%d.%d.%d.%d", OM_hits[hit]->get_SWCR('s'), OM_hits[hit]->get_SWCR('w'), OM_hits[hit]->get_SWCR('c'), OM_hits[hit]->get_SWCR('r')), Vacuum, mw_sizex/2.0, mw_sizey/2.0, mw_sizez/2.0);
		}
		else if(OM_hits[hit]->get_OM_num() < 648) // Xwall
		{
			box = gGeoManager->MakeBox(Form("OM:%d.%d.%d.%d", OM_hits[hit]->get_SWCR('s'), OM_hits[hit]->get_SWCR('w'), OM_hits[hit]->get_SWCR('c'), OM_hits[hit]->get_SWCR('r')), Vacuum, xw_sizex/2.0, xw_sizey/2.0, xw_sizez/2.0);
		}
		else // Gveto
		{
			box = gGeoManager->MakeBox(Form("OM:%d.%d.%d.%d", OM_hits[hit]->get_SWCR('s'), OM_hits[hit]->get_SWCR('w'), OM_hits[hit]->get_SWCR('c'), OM_hits[hit]->get_SWCR('r')), Vacuum, gv_sizex/2.0, gv_sizey/2.0, gv_sizez/2.0);
		}	
		
		TGeoHMatrix *trans = new TGeoHMatrix("Trans");
		
		trans->SetDx(OM_hits[hit]->get_xyz('x'));
		trans->SetDy(OM_hits[hit]->get_xyz('y'));
		trans->SetDz(OM_hits[hit]->get_xyz('z'));
		
		if(OM_hits[hit]->is_HT())
		{
			box->SetLineColor(kRed);
			box->SetLineWidth(3);						
		}
		else
		{
			box->SetLineColor(kOrange);
			box->SetLineWidth(2);						
		}

		top->AddNode(box, object_counter, trans);
		object_counter++;
	}

	
 	// Adding tracker hits
	for(int hit = 0; hit < tr_hits.size(); hit++)
	{
		double radius = tc_radius; // default value
		double sigma_R = 2.0; //default value
		double sigma_Z = 17.0; //tr_hits[hit]->get_sigma_Z();
		bool is_broken = false;
		
		if( tr_hits[hit]->get_r() > 35.0 || tr_hits[hit]->get_r() == -1.0 ) 
		{
			is_broken = true;
		}
		
		if( !is_broken )
		{
			radius = tr_hits[hit]->get_r();
			sigma_R = tr_hits[hit]->get_sigma_R();
		}
		
		double radius_min = radius - sigma_R; 
		if( radius_min <= 0 )
		{
			radius_min = 0.0;
		}
	
		TGeoVolume *tracker_cell = geom->MakeTube(Form("cell:%d.%d.%d", tr_hits[hit]->get_SRL('s'), tr_hits[hit]->get_SRL('r'), tr_hits[hit]->get_SRL('l')), Vacuum, radius_min, radius + sigma_R, sigma_Z);
		TGeoHMatrix *trans = new TGeoHMatrix("Trans");
		    	
		trans->SetDx( tr_hits[hit]->get_xy('x') );
		trans->SetDy( tr_hits[hit]->get_xy('y') );
		trans->SetDz( tr_hits[hit]->get_h() );
		
		if( is_broken )
		{
		    	tracker_cell->SetLineColor(kOrange);
			tracker_cell->SetLineWidth(1);		
		}
		else if( tr_hits[hit]->get_h() == 0.0 )
		{
			tracker_cell->SetLineColor(kMagenta);
			tracker_cell->SetLineWidth(1);
		}
		else
		{
			tracker_cell->SetLineColor(kRed);
			tracker_cell->SetLineWidth(1);
		}
		
		top->AddNode(tracker_cell, object_counter, trans);
		object_counter++;
	}
	
	geom->CloseGeometry();     

	TFile *file = new TFile(Form("./Events_visu/Run-%d_event-%d_3D.root", run_number, event_number), "RECREATE");
	file->WriteObject(top, "demonstrator");

	vector<TKtrack*> all_tracks = this->get_tracks(); 
	for(int i = 0; i < this->get_no_tracks(); i++)
	{
		TPolyLine3D *track = new TPolyLine3D();
		track->SetPoint(0, 0.0, all_tracks[i]->get_b(), all_tracks[i]->get_d());
		double x,y,z;
		
		x = 435.0;
		if( all_tracks[i]->get_side() == 0 ) 
		{
			x = -x;
		}		
		y = all_tracks[i]->get_a()*x + all_tracks[i]->get_b();
		
		if( y > 2505.5 )
		{
			x = (2505.5-all_tracks[i]->get_b())/all_tracks[i]->get_a();
			y = all_tracks[i]->get_a()*x + all_tracks[i]->get_b();
			
		}
		else if( y < -2505.5 )
		{
			x = (-2505.5-all_tracks[i]->get_b())/all_tracks[i]->get_a();
			y = all_tracks[i]->get_a()*x + all_tracks[i]->get_b();
		}
		z = all_tracks[i]->get_c()*x + all_tracks[i]->get_d();
		
		if( z > 1550.0 )
		{
			x = (1550.0-all_tracks[i]->get_d())/all_tracks[i]->get_c();
			y = all_tracks[i]->get_a()*x + all_tracks[i]->get_b();
			z = all_tracks[i]->get_c()*x + all_tracks[i]->get_d();
		}
		else if( z < -1550.0 )
		{
			x = (-1550.0-all_tracks[i]->get_d())/all_tracks[i]->get_c();
			y = all_tracks[i]->get_a()*x + all_tracks[i]->get_b();
			z = all_tracks[i]->get_c()*x + all_tracks[i]->get_d();
		}
		track->SetPoint(1, x, y, z);			

		track->SetLineColor(kBlue);
		track->SetLineWidth(2);
		file->WriteObject(track, Form("track-%d", i));
	}

	file->Close();			
}
