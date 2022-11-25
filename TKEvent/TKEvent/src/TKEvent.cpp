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
static double foil_spacex = 58.0; // probably wrong
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
	OM_hits = std::vector<TKOMhit*>();
	tr_hits = std::vector<TKtrhit*>();
}

TKEvent::TKEvent(int _run_number ,int _event_number)
{
	event_number = _event_number;
	run_number = _run_number;

	OM_hits = std::vector<TKOMhit*>();
	tr_hits = std::vector<TKtrhit*>();
}

TKEvent::~TKEvent()
{
	OM_hits.clear();
	tr_hits.clear();
	tracks.clear();
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
	return tracks;
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
	return tracks[_i];
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
	return tracks.size();
}
		
void TKEvent::print()
{
	std::cout << "*** Run: " << run_number << ", event: " << event_number << " ***" << std::endl << std::endl;
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
}

void TKEvent::print_tracks()
{
	std::cout << "*** Run: " << run_number << ", event: " << event_number << " ***" << std::endl;
		
	for (int i = 0; i < tracks.size(); i++)
	{
		tracks[i]->print();
	}
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

void TKEvent::reconstruct_track(bool save_sinograms)
{
	// first basic tracking - work in progress
	const int resolution = 250;

	for(int side = 0; side < 2; side++)
	{
		vector<double> hits_x;
		vector<double> hits_y;
		vector<double> hits_z;
		vector<double> hits_r;	
		
		for(int i = 0; i < tr_hits.size(); i++)
		{
			if( side == tr_hits[i]->get_SRL('s'))
			{
				// not using broken or too big (wrongly associated) tracker hits
				if( tr_hits[i]->get_r() != -1.0 && tr_hits[i]->get_r() < 35.0 )
				{
					hits_x.push_back( tr_hits[i]->get_xy('x') );
					hits_y.push_back( tr_hits[i]->get_xy('y') );
					hits_z.push_back( tr_hits[i]->get_h() );
					hits_r.push_back( tr_hits[i]->get_r() );
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
		
		for(int q = 0; q < 2; q++)
		{	
			double offset = (phi2 - phi1)/(2.0*resolution);
			TH2F *sinograms = new TH2F("sinograms", "sinograms; theta; r", resolution, phi1+offset, phi2+offset, resolution, R1, R2);		
			for(int i = 0; i < hits_x.size(); i++)
			{
				for(int k = 0; k < resolution; k++)
				{
					theta = phi1 + (k * (phi2 - phi1) / resolution);
					r = ( -hits_x[i]*cos(theta) ) - ( -hits_y[i]*sin(theta) );
					sinograms->Fill( theta, r - hits_r[i] );
					sinograms->Fill( theta, r + hits_r[i] );
				}	
			}
								
			double maximum = 0.0;
			for(int i = 1; i < resolution; i++)
			{
				for(int j = 1; j < resolution; j++)
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
				sinograms->Draw("COLZ");
				c2->SaveAs(Form("sinograms-run-%d_event-%d_side-%d_zoom-%d.png", run_number, event_number, side, q));
				c2->Close();
			}
			delete sinograms;
		}
			
		double a = 1.0 / tan(peak_Theta);
		double b = peak_R / sin(peak_Theta);
		
			
		// associating hits to track
		vector<bool> hits_associated;
		double denominator = sqrt((a*a) + 1);	
		double distance_from_wire;
		for(int i = 0; i < hits_x.size(); i++)
		{
			distance_from_wire = abs(hits_y[i] - a*hits_x[i] - b) / denominator;
			if( abs(distance_from_wire - hits_r[i]) < 10.0 )
			{
				hits_associated.push_back(true);
			}
			else
			{
				hits_associated.push_back(false);
			}
		
		}					
		
		// fitting z coordinates of hits
		vector<double> hits_p;
		double projection_distance;
		for(int i = 0; i < hits_x.size(); i++)
		{
			distance_from_wire = abs(hits_y[i] - a*hits_x[i] - b) / denominator;
			projection_distance = pow(hits_x[i], 2.0) + pow((hits_y[i] - b), 2.0) - pow(distance_from_wire, 2.0);
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
		
		tracks.push_back(new TKtrack(side, a, b, c, d));
		
		hits_x.erase(hits_x.begin(), hits_x.end());
		hits_y.erase(hits_y.begin(), hits_y.end());
		hits_z.erase(hits_z.begin(), hits_z.end());
		hits_r.erase(hits_r.begin(), hits_r.end());
	}
}
		
void TKEvent::make_top_projection()
{
	gROOT->SetBatch(true);
	TCanvas *canvas = new TCanvas("canvas","", 2000, 500);	
	canvas->Range(-2900.0, -700.0, 2900.0, 700.0);

	// Drawing mainwall and mainwall calo hits
	for(int om_side = 0; om_side < 2; om_side++) 
	{
		for(int om_column = 0; om_column < 20; om_column++) 
		{
			int swcr[4] = {om_side, -1, om_column, 0};
			TKOMhit* ohit = new TKOMhit(swcr); 
			
			TBox *calo = new TBox(ohit->get_xyz('y') - mw_sizey/2.0, 
					    -(ohit->get_xyz('x') - mw_sizex/2.0), 
					      ohit->get_xyz('y') + mw_sizey/2.0, 
					    -(ohit->get_xyz('x') + mw_sizex/2.0));
			
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
			if(om_side == 0 && om_column == 0) calo->Draw();
			else calo->Draw("same");
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
				          	    -(ohit->get_xyz('x') - xw_sizex/2.0), 
				          	      ohit->get_xyz('y') + xw_sizey/2.0, 
				          	    -(ohit->get_xyz('x') + xw_sizex/2.0));
				
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
		bool is_hit = false;
		bool is_broken = false;
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
				}						
				break;
			}
		}
		
		TEllipse *tracker_cell = new TEllipse(thit->get_xy('y'), -thit->get_xy('x'), radius, radius);
		if(is_hit)
		{
			if(is_broken)
			{
				tracker_cell->SetFillColor(kOrange);
				tracker_cell->SetLineWidth(1);					
			}
			else
			{
				tracker_cell->SetFillColor(kRed);
				tracker_cell->SetLineWidth(1);	
			}
		}
		else
		{
			tracker_cell->SetLineWidth(1);
		}
		
		tracker_cell->Draw("same");
	}
	
	// Drawing Bi sources
	for(int column = 0; column < 6; column++)
	{
		TBox *Bi_source = new TBox((column-2.5)*Bi_source_dist_y - Bi_source_y/2.0, -(-Bi_source_x/2.0), (column-2.5)*Bi_source_dist_y + Bi_source_y/2.0, -(Bi_source_x/2.0));
			
		Bi_source->SetFillColor(kBlue);	
	
		Bi_source->SetLineWidth(3);
		Bi_source->Draw("same");
	}
	
	// Drawing of tracks
	for (int i = 0; i < tracks.size(); i++)
	{
		TLine *track;
		double x, y;
		
		x = 435.0;
		if( tracks[i]->get_side() == 0) 
		{
			x = -x;
		}		
		y = tracks[i]->get_a()*x + tracks[i]->get_b();
		
		if( y > 2505.5 )
		{
			x = (2505.5-tracks[i]->get_b())/tracks[i]->get_a();
			y = tracks[i]->get_a()*x + tracks[i]->get_b();
			
		}
		else if( y < -2505.5 )
		{
			x = (-2505.5-tracks[i]->get_b())/tracks[i]->get_a();
			y = tracks[i]->get_a()*x + tracks[i]->get_b();
		}
		
		track = new TLine(tracks[i]->get_b(), 0.0, y, -x);			
		track->SetLineColor(kBlue);
		track->SetLineWidth(1);
		track->Draw("same");
	}
	
	canvas->SaveAs(Form("./Events_visu/Run-%d_event-%d_2D.png", run_number, event_number));
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
		double radius = tc_radius;
		bool is_broken = false;
		
		if( tr_hits[hit]->get_r() > 35.0 || tr_hits[hit]->get_r() == -1.0 ) 
		{
			is_broken = true;
		}
		
		if( !is_broken )
		{
			radius = tr_hits[hit]->get_r();
		}
		
		TGeoVolume *tracker_cell = geom->MakeTube(Form("cell:%d.%d.%d", tr_hits[hit]->get_SRL('s'), tr_hits[hit]->get_SRL('r'), tr_hits[hit]->get_SRL('l')), Vacuum, radius, radius+0.1, 0.0);
		TGeoHMatrix *trans = new TGeoHMatrix("Trans");
		    	
		trans->SetDx( tr_hits[hit]->get_xy('x') );
		trans->SetDy( tr_hits[hit]->get_xy('y') );
		trans->SetDz( tr_hits[hit]->get_h() );
		
		if( is_broken )
		{
		    	tracker_cell->SetLineColor(kOrange);
			tracker_cell->SetLineWidth(2);		
		}
		else
		{
			tracker_cell->SetLineColor(kRed);
			tracker_cell->SetLineWidth(2);
		}
		
		top->AddNode(tracker_cell, object_counter, trans);
		object_counter++;
	}
	
	geom->CloseGeometry();     

	TFile *file = new TFile(Form("./Events_visu/Run-%d_event-%d_3D.root", run_number, event_number), "RECREATE");
	file->WriteObject(top, "demonstrator");

	for(int i = 0; i < tracks.size(); i++)
	{
		TPolyLine3D *track = new TPolyLine3D();
		track->SetPoint(0, 0.0, tracks[i]->get_b(), tracks[i]->get_d());
		double x,y,z;
		
		x = 435.0;
		if( tracks[i]->get_side() == 0 ) 
		{
			x = -x;
		}		
		y = tracks[i]->get_a()*x + tracks[i]->get_b();
		
		if( y > 2505.5 )
		{
			x = (2505.5-tracks[i]->get_b())/tracks[i]->get_a();
			y = tracks[i]->get_a()*x + tracks[i]->get_b();
			
		}
		else if( y < -2505.5 )
		{
			x = (-2505.5-tracks[i]->get_b())/tracks[i]->get_a();
			y = tracks[i]->get_a()*x + tracks[i]->get_b();
		}
		z = tracks[i]->get_c()*x + tracks[i]->get_d();
		
		if( z > 1550.0 )
		{
			x = (1550.0-tracks[i]->get_d())/tracks[i]->get_c();
			y = tracks[i]->get_a()*x + tracks[i]->get_b();
			z = tracks[i]->get_c()*x + tracks[i]->get_d();
		}
		else if( z < -1550.0 )
		{
			x = (-1550.0-tracks[i]->get_d())/tracks[i]->get_c();
			y = tracks[i]->get_a()*x + tracks[i]->get_b();
			z = tracks[i]->get_c()*x + tracks[i]->get_d();
		}
		track->SetPoint(1, x, y, z);			

		track->SetLineColor(kBlue);
		track->SetLineWidth(1);
		file->WriteObject(track, Form("track-%d", i));
	}

	file->Close();			
}

void TKEvent::set_r(std::string _model_n)
{
	
	for(int tr_hit = 0; tr_hit < tr_hits.size(); tr_hit++)
	{
		double r;
		if(tr_hits[tr_hit]->get_tsp('0') != -1)
		{
			double min_time = 1e10;
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

			if(_model_n == "Manchester")
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
			else if(_model_n == "Betsy")
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
			else cout << "invalid drift time model: choose \"Betsy\" or \"Manu\"." << endl;
		}
		else
		{
			r = -1.0;
		}
		
		tr_hits[tr_hit]->set_r(r);
	}
}
