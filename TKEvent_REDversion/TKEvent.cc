#include "TKEvent.hh"

using namespace std;

TKEvent::TKEvent(int _run_number ,int _event_number)
{
	event_number = _event_number;
	run_number = _run_number;
}

TKEvent::~TKEvent()
{
	hit_om_num.clear();
	hit_is_HT.clear();
	hit_om_tdc.clear();
	hit_om_peak_cell.clear();
	
	hit_cell_num.clear();
	hit_height.clear();
	hit_radius.clear();
	
	rec_tracks_side.clear();
	rec_tracks_a.clear();
	rec_tracks_b.clear();
	rec_tracks_c.clear();
	rec_tracks_d.clear();
}

int TKEvent::get_run_number()
{
	return run_number;
}

int TKEvent::get_event_number()
{	
	return event_number;
}

int TKEvent::get_track_side(int index)
{
	return rec_tracks_side.at(index);
}

double TKEvent::get_track_a(int index)
{
	return rec_tracks_a.at(index);
}

double TKEvent::get_track_b(int index)
{
	return rec_tracks_b.at(index);
}

double TKEvent::get_track_c(int index)
{
	return rec_tracks_c.at(index);
}

double TKEvent::get_track_d(int index)
{
	return rec_tracks_d.at(index);
}

int TKEvent::get_no_tracks()
{
	return rec_tracks_a.size();
}

void TKEvent::add_event(snemo::datamodel::raw_event_data red)
{
	const std::vector<snemo::datamodel::calo_digitized_hit> red_calo_hits = red.get_calo_hits();
	for (const snemo::datamodel::calo_digitized_hit & red_calo_hit : red_calo_hits)
	{
		add_calo_hit(red_calo_hit);
	}

	const std::vector<snemo::datamodel::tracker_digitized_hit> red_tracker_hits = red.get_tracker_hits();
	for (const snemo::datamodel::tracker_digitized_hit & red_tracker_hit : red_tracker_hits)
	{
	  	add_tracker_hit(red_tracker_hit);
	} 
}

void TKEvent::add_calo_hit(const snemo::datamodel::calo_digitized_hit & red_calo_hit)
{
	if(red_calo_hit.is_high_threshold())
	{
		hit_is_HT.push_back(true);
	}
	else if(red_calo_hit.is_low_threshold_only())
	{
		hit_is_HT.push_back(false);
	}
	else return;	
	
	const snemo::datamodel::timestamp & reference_time = red_calo_hit.get_reference_time();
	const int64_t calo_tdc = reference_time.get_ticks();
	const int16_t peak_cell = red_calo_hit.get_fwmeas_peak_cell();
	
	hit_om_tdc.push_back( calo_tdc );
	hit_om_peak_cell.push_back( peak_cell );
	
	const sncabling::om_id om_id = red_calo_hit.get_om_id();
	int om_side, om_wall, om_column, om_row, om_num;
	if(om_id.is_main())
	{
		om_side   = om_id.get_side();
		om_column = om_id.get_column();
		om_row    = om_id.get_row();
		om_num = om_side*20*13 + om_column*13 + om_row;
	}

	else if(om_id.is_xwall())
	{
		om_side   = om_id.get_side();
		om_wall   = om_id.get_wall();
		om_column = om_id.get_column();
		om_row    = om_id.get_row();
		om_num = 520 + om_side*64 +  om_wall*32 + om_column*16 + om_row;
	}

	else if(om_id.is_gveto())
	{
		om_side = om_id.get_side();
		om_wall = om_id.get_wall();
		om_column = om_id.get_column();
		om_num = 520 + 128 + om_side*32 + om_wall*16 + om_column;
	}
	
	hit_om_num.push_back(om_num);
}

void TKEvent::add_calo_hit(int _hit_om_num, bool _hit_is_HT, int64_t calo_tdc, int16_t peak_cell)
{
	hit_om_num.push_back(_hit_om_num);
	hit_is_HT.push_back(_hit_is_HT);
	hit_om_tdc.push_back(calo_tdc);
	hit_om_peak_cell.push_back(peak_cell);
}
		
void TKEvent::add_calo_hit(int cell_side, int cell_wall, int cell_column, int cell_row, bool is_HT, int64_t calo_tdc, int16_t peak_cell)
{
	add_calo_hit(SWCR_to_OMnum(cell_side, cell_wall, cell_column, cell_row), is_HT, calo_tdc, peak_cell);
}

void TKEvent::add_tracker_hit(const snemo::datamodel::tracker_digitized_hit & red_tracker_hit)
{
	const sncabling::gg_cell_id gg_id = red_tracker_hit.get_cell_id();
	int cell_side  = gg_id.get_side();
	int cell_row   = gg_id.get_row();
	int cell_layer = gg_id.get_layer();
	hit_cell_num.push_back(113*9*cell_side + 9*cell_row + cell_layer);

	double height = tc_sizez;
	double radius = tc_radius;
	for(const snemo::datamodel::tracker_digitized_hit::gg_times & gg_timestamps : red_tracker_hit.get_times())
	{
		height = get_tracker_hit_height(gg_timestamps);
	
		const int64_t anode_tdc = gg_timestamps.get_anode_time(0).get_ticks();
		if(anode_tdc != snfee::data::INVALID_TICKS)
		{
			double min_time = 1e10;
			//double calo_hit = (calo_tdc * 6.25) - 400.0 + (400.0 * peak_cell / 1024.0);
			
			for(int om_hit = 0; om_hit < hit_om_num.size(); om_hit++)
			{
				if(2*anode_tdc - (hit_om_tdc[om_hit]-44) < 800 && 2*anode_tdc - (hit_om_tdc[om_hit]-44) > 0 )
				{
					if(min_time > double(2*anode_tdc - (hit_om_tdc[om_hit]-44)))
					{
						min_time = 6.25 * (2*anode_tdc - (hit_om_tdc[om_hit]-44));
					}
				}
			}
			radius = drift_time_to_radius(min_time, "Manu");
		}		
	}
	hit_radius.push_back(radius);
	hit_height.push_back(height);
}
		
void TKEvent::add_tracker_hit(int _hit_cell_num, double _hit_height, double _hit_radius)
{
	hit_cell_num.push_back(_hit_cell_num);
	hit_height.push_back(_hit_height);
	hit_radius.push_back(_hit_radius);
}
		
void TKEvent::add_tracker_hit(int hit_cell_side, int hit_cell_row, int hit_cell_layer, double hit_height, double hit_radius)
{
	add_tracker_hit(SRL_to_Cellnum(hit_cell_side, hit_cell_row, hit_cell_layer), hit_height, hit_radius);
}

void TKEvent::reconstruct_track()
{
	// first basic tracking - work in progress
	const int resolution = 250;

	for(int side = 0; side < 2; side++)
	{
		vector<double> hits_x;
		vector<double> hits_y;
		vector<double> hits_z;
		vector<double> hits_r;	
		
		for(int i = 0; i < hit_cell_num.size(); i++)
		{
			int hit_SRL[3];
			Cellnum_to_SRL(hit_cell_num[i] , hit_SRL);
			
			if( side == hit_SRL[0] )
			{
				// not using broken or too big (wrongly associated) tracker hits
				if( hit_radius[i] != tc_radius && hit_radius[i] < 35.0 )
				{
					hits_x.push_back( (2.0*hit_SRL[0]-1.0) * (hit_SRL[2]*2.0*tc_radius + tc_radius + foil_spacex/2.0) );
					hits_y.push_back( (hit_SRL[1]-56.0) * 2.0 * tc_radius );
					hits_z.push_back( hit_height[i] );
					hits_r.push_back( hit_radius[i] );
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


			TCanvas* c2 = new TCanvas("sinograms");
			sinograms->Draw("COLZ");
			//c2->SaveAs(Form("sinograms-run-%d_event-%d_side-%d_zoom-%d.png", run_number, event_number, side, q));
			c2->Close();
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
			if( abs(distance_from_wire - hits_r[i]) < 5.0 )
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
		
		
		rec_tracks_side.push_back(side);
		rec_tracks_a.push_back(a);
		rec_tracks_b.push_back(b);
		rec_tracks_c.push_back(c);
		rec_tracks_d.push_back(d);
		
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
			
			double calo_position[3];
			SWCR_to_xyz(om_side, -1, om_column, 0, calo_position);
			
			TBox *calo = new TBox(calo_position[1] - mw_sizey/2.0, -(calo_position[0] - mw_sizex/2.0), calo_position[1] + mw_sizey/2.0, -(calo_position[0] + mw_sizex/2.0));
			
			bool is_hit = false;
			for(int hit = 0; hit < hit_om_num.size(); hit++)
			{
				int SWCR_hit[4];
				OMnum_to_SWCR(hit_om_num[hit], SWCR_hit);
				if(om_side == SWCR_hit[0] && -1 == SWCR_hit[1] && om_column == SWCR_hit[2]) is_hit = true;
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
				double calo_position[3];
				SWCR_to_xyz(om_side, om_wall, om_column, 0, calo_position);
				
				TBox *calo = new TBox(calo_position[1] - xw_sizey/2.0, -(calo_position[0] - xw_sizex/2.0), calo_position[1] + xw_sizey/2.0, -(calo_position[0] + xw_sizex/2.0));
				
				bool is_hit = false;
				for(int hit = 0; hit < hit_om_num.size(); hit++)
				{
					int SWCR_hit[4];
					OMnum_to_SWCR(hit_om_num[hit], SWCR_hit);
					if(om_side == SWCR_hit[0] && om_wall == SWCR_hit[1] && om_column == SWCR_hit[2]) is_hit = true;
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
		double cell_position[2];
		Cellnum_to_xy(cell_num, cell_position);
		
		double radius = tc_radius;
		bool is_hit = false;
		bool is_broken = false;
		for(int hit = 0; hit < hit_cell_num.size(); hit++)
		{
			if(cell_num == hit_cell_num[hit])
			{
				is_hit = true;
				if( hit_radius[hit] > 35.0 || hit_radius[hit] == tc_radius )
				{
					is_broken = true;
				}
				else
				{
					radius = hit_radius[hit];
				}						
				break;
			}
		}
		
		TEllipse *tracker_cell = new TEllipse(cell_position[1], -cell_position[0], radius, radius);
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
	for (int i = 0; i < rec_tracks_side.size(); i++)
	{
		TLine *track;
		double x, y;
		
		x = 435.0;
		if( rec_tracks_side[i] == 0) 
		{
			x = -x;
		}		
		y = rec_tracks_a[i]*x + rec_tracks_b[i];
		
		if( y > 2505.5 )
		{
			x = (2505.5-rec_tracks_b[i])/rec_tracks_a[i];
			y = rec_tracks_a[i]*x + rec_tracks_b[i];
			
		}
		else if( y < -2505.5 )
		{
			x = (-2505.5-rec_tracks_b[i])/rec_tracks_a[i];
			y = rec_tracks_a[i]*x + rec_tracks_b[i];
		}
		
		track = new TLine(rec_tracks_b[i], 0.0, y, -x);			
		track->SetLineColor(kBlue);
		track->SetLineWidth(1);
		track->Draw("same");
	}
	
	canvas->SaveAs(Form("Run-%d_event-%d_2D.png", run_number, event_number));
}
		
void TKEvent::build_event_3D()
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
		for(int hit = 0; hit < hit_om_num.size(); hit++)
		{
			if(omnum == hit_om_num[hit]) is_hit = true;
		}
		if(is_hit) continue;
		
		TGeoVolume *calo;
		
		double calo_position[3];
		int SWCR[4];
		OMnum_to_SWCR(omnum, SWCR);
		OMnum_to_xyz(omnum, calo_position);
		
		if(omnum < 520) // Mainwall
		{
			calo = gGeoManager->MakeBox(Form("OM:%d.%d.%d.%d", SWCR[0], SWCR[1], SWCR[2], SWCR[3]), Vacuum, mw_sizex/2.0, mw_sizey/2.0, mw_sizez/2.0);
		}
		else if(omnum < 648) // Xwall
		{
			calo = gGeoManager->MakeBox(Form("OM:%d.%d.%d.%d", SWCR[0], SWCR[1], SWCR[2], SWCR[3]), Vacuum, xw_sizex/2.0, xw_sizey/2.0, xw_sizez/2.0);
		}
		else // Gveto
		{
			calo = gGeoManager->MakeBox(Form("OM:%d.%d.%d.%d", SWCR[0], SWCR[1], SWCR[2], SWCR[3]), Vacuum, gv_sizex/2.0, gv_sizey/2.0, gv_sizez/2.0);
		}	
		
		TGeoHMatrix *trans = new TGeoHMatrix("Trans");
		
		trans->SetDx(calo_position[0]);
		trans->SetDy(calo_position[1]);
		trans->SetDz(calo_position[2]);
		
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
	for(int hit = 0; hit < hit_om_num.size(); hit++)
	{
		int omnum = hit_om_num[hit];
		double calo_position[3];
		int SWCR[4];
		OMnum_to_SWCR(omnum, SWCR);
		SWCR_to_xyz(SWCR[0], SWCR[1], SWCR[2], SWCR[3], calo_position);
		
		TGeoVolume *box;
		if(omnum < 520) // Mainwall
		{
			box = gGeoManager->MakeBox(Form("OM:%d.%d.%d.%d", SWCR[0], SWCR[1], SWCR[2], SWCR[3]), Vacuum, mw_sizex/2.0, mw_sizey/2.0, mw_sizez/2.0);
		}
		else if(omnum < 648) // Xwall
		{
			box = gGeoManager->MakeBox(Form("OM:%d.%d.%d.%d", SWCR[0], SWCR[1], SWCR[2], SWCR[3]), Vacuum, xw_sizex/2.0, xw_sizey/2.0, xw_sizez/2.0);
		}
		else // Gveto
		{
			box = gGeoManager->MakeBox(Form("OM:%d.%d.%d.%d", SWCR[0], SWCR[1], SWCR[2], SWCR[3]), Vacuum, gv_sizex/2.0, gv_sizey/2.0, gv_sizez/2.0);
		}	
		
		TGeoHMatrix *trans = new TGeoHMatrix("Trans");
		
		trans->SetDx(calo_position[0]);
		trans->SetDy(calo_position[1]);
		trans->SetDz(calo_position[2]);
		
		if(hit_is_HT[hit] == 1)
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
	for(int hit = 0; hit < hit_cell_num.size(); hit++)
	{
		int SRL[3];
		Cellnum_to_SRL(hit_cell_num[hit], SRL);	
		
		double radius = tc_radius;
		bool is_broken = false;
		
		if( hit_radius[hit] > 35.0 || hit_radius[hit] == tc_radius ) 
		{
			is_broken = true;
		}
		
		if( !is_broken )
		{
			radius = hit_radius[hit];
		}
		
		TGeoVolume *tracker_cell = geom->MakeTube(Form("cell:%d.%d.%d", SRL[0], SRL[1], SRL[2]), Vacuum, radius, radius+0.1, 0.0);
		TGeoHMatrix *trans = new TGeoHMatrix("Trans");
		
		double position_xy[2];
		Cellnum_to_xy(hit_cell_num[hit], position_xy);
    	
		trans->SetDx( position_xy[0] );
		trans->SetDy( position_xy[1] );
		trans->SetDz( hit_height[hit] );
		
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

	TFile *file = new TFile(Form("Run-%d_event-%d_3D.root", run_number, event_number), "RECREATE");
	file->WriteObject(top, "demonstrator");

	for(int i = 0; i < rec_tracks_side.size(); i++)
	{
		TPolyLine3D *track = new TPolyLine3D();
		track->SetPoint(0, 0.0, rec_tracks_b[i], rec_tracks_d[i]);
		double x,y,z;
		
		x = 435.0;
		if( rec_tracks_side[i] == 0 ) 
		{
			x = -x;
		}		
		y = rec_tracks_a[i]*x + rec_tracks_b[i];
		
		if( y > 2505.5 )
		{
			x = (2505.5-rec_tracks_b[i])/rec_tracks_a[i];
			y = rec_tracks_a[i]*x + rec_tracks_b[i];
			
		}
		else if( y < -2505.5 )
		{
			x = (-2505.5-rec_tracks_b[i])/rec_tracks_a[i];
			y = rec_tracks_a[i]*x + rec_tracks_b[i];
		}
		z = rec_tracks_c[i]*x + rec_tracks_d[i];
		
		if( z > 1550.0 )
		{
			x = (1550.0-rec_tracks_d[i])/rec_tracks_c[i];
			y = rec_tracks_a[i]*x + rec_tracks_b[i];
			z = rec_tracks_c[i]*x + rec_tracks_d[i];
		}
		else if( z < -1550.0 )
		{
			x = (-1550.0-rec_tracks_d[i])/rec_tracks_c[i];
			y = rec_tracks_a[i]*x + rec_tracks_b[i];
			z = rec_tracks_c[i]*x + rec_tracks_d[i];
		}
		track->SetPoint(1, x, y, z);			

		track->SetLineColor(kBlue);
		track->SetLineWidth(1);
		file->WriteObject(track, Form("track-%d", i));
	}

	file->Close();			
}



//-------------------------------------------------------------------------------------------------------------------------------------
  
double TKEvent::get_tracker_hit_height(const snemo::datamodel::tracker_digitized_hit::gg_times & gg_timestamps)
{
	const int64_t anode_tdc = gg_timestamps.get_anode_time(0).get_ticks();
	double hit_height = 0.0;
	if (anode_tdc != snfee::data::INVALID_TICKS)
	{
		const int64_t top_cathode_tdc = gg_timestamps.get_top_cathode_time().get_ticks();
		const int64_t bottom_cathode_tdc = gg_timestamps.get_bottom_cathode_time().get_ticks();
		if( top_cathode_tdc != snfee::data::INVALID_TICKS && bottom_cathode_tdc != snfee::data::INVALID_TICKS )
		{
			int64_t propagation_time_tdc = bottom_cathode_tdc + top_cathode_tdc - 2*anode_tdc;
			hit_height = tc_sizez*( double(bottom_cathode_tdc - anode_tdc) / double(propagation_time_tdc) ) - tc_sizez/2.0;
		}
	}
	//else cout << "missing r0" << endl;
	
	return hit_height;	
}

//-------------------------------------------------------------------------------------------------------------------------------------
  
// drift model
double TKEvent::drift_time_to_radius(double time_, string type) 
{
	double r = tc_radius;
	if(type == "Manu")
	{
		const double A1 = 0.570947153108633;
		const double B1 = 0.580148313540993;
		const double C1 = 1.6567483468611;
		const double A2 = 1.86938462695651;
		const double B2 = 0.949912427483918;
		const double t_usec = time_ / CLHEP::microsecond;
		const double ut = 10. * t_usec;
		
		r = A1 * ut / (std::pow(ut, B1) + C1);
		if (r > A2 * ut / (std::pow(ut, B2))) 
		{
			r = A2 * ut / (std::pow(ut, B2));
		}
		r *= CLHEP::cm;
		
	}
	// modification of Betsy's model
	else if(type == "Betsy")
	{
		const double a1 = 0.828;
		const double b1 = -0.907;
		const double a2 = 0.402;
		const double b2 = -1.955;
		
		time_ = time_ / 1000.0;
		if(time_ < 3.0845 && time_ > 0.0)
		{
			r = pow(time_/a1, 1.0/(1.0-b1)) * 10.0;
		}
		else
		{
			r = pow(time_/a2, 1.0/(1.0-b2)) * 10.0;
		}
	}
	else cout << "invalid drift time model: choose \"Betsy\" or \"Manu\"." << endl;

	return r;
 }

//-------------------------------------------------------------------------------------------------------------------------------------

// tracker: (side, row, layer) to x,y coordinates (centers of tracker cells)
void TKEvent::SRL_to_xy(int cell_side, int cell_row, int cell_layer, double position_xy[])
{
	position_xy[0] = (2.0*cell_side-1.0) * (cell_layer*2.0*tc_radius + tc_radius + foil_spacex/2.0);
	position_xy[1] = (cell_row-56.0) * 2.0 * tc_radius;
	
	return;
}

//-------------------------------------------------------------------------------------------------------------------------------------

// tracker: Cell number to x,y coordinates (centers of tracker cells)
void TKEvent::Cellnum_to_xy(int cell_num, double position_xy[])
{	
	int SRL[3];
	Cellnum_to_SRL(cell_num, SRL);
	SRL_to_xy(SRL[0], SRL[1], SRL[2], position_xy);
	
	return;
}
	
//-------------------------------------------------------------------------------------------------------------------------------------

// tracker: (side, row, layer) to Cell number
int TKEvent::SRL_to_Cellnum(int cell_side, int cell_row, int cell_layer)
{
	return 113*9*cell_side + 9*cell_row + cell_layer;
} 

//-------------------------------------------------------------------------------------------------------------------------------------

// tracker: Cell number to (side, row, layer)
void TKEvent::Cellnum_to_SRL(int cell_num, int SRL[])
{
	if( cell_num > 2033 || cell_num < 0) std::cout << "warning: " << cell_num << " is not valid tracker cell number" << std::endl; 
	
	SRL[0] = cell_num / 1017;
	SRL[1] = (cell_num % 1017)/ 9; 
	SRL[2] = cell_num % 9;
	return;
}

//-------------------------------------------------------------------------------------------------------------------------------------

// calo: (side, wall, column, row) to OM number
int TKEvent::SWCR_to_OMnum(int om_side, int om_wall, int om_column, int om_row)
{
	// auto detect MW
	if ((om_side!=-1) && (om_wall==-1) && (om_column!=-1) && (om_row!=-1))
		return 260*om_side + 13*om_column + om_row;

	// auto detect XW
	else if ((om_side!=-1) && (om_wall!=-1) && (om_column!=-1) && (om_row!=-1))
		return 520 + 64*om_side + 32*om_wall + 16*om_column + om_row;

	// auto detect GV
	else if ((om_side!=-1) && (om_wall!=-1) && (om_column!=-1) && (om_row==-1))
		return 520 + 128 + 32*om_side + 16*om_wall + om_column;

	else 
	{
		std::cout << "warning: " << om_side << "." << om_wall << "." << om_column << "." << om_row << " is not valid OM" << std::endl;
		return -1;
	}		
}
	
//-------------------------------------------------------------------------------------------------------------------------------------

// calo: OM number to (side, wall, column, row)
void TKEvent::OMnum_to_SWCR(int om_num, int OM[])
{
	//mainwall IT
	if(om_num < 260) 
	{
		OM[0] = 0;
		OM[1] = -1;
		OM[2] = om_num / 13;
		OM[3] = om_num % 13;
	}
	//mainwall FR
	else if(om_num < 520)
	{
		OM[0] = 1;
		OM[1] = -1;
		OM[2] = (om_num - 260) / 13;
		OM[3] = (om_num - 260) % 13;
	}
	//Xcalo IT	
	else if(om_num < 584)
	{
		OM[0] = 0;
		OM[1] = (om_num - 520) / 32;
		OM[2] = ((om_num - 520) / 16) % 2;
		OM[3] = (om_num -520) % 16;
	}
	//Xcalo FR
	else if(om_num < 648)
	{
		OM[0] = 1;
		OM[1] = (om_num - 520 - 64) / 32;
		OM[2] = ((om_num - 520 - 64) / 16) % 2;
		OM[3] = (om_num -520 - 64) % 16;
	}
	//GVeto IT
	else if(om_num < 680)
	{
		OM[0] = 0;
		OM[1] = (om_num - 520 - 128) / 16;
		OM[2] = (om_num - 520 - 128) % 16;
		OM[3] = -1;
	}
	//GVeto FR
	else if(om_num < 712)
	{
		OM[0] = 1;
		OM[1] = (om_num - 520 - 128 - 32) / 16;
		OM[2] = (om_num - 520 - 128 - 32) % 16;
		OM[3] = -1;
	}
	return;
}
    
//-------------------------------------------------------------------------------------------------------------------------------------

// calo: OM number to x,y,z coordinates (centers of calos)
// warning: dimensions are approximately taken from Falaise and might not be exactly right
void TKEvent::OMnum_to_xyz(int om_num, double position_xyz[])
{
	int SWCR[4];
	OMnum_to_SWCR(om_num, SWCR);
	SWCR_to_xyz(SWCR[0], SWCR[1], SWCR[2], SWCR[3], position_xyz);

	return;
}

//-------------------------------------------------------------------------------------------------------------------------------------

// calo: (side, wall, column, row) to x,y,z coordinates (centers of calos)
// warning: dimensions are approximately taken from Falaise and might not be exactly right
void TKEvent::SWCR_to_xyz(int om_side, int om_wall, int om_column, int om_row, double position_xyz[])
{
	int om_num = SWCR_to_OMnum(om_side, om_wall, om_column, om_row);
	int om_type;
	
	if(om_num < 520)
	{
		om_type = 1302;
	}
	else if(om_num < 648)
	{
		om_type = 1232;
	}
	else
	{
		om_type = 1252;
	}

	switch(om_type)
	{
		case 1302: //MW
			if(om_side == 1)
				position_xyz[0] = 532.0;
			else
				position_xyz[0] = -532.0;
			position_xyz[1] = ((double)om_column - 9.5) * 259.0;
			position_xyz[2] = ((double)om_row - 6) * 259.0;
			break;
			
		case 1232: //XW
			if(om_wall == 1)
				position_xyz[1] = 2580.5;
			else
				position_xyz[1] = -2580.5;
			if(om_side == 1)
			{
				if(om_column == 1)
					position_xyz[0] = 333.0;
				else
					position_xyz[0] = 130.0;
			}
			else
			{
				if(om_column == 1)
					position_xyz[0] = -333.0;
				else
					position_xyz[0] = -130.0;
			}
			position_xyz[2] = ((double)om_row - 7.5) * 212.0;
			break;
			
		case 1252: //GV
			if(om_side == 1)
				position_xyz[0] = 213.5;
			else
				position_xyz[0] = -213.5;
			if(om_wall == 1)
				position_xyz[2] = 1625.0;
			else
				position_xyz[2] = -1625.0;
			if(om_column > 7)
				position_xyz[1] = 161.0 + (((double)om_column-8) * 311.5);
			else
				position_xyz[1] = -161.0 + (((double)om_column-7) * 311.5);
			break;	
	}
	return;	
}

