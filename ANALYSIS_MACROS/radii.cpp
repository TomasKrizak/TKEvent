#include "../TKEvent/include/TKEvent.h" 

R__LOAD_LIBRARY(../TKEvent/lib/libTKEvent.so);

void radii()
{
	int run_number;
	cout << "enter run number: ";
	cin >> run_number;
	
	int lower_limit, upper_limit;
	cout << "enter starting event number: ";
	cin >> lower_limit;
		
	cout << "enter final event number: (-1 for entire run)";
	cin >> upper_limit;
	
	TFile* file_data = new TFile(Form("../runs/Run-%d.root", run_number));
	TTree* tree = (TTree*) file_data->Get("Event");

	if( upper_limit == -1 ) 
	{	
		upper_limit = tree->GetEntries();
	}

	TKEvent* event = new TKEvent();
	tree->SetBranchAddress("Eventdata", &event);

	cout << "Run number " << run_number << ", " << tree->GetEntries() << " events available." << endl << endl;


	TFile *file = new TFile(Form("tracker_hits_map_chi-%d.root", run_number), "RECREATE");
	
	
	TH2F* all_hits = new TH2F("all_hits", "all_hits", 5200, -2600.0, 2600.0, 1000, -500.0 , 500.0); 
	TH2F* rec_hits = new TH2F("rec_hits", "rec_hits", 5200, -2600.0, 2600.0, 1000, -500.0 , 500.0); 
	

	for(UInt_t ev = lower_limit; ev < upper_limit; ev++)	// Loop over events
	{
		tree->GetEntry(ev);
		event->set_r("Manchester","distance");
		event->reconstruct_ML_3D(0);	
		
		for(int i = 0; i < event->get_tr_hits().size(); i++)
		{		
			double radius = event->get_tr_hit(i)->get_r();
			if(radius == -1) continue;
			double x = event->get_tr_hit(i)->get_xy('x');
			double y = event->get_tr_hit(i)->get_xy('y');
			
			double dist;
			double x_j, y_k;
			for(int j = 0; j < 44; j++)
			{
				x_j = (double)j*1.0 - 21.5;
				for(int k = 0; k < 44; k++)
				{
					y_k = (double)k*1.0 - 21.5;
					dist = abs(sqrt(x_j*x_j + y_k*y_k) - radius);
					if( dist < sqrt(2.0)/2.0 )
					{
						all_hits->Fill(y+y_k, -(x+x_j));
					}
				}
			}
		}
		
		for(int i = 0; i < event->get_tracks().size(); i++)
		{
			TKtrack* track = event->get_tracks()[i];
			if( track->get_chi_squared() < 0.1 ) continue;
			if( track->get_chi_squared() > 1.0 ) continue;
			double a = track->get_a();
			double b = track->get_b();
			double c = track->get_c();
			double d = track->get_d();
			int side = track->get_side();
			
			
			if( track->get_ambiguity_type() == 2 )
			{
				if(abs(b- 417.5) > 25 && abs(b+ 417.5) > 25 && 
				   abs(b-1252.5) > 25 && abs(b+1252.5) > 25 &&
				   abs(b-2087.5) > 25 && abs(b+2087.5) > 25) continue;
			}
			if( track->get_ambiguity_type() == 3 )
			{
				 continue;
			}
			if( track->get_ambiguity_type() == 4 )
			{
				 continue;
			}
			
			//if(c == 0 && d == 0) continue;
			for(int i = 0; i < track->get_associated_tr_hits().size(); i++)
			{		
				TKtrhit* hit = track->get_associated_tr_hits()[i];
				double radius = hit->get_r();
				if(radius == -1) continue;
				double x = hit->get_xy('x');
				double y = hit->get_xy('y');
				
				double dist;
				double x_j, y_k;
				for(int j = 0; j < 44; j++)
				{
					x_j = (double)j*1.0 - 21.5;
					for(int k = 0; k < 44; k++)
					{
						y_k = (double)k*1.0 - 21.5;
						dist = abs(sqrt(x_j*x_j + y_k*y_k) - radius);
						if( dist < sqrt(2.0)/2.0 )
						{
							rec_hits->Fill(y+y_k, -(x+x_j));
						}
					}
				}
			}
		}	
			
		if(ev % 1000 == 0) std::cout << ev << std::endl;
	}
	/*
	for (int i = -2486; i < 2487; i=+44) 
	{
		TLine *vLine = new TLine(double(i), -424.0, double(i), 424.0);
		vLine->SetLineWidth(2);
		vLine->SetLineStyle(1); // Set line style
		vLine->Draw();
	}
	
	for (int i = -424; i < 0 ; i=+44) 
	{
		TLine *hLine = new TLine(-2486.0, double(i), 2486.0, double(i));
		hLine->SetLineWidth(2);
		hLine->SetLineStyle(1); // Set line style
		hLine->Draw();
	}
	for (int i = 28; i < 425 ; i=+44) 
	{
		TLine *hLine = new TLine(-2486.0, double(i), 2486.0, double(i));
		hLine->SetLineWidth(2);
		hLine->SetLineStyle(1); // Set line style
		hLine->Draw();
	}*/
	
	file_data->Close();
	file->WriteObject(rec_hits, "rec_hits");
	file->WriteObject(all_hits, "all_hits");
	file->Close();
	
}






