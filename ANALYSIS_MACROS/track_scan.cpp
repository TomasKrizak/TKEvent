#include "../TKEvent/include/TKEvent.h" 

R__LOAD_LIBRARY(../TKEvent/lib/libTKEvent.so);

void track_scan()
{
	int run_number;
	cout << "enter run number: ";
	cin >> run_number;
	
	int lower_limit, upper_limit;
	cout << "enter starting event number: ";
	cin >> lower_limit;
		
	cout << "enter final event number: (-1 for entire run)";
	cin >> upper_limit;
	
	bool ambiguity_cut;
	cout << "Do you want to filter ambiguities? (good for calibration runs)" << endl;
	cout << "Enter 1 for yes, 0 for no: ";
	cin >> ambiguity_cut;
	
	TFile* file_data = new TFile(Form("../runs/Run-%d.root", run_number));
	TTree* tree = (TTree*) file_data->Get("Event");

	if( upper_limit == -1 ) 
	{	
		upper_limit = tree->GetEntries();
	}

	TKEvent* event = new TKEvent();
	tree->SetBranchAddress("Eventdata", &event);

	cout << "Run number " << run_number << ", " << tree->GetEntries() << " events available." << endl << endl;


	TFile *file = new TFile(Form("tracks_run-%d.root", run_number), "RECREATE");
	
	TH2F* YX_view = new TH2F("YX_view", "YX_view", 5200, -2600.0, 2600.0, 1000, -500.0 , 500.0); 
	TH2F* YZ_view = new TH2F("YZ_view", "YZ_view", 5200, -2600.0, 2600.0, 3000, -1500.0 , 1500.0); 
	TH2F* XZ_view = new TH2F("XZ_view", "XZ_view", 1000, -500.0, 500.0, 3000, -1500.0 , 1500.0); 


	for(UInt_t ev = lower_limit; ev < upper_limit; ev++)	// Loop over events
	{
		tree->GetEntry(ev);
		
		event->set_r("Manchester","distance");
		event->set_h();
		event->reconstruct_ML(0);	
		for(int i = 0; i < event->get_tracks().size(); i++)
		{
			if( event->get_tracks()[i]->get_chi_squared() < 0.1 ) continue;
			if( event->get_tracks()[i]->get_chi_squared() > 1.0 ) continue;
			double x,y,z;
			double a = event->get_tracks()[i]->get_a();
			double b = event->get_tracks()[i]->get_b();
			double c = event->get_tracks()[i]->get_c();
			double d = event->get_tracks()[i]->get_d();
			int side = event->get_tracks()[i]->get_side();
			
			if(ambiguity_cut == true)
			{	
				switch( event->get_tracks()[i]->get_ambiguity_type() )
				{
					case 0: 
						break;
					case 1: 
						if(abs(b- 417.5) > 25 && abs(b+ 417.5) > 25 && 
						   abs(b-1252.5) > 25 && abs(b+1252.5) > 25 &&
						   abs(b-2087.5) > 25 && abs(b+2087.5) > 25) continue;
						break;
					case 2:
						if(abs(b- 417.5) > 25 && abs(b+ 417.5) > 25 && 
						   abs(b-1252.5) > 25 && abs(b+1252.5) > 25 &&
						   abs(b-2087.5) > 25 && abs(b+2087.5) > 25) continue;
						break;
					case 3: 
						continue;
						break;
					case 4:
						continue;
						break;
				}
			}
			

			if(c == 0 && d == 0) continue;
			
			if(side == 1)
			{			
				if( abs(a) <= 1.0 )
				{
					for(int j = 0; j < 500; j++)
					{
						x = j + 0.5;
						y = (a*x + b);
						z = (c*x + d);
						YX_view->Fill(y, -x);
						YZ_view->Fill(y, z);	
					}
				}
				else
				{
					for(int j = 0; j < 5200; j++)
					{
						y = j + 0.5 - 2600.0;
						x = (y - b) / a;
						z = (c*x + d);
						if(x > 0)
						{
							YX_view->Fill(y, -x);
							YZ_view->Fill(y, z);
						}
					}	
				}
				
				if( abs(c) <= 1.0 )
				{
					for(int j = 0; j < 500; j++)
					{
						x = j + 0.5;
						z = (c*x + d);
						XZ_view->Fill(x, z);		
					}
				}
				else
				{
					for(int j = 0; j < 3000; j++)
					{
						z = j + 0.5 - 1500.0;
						x = (z - d) / c;
						if(x > 0)
						{
							XZ_view->Fill(x, z);
						}
					}	
				}
			}
			else
			{
				if( abs(a) <= 1.0 )
				{
					for(int j = 0; j < 500; j++)
					{
						x = j + 0.5 - 500.0;
						y = (a*x + b);
						z = (c*x + d);
						YX_view->Fill(y, -x);
						YZ_view->Fill(y, z);
					}
				}
				else
				{
					for(int j = 0; j < 5200; j++)
					{
						y = j + 0.5 - 2600.0;
						x = (y - b) / a;
						z = (c*x + d);
						if(x < 0)
						{
							YX_view->Fill(y, -x);
							YZ_view->Fill(y, z);
						}
					}	
				}
				
				
				if( abs(c) <= 1.0 )
				{
					for(int j = 0; j < 500; j++)
					{
						x = j + 0.5 - 500.0;
						z = (c*x + d);
						XZ_view->Fill(x, z);
					}
				}
				else
				{
					for(int j = 0; j < 3000; j++)
					{
						z = j + 0.5 - 1500.0;
						x = (z - d) / c;
						if(x < 0)
						{
							XZ_view->Fill(x, z);
						}
					}	
				}
			}	
			
			
		}

		if(ev % 1000 == 0) std::cout << ev << std::endl;
	}
	
	file_data->Close();
	
	file->WriteObject(YX_view, "YX_view");
	file->WriteObject(YZ_view, "YZ_view");
	file->WriteObject(XZ_view, "XZ_view");

	
	file->Close();
	
}






