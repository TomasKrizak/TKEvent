#include "../TKEvent/include/TKEvent.h" 

R__LOAD_LIBRARY(../TKEvent/lib/libTKEvent.so);

void show()
{
	int run_number;
	cout << "enter run number: ";
	cin >> run_number;
	
	int lower_limit, upper_limit;
	cout << "enter starting event number: ";
	cin >> lower_limit;
		
	cout << "enter final event number: (-1 for entire run)";
	cin >> upper_limit;
	
	TFile* file = new TFile(Form("../runs/Run-%d.root", run_number));
	TTree* tree = (TTree*) file->Get("Event");

	if( upper_limit == -1 ) 
	{	
		upper_limit = tree->GetEntries();
	}

	TKEvent* event;
	tree->SetBranchAddress("Eventdata", &event);

	cout << "Run number " << run_number << ", " << tree->GetEntries() << " events available." << endl << endl;
	
	//clock_t start = clock();
	int broken = 0;
	int all = 0;
	
	TKtimer timer;
	for(UInt_t i = lower_limit; i < upper_limit+1; i++)	// Loop over events
	{
		event = new TKEvent();
		if(i % 100 == 0) cout << "proccesing event number: " << i <<endl;
		tree->GetEntry(i);
		
		event->set_r("Manchester", "distance");
		//event->set_sigma_R();
		event->set_h();
		//event->set_sigma_Z();
		
		event->reconstruct(timer);
		//event->reconstruct_simple();
		//event->make_top_projection(3, 2);

		//event->print();
		//event->build_event();
		
		/*
		for(int j = 0; j < event->get_tr_hits().size(); j++)
		{
			TKtrhit* hit = event->get_tr_hits()[j]; 
			if( hit->get_SRL('S') == 1 )
			{
				if( hit->get_SRL('R') == 85 && hit->get_SRL('L') == 4 )
				{
					all++;
					if(hit->get_r() == -1)
					{
						broken++;
					}
					//event->make_top_projection(2);
					//event->print();
					break;					
				} 
			}
		}*/
		/*
		for(int j = 0; j < event->get_clusters().size(); j++)
		{

			if( event->get_cluster(j)->get_track()->get_side() == 1 )
			{
				TKtrack* track = event->get_cluster(j)->get_track();
				for(int k = 0; k < track->get_associated_tr_hits().size(); k++)
				{
					TKtrhit* hit = track->get_associated_tr_hits()[k]; 
					if( hit->get_SRL('R') == 85 && hit->get_SRL('L') == 4 )
					{
						event->make_top_projection(2);
						//event->print();
						continue;					
					} 
					
				}
			}
		}
		*/
		
		/*
		for(int j = 0; j < event->get_clusters().size(); j++)
		{
			int ambiguity;
			ambiguity = event->get_cluster(j)->get_ambiguity_type();
			if(ambiguity == 2)
			{
				if( event->get_cluster(j)->get_track()->get_side() == 1 )
				{
					if(abs(event->get_cluster(j)->get_track()->get_b() - 1252.5) < 70.0 )
					{
						//event->get_cluster(j)->get_track()->print();
						//event->get_clusters()[j]->draw_ML_vertical(run_number, i);
						//cout << ambiguity << endl;
						event->make_top_projection(2);
						event->print();
						//event->build_event();
						
					}
				}
			}
		}
			*/	
		/*
		for(int j = 0; j < event->get_no_tracks(); j++)
		{
			if( event->get_tracks()[j]->get_chi_squared() > 0.5 ) continue;
			if( event->get_tracks()[j]->get_ambiguity_type() != 0 )
			{
				ambiguous = ambiguous + 0.5;
				all = all + 0.5;
			}
			else
			{
				all++;
			}
			
			if(event->get_tracks()[j]->get_quality_Z() < 10e-10 && event->get_tracks()[j]->get_c() != 0.0)	
			{	
				cout << event->get_event_number() << endl;
				event->get_tracks()[j]->print();
				event->make_top_projection(2);
			}			
		
		}*/
		delete event;
		
	}
	//cout<< "all hits 1.85.4: " << all << "	broken hits 1.85.4: " << broken << endl;
	//cout<<"Finished after " << (clock() - start) / (double) CLOCKS_PER_SEC << " seconds." <<endl;
}

