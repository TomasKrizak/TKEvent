#include "../TKEvent/include/TKEvent.h" 

R__LOAD_LIBRARY(../TKEvent/lib/libTKEvent.so);

void reco_ev()
{
	int run_number;
	cout << "enter run number: ";
	cin >> run_number;
	
	int lower_limit, upper_limit;
	cout << "enter starting event number: ";
	cin >> lower_limit;
		
	cout << "enter final event number: (-1 for entire run)";
	cin >> upper_limit;
	
   	TFile* reco_file = new TFile(Form("../runs/Run-%d_reco.root", run_number), "RECREATE");
	TTree* trec = new TTree("Event","All data from event");
	TFile* file = new TFile(Form("../runs/Run-%d.root", run_number));
	TTree* tree = (TTree*) file->Get("Event");
    
	if( upper_limit == -1 ) 
	{	
		upper_limit = tree->GetEntries();
	}
 
	TKEvent* event; 
   	TKEvent* reco_event = new TKEvent();
	tree->SetBranchAddress("Eventdata", &event);
	trec->Branch("Eventdata", &reco_event);	

	cout << "Run number " << run_number << ", " << tree->GetEntries() << " events available." << endl << endl;
 
	for(UInt_t i = lower_limit; i < upper_limit + 1; i++)	// Loop over events
	{
	 
		event = new TKEvent();
		tree->GetEntry(i);
		
		event->set_r("Manchester", "distance");	// calculates tracker hit radii
		//event->set_sigma_R(); 	// space for implementation of better R uncertainty model
		
		event->set_h();		// calculates tracker hit heights
		//event->set_sigma_Z(); 	// space for implementation of better Z uncertainty model
		
		event->reconstruct(0); 	// currently best reconstuction method
        
		//reco_event = new TKEvent();
		reco_event = event;
		trec->Fill();
			
		delete event;
        
       	if(i % 100 == 0)
		{
			std::cout << "Event No. " << i << " reconstructed!" << std::endl;
		}
	}
    
	reco_file->cd();
	trec->Write();

	reco_file->Close();
	file->Close();
}

