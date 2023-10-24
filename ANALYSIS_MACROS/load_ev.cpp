#include "../TKEvent/include/TKEvent.h" 

R__LOAD_LIBRARY(../TKEvent/lib/libTKEvent.so);

void load_ev()
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

	TKEvent* event = new TKEvent();
	tree->SetBranchAddress("Eventdata", &event);

	cout << "Run number " << run_number << ", " << tree->GetEntries() << " events available." << endl << endl;

	for(UInt_t i = lower_limit; i < upper_limit + 1; i++)	// Loop over events
	{
		tree->GetEntry(i);
		
		
		event->set_r("Manchester", "distance");	// needed for reconstruction
		
		event->reconstruct_ML(0);

		event->print();
		
		event->make_top_projection(2);
		
		event->build_event();
	}
}

