#include "../TKEvent/include/TKEvent.h" 

R__LOAD_LIBRARY(../TKEvent/lib/libTKEvent.so);

void load_ev()
{
	TFile* f = new TFile("./Run-728.root");
	TTree* s = (TTree*) f->Get("Event");

	TKEvent* Eve = new TKEvent(-1,-1);
	s->SetBranchAddress("Eventdata", &Eve);

	for(UInt_t i=0; i < s->GetEntries(); i++)	// Loop over events
	{
		s->GetEntry(i);
		Eve->print();
		Eve->reconstruct_track();
		Eve->make_top_projection();
      		Eve->build_event();

	}	
}
