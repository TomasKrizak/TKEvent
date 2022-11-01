#include "../Event3D/include/Event_3D.hh" 

R__LOAD_LIBRARY(../Event3D/lib/libEvent3D.so);

void load_ev()
{
	TFile* f = new TFile("./Default.root");
	TTree* s = (TTree*) f->Get("Event");

	Event_3D* Eve = new Event_3D(-1,-1);
	s->SetBranchAddress("Eventdata", &Eve);

	for(UInt_t i=0; i < s->GetEntries(); i++)	// Loop over events
	{
		s->GetEntry(i);
		Eve->print();
		//Eve->make_top_projection();
      		//Eve->build_event_3D();

	}	
}
