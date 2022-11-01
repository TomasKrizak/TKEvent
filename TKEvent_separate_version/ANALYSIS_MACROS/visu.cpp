void visu()
{
	cout << endl << "warning: used dimensions and positions of demonstrator parts may not be fully correct" << endl << endl;

	int run_number;
	int event_number;
	
	cout << "Choose run (-1 to show legend): ";
	cin  >> run_number; 
	
	if(run_number == -1)
	{
		cout << endl << "3D model contains:" << endl;
		cout << "	Calorimeter hits: low treshold hits (yellow)" << endl;
		cout << "	                  high treshold hits (red)" << endl;
		cout << "	Tracker hits: hits with usable radii (red)" << endl;
		cout << "	              hits with missing or wrong radii (yellow)" << endl;
		cout << "	              hits with missing height are set to z = 0.0" << endl;
		cout << "	Bi calibration sources: light blue (positions not clear at the moment)" << endl;
		cout << "	Reconstrcution (if calculated): one line per side of detector (blue)" << endl << endl;
	}
	else
	{
		cout << "Choose event to visualize (-1 to exit): ";
		cin  >> event_number; 
		
		if(event_number != -1)
		{
			TCanvas* c = new TCanvas("3D_demonstrator", "", 1920, 1080);
			
			TLatex* title = new TLatex(-0.9, 0.9, Form("Run %d | Event %d", run_number, event_number));
			title->SetTextSize(0.03);
			title->Draw();
			
			TGLViewer* view = (TGLViewer*)gPad->GetViewer3D();
			TGeoManager* manager = new TGeoManager();
			
			TFile* file = new TFile(Form("./Run-%d_event-%d_3D.root", run_number, event_number));
			TGeoVolume* geo;
			file->GetObject("demonstrator", geo);
			
			manager->SetTopVolume(geo);
			manager->CloseGeometry(); 
			geo->Draw("same");
			
			TPolyLine3D *track1;
			file->GetObject(Form("track-%d", 0), track1);
			track1->Draw("same");
			
			TPolyLine3D *track2;
			file->GetObject(Form("track-%d", 1), track2);
			track2->Draw("same");
		}
	}
}


