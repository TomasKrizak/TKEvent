// Standard headers
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

// ROOT headers
#include "TTree.h"
#include "TFile.h"

// SNFEE headers
#include <snfee/snfee.h>
#include <snfee/io/multifile_data_reader.h>

// SNCabling headers
#include <sncabling/om_id.h>
#include <sncabling/gg_cell_id.h>
#include <sncabling/label.h>

// Bayeux headers
#include <snemo/datamodels/raw_event_data.h>
#include <snemo/datamodels/calo_digitized_hit.h>
#include <snemo/datamodels/tracker_digitized_hit.h>

// TK headers
#include "TKEvent.h"

int main (int argc, char *argv[])
{
	const char *red_path = getenv("RED_PATH");

	int run_number = -1;
	int event_number = -1;

	std::string input_filename = "";

	for (int iarg=1; iarg<argc; ++iarg)
	{
      		std::string arg (argv[iarg]);
		if (arg[0] == '-')
		{
			if (arg=="-i" || arg=="--input")
			input_filename = std::string(argv[++iarg]);

			else if (arg=="-r" || arg=="--run")
			run_number = atoi(argv[++iarg]);

			else if (arg=="-e" || arg=="--event")
			event_number = atoi(argv[++iarg]);
		}
	}

	if (event_number == -1)
	{
		event_number = INT_MAX;
	}

	if (input_filename.empty())
	{
		if (run_number == -1)
		{
			std::cerr << "*** missing run_number (-r/--run RUN_NUMBER)" << std::endl;
			return 1;
		}
      
		char input_filename_buffer[128];
		snprintf(input_filename_buffer, sizeof(input_filename_buffer),
		       "%s/snemo_run-%d_red.data.gz", red_path, run_number);
		input_filename = std::string(input_filename_buffer);
	}

	snfee::initialize();

	/// Configuration for raw data reader
	snfee::io::multifile_data_reader::config_type reader_cfg;
	reader_cfg.filenames.push_back(input_filename);

	// Instantiate a reader
	std::cout << "Opening " << input_filename << " ..." << std::endl;
	snfee::io::multifile_data_reader red_source (reader_cfg);

	// Working RED object
	snemo::datamodel::raw_event_data red;

	// RED counter
	std::size_t red_counter = 0;
	bool event_found = false;

//MIROOOOOOOOOOOOOOOOOO 1

	TFile* file;
	TTree* tree;

	TKEvent *event = new TKEvent();
	 
	file = new TFile(Form("Run-%d.root", run_number),"RECREATE");
	tree = new TTree("Event","All data from event");
	tree->Branch("Eventdata", &event);	

//MIROOOOOOOOOOOOOOOOOO 1 END

	while (red_source.has_record_tag())
	{
		// Check the serialization tag of the next record:
		DT_THROW_IF(!red_source.record_tag_is(snemo::datamodel::raw_event_data::SERIAL_TAG),
			  std::logic_error, "Unexpected record tag '" << red_source.get_record_tag() << "'!");

		// Load the next RED object:
		red_source.load(red);
		red_counter++;

		// Event number
		int32_t red_event_id = red.get_event_id();

		event = new TKEvent(run_number, red_event_id);

		if (event_number < red_event_id)
		break;

		event_found = true;

		// Container of merged TriggerID(s) by event builder
		const std::set<int32_t> & red_trigger_ids = red.get_origin_trigger_ids();
	     
	//-------------------------------------

		const std::vector<snemo::datamodel::calo_digitized_hit> red_calo_hits = red.get_calo_hits();
		for (const snemo::datamodel::calo_digitized_hit & red_calo_hit : red_calo_hits)
		{	
			bool is_HT;
			if(red_calo_hit.is_high_threshold())
			{
				is_HT = true;
			}
			else if(red_calo_hit.is_low_threshold_only())
			{
				is_HT = false;
			}
			else continue;
		
			const snemo::datamodel::timestamp & reference_time = red_calo_hit.get_reference_time();
				
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
		
			event->add_OM_hit(om_num, is_HT, reference_time.get_ticks(), red_calo_hit.get_fwmeas_peak_cell());
		}

		const std::vector<snemo::datamodel::tracker_digitized_hit> red_tracker_hits = red.get_tracker_hits();
		for (const snemo::datamodel::tracker_digitized_hit & red_tracker_hit : red_tracker_hits)
		{
		  	const sncabling::gg_cell_id gg_id = red_tracker_hit.get_cell_id();
			int srl[3] = {gg_id.get_side(), gg_id.get_row(), gg_id.get_layer()};
			for(const snemo::datamodel::tracker_digitized_hit::gg_times & gg_timestamps : red_tracker_hit.get_times())
			{
				
				int64_t tsp[7];
		
				for (int it = 0; it < 5; it++)
				{
					if(gg_timestamps.get_anode_time(it).get_ticks() != snfee::data::INVALID_TICKS)
					{
						tsp[it] = gg_timestamps.get_anode_time(it).get_ticks();
					}
					else
					{
						tsp[it] = -1;
					}
				}

				if(gg_timestamps.get_bottom_cathode_time().get_ticks() != snfee::data::INVALID_TICKS)
				{
					tsp[5] = gg_timestamps.get_bottom_cathode_time().get_ticks();
				}
				else
				{
					tsp[5] = -1;
				}

				if(gg_timestamps.get_top_cathode_time().get_ticks() != snfee::data::INVALID_TICKS)
				{
					tsp[6] = gg_timestamps.get_top_cathode_time().get_ticks();
				}
				else
				{
					tsp[6] = -1;
				}
		
				// Tracker hit is added here
				event->add_tracker_hit(srl, tsp);			
				
			}	
		} 

		//event->set_r("Manchester");
		//event->reconstruct_track(0);
		
		//event->print();
		//event->print_tracks();
		//event->make_top_projection();
		//event->build_event();
		
		tree->Fill();	
		
		if(red_event_id % 10000 == 0)
		{
			std::cout << "Event No. " << red_event_id << " converted!" << std::endl;
		}
		
	} // (while red_source.has_record_tag())

	tree->Write();
	file->Close();	
	//delete file;

	snfee::terminate();

	return 0;
}

