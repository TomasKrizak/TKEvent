#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

#include "TTree.h"
#include "TFile.h"

#include <snfee/snfee.h>
#include <snfee/io/multifile_data_reader.h>

#include <sncabling/om_id.h>
#include <sncabling/gg_cell_id.h>
#include <sncabling/label.h>

#include <snemo/datamodels/raw_event_data.h>
#include <snemo/datamodels/calo_digitized_hit.h>
#include <snemo/datamodels/tracker_digitized_hit.h>

#include "Event_3D.hh"

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
      std::cerr << "*** missing event_number (-e/--event EVENT_NUMBER)" << std::endl;
      return 1;
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

  std::cout << "Searching for event " << event_number << " ..." << std::endl;

//MIROOOOOOOOOOOOOOOOOO 1


TFile* subor;
TTree* strom;

Event_3D *event_3D = new Event_3D(-1, -1);
 
subor = new TFile("Default.root","NEW");
strom = new TTree("Event","All data from event");
strom->Branch("Eventdata", &event_3D);	

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

      event_3D = new Event_3D(run_number, red_event_id);

      if (event_number < red_event_id)
	break;

      event_found = true;

      // Container of merged TriggerID(s) by event builder
      const std::set<int32_t> & red_trigger_ids = red.get_origin_trigger_ids();
      
//-------------------------------------
      
//************************************************************************
//MIRO - REPLACEMENT CODE TO DIVIDE THE DATA AND THE WRITING PROCESS
	const std::vector<snemo::datamodel::calo_digitized_hit> red_calo_hits = red.get_calo_hits();
	for (const snemo::datamodel::calo_digitized_hit & red_calo_hit : red_calo_hits)
	{	
		const snemo::datamodel::timestamp & reference_time = red_calo_hit.get_reference_time();
				
		const sncabling::om_id om_id = red_calo_hit.get_om_id();
		int om_side, om_wall, om_column, om_row, om_num;
		if(om_id.is_main())
		{
			om_side   = om_id.get_side();
			om_column = om_id.get_column();
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
		
		event_3D->add_calo_hit(om_num, red_calo_hit.is_high_threshold(), reference_time.get_ticks(), red_calo_hit.get_fwmeas_peak_cell());
	}

	const std::vector<snemo::datamodel::tracker_digitized_hit> red_tracker_hits = red.get_tracker_hits();
	for (const snemo::datamodel::tracker_digitized_hit & red_tracker_hit : red_tracker_hits)
	{
	  	const sncabling::gg_cell_id gg_id = red_tracker_hit.get_cell_id();
		int cell_side  = gg_id.get_side();
		int cell_row   = gg_id.get_row();
		int cell_layer = gg_id.get_layer();

		double height = tc_sizez;
		double radius = tc_radius;
		for(const snemo::datamodel::tracker_digitized_hit::gg_times & gg_timestamps : red_tracker_hit.get_times())
		{

			const int64_t anode_tdc = gg_timestamps.get_anode_time(0).get_ticks();
			height = 0.0;
			if (anode_tdc != snfee::data::INVALID_TICKS)
			{
				const int64_t top_cathode_tdc = gg_timestamps.get_top_cathode_time().get_ticks();
				const int64_t bottom_cathode_tdc = gg_timestamps.get_bottom_cathode_time().get_ticks();
				if( top_cathode_tdc != snfee::data::INVALID_TICKS && bottom_cathode_tdc != snfee::data::INVALID_TICKS )
				{
					int64_t propagation_time_tdc = bottom_cathode_tdc + top_cathode_tdc - 2*anode_tdc;
					height = tc_sizez*( double(bottom_cathode_tdc - anode_tdc) / double(propagation_time_tdc) ) - tc_sizez/2.0;
				}
			}
			//else cout << "missing r0" << endl;

			const int64_t anode_tdc2 = gg_timestamps.get_anode_time(0).get_ticks();
			if(anode_tdc2 != snfee::data::INVALID_TICKS)
			{
				double min_time = 1e10;
				//double calo_hit = (calo_tdc * 6.25) - 400.0 + (400.0 * peak_cell / 1024.0);
				
				for(int om_hit = 0; om_hit < event_3D->get_hit_om_num_v().size(); om_hit++)
				{
					if(2*anode_tdc2 - (event_3D->get_hit_om_tdc_v()[om_hit]-44) < 800 && 2*anode_tdc2 - (event_3D->get_hit_om_tdc_v()[om_hit]-44) > 0 )
					{
						if(min_time > double(2*anode_tdc2 - (event_3D->get_hit_om_tdc_v()[om_hit]-44)))
						{
							min_time = 6.25 * (2*anode_tdc2 - (event_3D->get_hit_om_tdc_v()[om_hit]-44));
						}
					}
				}
				radius = event_3D->drift_time_to_radius(min_time, "Manu");
			}			
		}	

		// Tracker hit is added here
		event_3D->add_tracker_hit(113*9*cell_side + 9*cell_row + cell_layer, height, radius);
	} 

// END OF MIRO
//************************************************************************



      //event_3D->add_event(red);
	
    /*  event_3D->print();

      event_3D->reconstruct_track();

      event_3D->print_tracks();*/

//////////////////////////////////
	
	strom->Fill();	

//-------------------------------------
    
      //break;

    } // (while red_source.has_record_tag())

	strom->Write();
	subor->Close();	
	//delete subor;

  if (!event_found)
    std::cerr << "=> Event was not found ! (only " << red_counter <<  " RED in this file)" << std::endl;
    
  snfee::terminate();

  return 0;
}

