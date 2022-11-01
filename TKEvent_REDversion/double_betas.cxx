#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

#include <snfee/snfee.h>
#include <snfee/io/multifile_data_reader.h>

#include <sncabling/om_id.h>
#include <sncabling/gg_cell_id.h>
#include <sncabling/label.h>

#include <snemo/datamodels/raw_event_data.h>
#include <snemo/datamodels/calo_digitized_hit.h>
#include <snemo/datamodels/tracker_digitized_hit.h>

#include "build_event_3D.cc"
#include <TH2D.h>

int main (int argc, char *argv[])
{
  const char *red_path = getenv("RED_PATH");

  int run_number = -1;

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
	}
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
      event_found = true;

      // Container of merged TriggerID(s) by event builder
      const std::set<int32_t> & red_trigger_ids = red.get_origin_trigger_ids();
      
//-------------------------------------
      
      int event_number = red_event_id;
      Event_3D *event_3D = new Event_3D(run_number, event_number);
	event_3D->add_event(red);
	event_3D->reconstruct_track();
	
	if(event_3D->get_no_tracks() != 2) continue;
	
	double x,y,z;
	double a0 = event_3D->get_track_a(0);
	double b0 = event_3D->get_track_b(0);
	double c0 = event_3D->get_track_c(0);
	double d0 = event_3D->get_track_d(0);
	
	double a1 = event_3D->get_track_a(1);
	double b1 = event_3D->get_track_b(1);
	double c1 = event_3D->get_track_c(1);
	double d1 = event_3D->get_track_d(1);
	
	if(c0 == 0.0 && d0 == 0.0) continue;
	if(c1 == 0.0 && d1 == 0.0) continue;
	
	if(abs(b0 - b1) > 20.0) continue;
	if(abs(d0 - d1) > 30.0) continue;
	
	event_3D->make_top_projection();
	//event_3D->build_event_3D();
	std::cout << event_number << std::endl;


	if(event_number % 1000 == 0) std::cout << event_number << std::endl;
	if(event_number > 100000) break;

//-------------------------------------

    } // (while red_source.has_record_tag())
	
	

  snfee::terminate();

  return 0;
}

