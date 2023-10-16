/******************************************************************************/
// Simple event display for the LAr calorimeter with inclined modules
// author: Giovanni Marchiori (giovanni.marchiori@cern.ch)
//
// Compile with make and execute with ./display -h to see available options
/******************************************************************************/

#include "EventDisplay.h"

#include <filesystem>
#include <string>
#include <TRint.h>

using namespace std;

const std::string usage = "usage: display [--nog4] [-g|--geometry <geomFile>] [-e|--events <eventFile>] [-s|--skipEvents <nevents>]\n";

int main(int argc, char* argv[]) {

  EventDisplay display;
  int skipEvents=0;
  
  const std::vector<std::string_view> args(argv + 1, argv + argc);

  // check if -h option, in that case print usage and exit
  for (auto arg:args) {
    if (arg == "-h" || arg == "--help") {
      std::cout << usage << std::endl;
      exit(0);
    }
  }

  // parse the arguments
  try {
    for (unsigned int i=0; i<args.size(); i++) {
      std::string_view arg = args[i];

      if (arg == "--nog4") {
	display.useG4geom = false;
	continue;
      }      
      else if (arg == "-g" || arg == "--geometry") {
	if (i==args.size()-1) {
	  throw std::runtime_error("missing geometry file after -g/--geometry option!");
	}
	else {
	  i++;
	  display.geomFile = args[i];
	}
	continue;
      }
      
      else if (arg == "-e" || arg == "--events") {
	if (i==args.size()-1) {
	  throw std::runtime_error("missing event file after -e/--events option!");
	}
	else {
	  i++;
	  display.evtFile = args[i];
	}
	continue;
      }
      
      else if (arg == "-s" || arg == "--skipEvents") {
	if (i==args.size()-1) {
	  throw std::runtime_error("missing number of events to skip after -s/--skipEvents option!");
	}
	else {
	  i++;
	  skipEvents = stoi(string(args[i]));
	  if (skipEvents<0) skipEvents = 0;
	}
	continue;
      }
      
      else {
	std::cerr << "Unknown option " << arg << std::endl;
	std::cerr << usage << std::endl;
	return EXIT_FAILURE;
      }
    }
  } catch (const std::exception &x) {
    std::cerr << "Exception: " << x.what() << '\n';
    std::cerr << usage << std::endl;
    return EXIT_FAILURE;
  }
  
  // check if the files exist
  if (!std::filesystem::exists(display.geomFile)) {
    std::cerr << "Geometry file " << display.geomFile << " does not exist!!!"  << std::endl;
    return EXIT_FAILURE;
  }
  if (!std::filesystem::exists(display.evtFile)) {
    std::cerr << "Event file " << display.evtFile << " does not exist!!!"  << std::endl;
    return EXIT_FAILURE;
  }
  
  argc = 0;
  argv = nullptr;
  TRint* gApplication = new TRint("display", &argc, argv);
  display.startDisplay(skipEvents);
  gApplication->Run();
  return 0;
}
