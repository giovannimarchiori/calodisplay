/******************************************************************************/
// Simple event display for the ALLEGRO detector with ECAL with inclined modules
// author: Giovanni Marchiori (giovanni.marchiori@cern.ch)
//
// main.cpp: source for executable
// 
// Compile with make and execute with ./calodisplay -h to see available options
//
/******************************************************************************/

#include "EventDisplay.h"
#include "Globals.h"

#include <filesystem>
#include <string>
#include <TRint.h>

using namespace std;

std::string configFile = "config_ddsim.json";

std::string usage(bool verbose = false)
{
  EventDisplay aDisplay;
  std::string use = "Usage:\n";
  use += "calodisplay [--nog4] [--fulldet] [--dohcal] [--doendcap] [--sw] [--notopo] [--debug] [-g|--geometry <geomFile.root>] [-e|--events <eventFile.root>] [-c|--config <configFile.json>] [-s|--skipEvents <nevents>]\n\n";
  use += "--nog4      : use simplified geometry rather than G4 geometry in file <geomFile.root>\n";
  use += "--fulldet   : show the full detector (including tracker, muon detector, interaction region)\n";
  use += "--dohcal    : show HCAL volume, hits and cells\n";
  use += "--doendcap  : show ECAL and (optionally) HCAL endcaps\n";
  use += "--sw        : show sliding-window clusters\n";
  use += "--notopo    : do not show topoclusters\n\n";
  use += "--debug     : display extended information on the console during execution\n\n";

  if (verbose)
  {
    use += "Default values for options with parameters:\n";
    use += "geomFile: ";
    use += aDisplay.geomFile;
    use += "\n";
    use += "eventFile: ";
    use += aDisplay.evtFile;
    use += "\n";
    use += "configFile: ";
    use += configFile;
    use += "\n";
    use += "\n";
  }
  return use;
}

int main(int argc, char *argv[])
{

  EventDisplay display;
  int skipEvents = 0;

  const std::vector<std::string_view> args(argv + 1, argv + argc);

  // check if -h option, in that case print usage and exit
  for (auto arg : args)
  {
    if (arg == "-h" || arg == "--help")
    {
      std::cout << usage(true) << std::endl;
      exit(0);
    }
  }

  // check if -c or --config options, to load default configuration
  for (unsigned int i = 0; i < args.size(); i++)
  {
    std::string_view arg = args[i];
    if (arg == "-c" || arg == "--config")
    {
      if (i == args.size() - 1)
      {
	throw std::runtime_error("missing config file after -c/--config option!");
      }
      else
      {
	configFile = args[i+1];
      }
      break;
    }
  }

  // check that config file exists
  if (configFile != "" && !std::filesystem::exists(configFile))
  {
    std::cerr << "Config file " << configFile << " does not exist!!!" << std::endl;
    return EXIT_FAILURE;
  }

  // load default configuration
  displayConfig.ReadFromFile(configFile);
  
 
  // parse the other arguments
  try
  {
    for (unsigned int i = 0; i < args.size(); i++)
    {
      std::string_view arg = args[i];

      if (arg == "--nog4")
      {
        display.useG4geom = false;
        continue;
      }

      else if (arg == "--dohcal")
      {
        display.doHCal = true;
        continue;
      }

      else if (arg == "--doendcap")
      {
        display.doEndcaps = true;
        continue;
      }

      else if (arg == "--fulldet")
      {
        display.showFullDetector = true;
        continue;
      }

      else if (arg == "--sw")
      {
        displayConfig.setBoolConfig("drawCaloClusters", true);
        continue;
      }

      else if (arg == "--notopo")
      {
        displayConfig.setBoolConfig("drawTopoClusters", false);
        continue;
      }

      else if (arg == "--debug")
      {
        display.debug = true;
        continue;
      }

      else if (arg == "-g" || arg == "--geometry")
      {
        if (i == args.size() - 1)
        {
          throw std::runtime_error("missing geometry file after -g/--geometry option!");
        }
        else
        {
          i++;
          display.geomFile = args[i];
        }
        continue;
      }

      else if (arg == "-e" || arg == "--events")
      {
        if (i == args.size() - 1)
        {
          throw std::runtime_error("missing event file after -e/--events option!");
        }
        else
        {
          i++;
          display.evtFile = args[i];
        }
        continue;
      }

      else if (arg == "-c" || arg == "--config")
      {
	/*
        if (i == args.size() - 1)
        {
          throw std::runtime_error("missing config file after -c/--config option!");
        }
        else
        {
          i++;
          configFile = args[i];
        }
        continue;
	*/
	i++; continue;
      }

      else if (arg == "-s" || arg == "--skipEvents")
      {
        if (i == args.size() - 1)
        {
          throw std::runtime_error("missing number of events to skip after -s/--skipEvents option!");
        }
        else
        {
          i++;
          skipEvents = stoi(string(args[i]));
          if (skipEvents < 0)
            skipEvents = 0;
        }
        continue;
      }

      else
      {
        std::cerr << "Unknown option " << arg << std::endl;
        std::cerr << usage() << std::endl;
        return EXIT_FAILURE;
      }
    }
  }
  catch (const std::exception &x)
  {
    std::cerr << "Exception: " << x.what() << '\n';
    std::cerr << usage() << std::endl;
    return EXIT_FAILURE;
  }

  // check if the files exist
  if (display.useG4geom && !std::filesystem::exists(display.geomFile))
  {
    std::cerr << "G4 geometry file " << display.geomFile << " does not exist!!!" << std::endl;
    return EXIT_FAILURE;
  }
  if (display.evtFile != "" && !std::filesystem::exists(display.evtFile))
  {
    std::cerr << "Event file " << display.evtFile << " does not exist!!!" << std::endl;
    return EXIT_FAILURE;
  }

  // check if there are some conflicting options
  if (!display.useG4geom && display.showFullDetector)
  {
    std::cerr << "Conflicting options --nog4 and --fulldet selected: simplified geometry not implemented for full detector!!!" << std::endl;
    return EXIT_FAILURE;
  }

  // print the configuration
  displayConfig.Print();

  // start the application
  argc = 0;
  argv = nullptr;
  TRint *gApplication = new TRint("display", &argc, argv);
  display.startDisplay(skipEvents);
  gApplication->Run();
  return 0;
}
