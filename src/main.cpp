/******************************************************************************/
// Event display for the ALLEGRO detector for FCC-ee
// author: Giovanni Marchiori (giovanni.marchiori@cern.ch)
//
// main.cpp: source for executable
// 
// Compile with make and execute with ./calodisplay -h to see available options
//
/******************************************************************************/

#include "EventDisplay.h"
#include "Globals.h"
#include "PathUtils.h"

#include <filesystem>
#include <string>
#include <TRint.h>

using namespace std;

std::string configFile = "config.json";

std::string usage(bool verbose = false)
{
  EventDisplay aDisplay;
  std::string use = "Usage:\n";
  use += "calodisplay [--nog4] [--fulldet] [--dohcal] [--doendcap] [--nocache] [--debug] [-g|--geometry <geomFile.root>] [-e|--events <eventFile.root>] [-c|--config <configFile.json>] [--set <key>=<value>] [-s|--skipEvents <nevents>]\n\n";
  use += "--nog4      : use simplified geometry rather than G4 geometry in file <geomFile.root>\n";
  use += "--fulldet   : show the full detector (including tracker, muon detector, interaction region),\n";
  use += "              with HCAL and endcaps\n";
  use += "--dohcal    : show the HCAL volumes and readout\n";
  use += "--doendcap  : show the ECAL and (with --dohcal) HCAL endcaps\n";
  use += "--nocache   : read the full G4 geometry, not the one without the volumes that are never drawn\n";
  use += "              from the geometry cache (slower, but all the volumes can be turned on in the browser)\n";
  use += "--debug     : display extended information on the console during execution\n";
  use += "--set <key>=<value> : set a parameter of the configuration (overriding the configuration file),\n";
  use += "              e.g. --set drawCaloClusters=true (can be repeated)\n\n";

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

  // to find the files distributed with calodisplay (icons, field map, ..)
  PathUtils::setInstallDir(argv[0]);

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

  // check that config file exists (also in the installation directory)
  configFile = PathUtils::findFile(configFile);
  if (configFile != "" && !std::filesystem::exists(configFile))
  {
    std::cerr << "Config file " << configFile << " does not exist!!!" << std::endl;
    return EXIT_FAILURE;
  }

  // load default configuration
  if (configFile != "" && !displayConfig.readFromFile(configFile))
  {
    return EXIT_FAILURE;
  }
  
 
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

      else if (arg == "--set")
      {
        if (i == args.size() - 1)
          throw std::runtime_error("missing <key>=<value> after --set option!");
        i++;
        std::string setting(args[i]);
        size_t eq = setting.find('=');
        if (eq == std::string::npos || eq == 0)
          throw std::runtime_error("expected <key>=<value> after --set option, got " + setting);
        std::string key = setting.substr(0, eq);
        if (displayConfig.configMap.count(key) == 0)
          std::cout << "WARNING: --set: " << key << " is not in the configuration, adding it" << std::endl;
        displayConfig.setStringConfig(key, setting.substr(eq + 1));
        continue;
      }

      else if (arg == "--nocache")
      {
        display.useGeometryCache = false;
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
          std::string s(args[i]);
          displayConfig.setStringConfig("eventFile", s);
        }
        continue;
      }

      else if (arg == "-c" || arg == "--config")
      {
        i++;
        continue;
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

  // check if the files exist (the geometry also in the installation directory)
  display.geomFile = PathUtils::findFile(display.geomFile);
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

  // the full detector includes the HCAL and the endcaps
  if (display.showFullDetector)
  {
    display.doHCal = true;
    display.doEndcaps = true;
  }

  // check if there are some conflicting options
  if (!display.useG4geom && display.showFullDetector)
  {
    std::cerr << "Conflicting options --nog4 and --fulldet selected: simplified geometry not implemented for full detector!!!" << std::endl;
    return EXIT_FAILURE;
  }

  // store event file and geom file in config
  displayConfig.setStringConfig("geometryFile", display.geomFile);
  displayConfig.setStringConfig("eventFile", display.evtFile);
  displayConfig.setStringConfig("configFile", configFile);
  
  // print the configuration
  displayConfig.print();

  // start the application
  argc = 0;
  argv = nullptr;
  TRint *app = new TRint("display", &argc, argv);
  display.startDisplay(skipEvents);
  app->Run();
  return 0;
}
