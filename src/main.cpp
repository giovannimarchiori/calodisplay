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

#include <csignal>
#include <filesystem>
#include <sstream>
#include <string>
#include <vector>
#include <TRint.h>

using namespace std;

std::string configFile = "config.json";

std::string usage(bool verbose = false)
{
  EventDisplay aDisplay;
  std::string use = "Usage:\n\n";
  use += "calodisplay [-h] [--nog4] [--fulldet] [--dohcal] [--doendcap] [--nocache] [--debug] [-g|--geometry <geomFile.root>] [-e|--events <eventFile.root>] [-c|--config <configFile.json>] [--set <key>=<value>] [-s|--skipEvents <nevents>]\n";
  use += "            [--screenshots <prefix> [--nevents <n>] [--views <views>] [--width <pixels>]\n";
  use += "            [--focus <cluster|track|track+cluster|all>[:<scale>]] [--camera <view>=<values>]]\n\n";
  use += "-h|--help           : print this help, with the default values of the options, and exit\n";
  use += "--nog4              : use a simplified geometry rather than the G4 geometry\n";
  use += "--fulldet           : show the full detector (otherwise only the ECAL is shown)\n";
  use += "--dohcal            : show the HCAL volumes and readout (not necessary if --fulldet is set)\n";
  use += "--doendcap          : show the ECAL and (with --dohcal) HCAL endcaps (not necessary if --fulldet is set)\n";
  use += "--nocache           : read the full G4 geometry, not the one without the volumes that are never drawn\n";
  use += "                      from the geometry cache (slower, but all the volumes can be turned on in the browser)\n";
  use += "-g|--geometry <geomFile.root> : the G4 geometry of the detector (searched also in the installation\n";
  use += "                      directory)\n";
  use += "-e|--events <eventFile.root> : the events to display (without it, only the detector is shown)\n";
  use += "-c|--config <configFile.json> : the configuration of the display (collections to draw, energy\n";
  use += "                      thresholds, ...; searched also in the installation directory)\n";
  use += "--set <key>=<value> : set a parameter of the configuration (overriding the configuration file),\n";
  use += "                      e.g. --set drawCaloClusters=true (can be repeated)\n";
  use += "-s|--skipEvents <nevents> : start from event <nevents> (the first one is 0)\n";
  use += "--debug             : display extended information on the console during execution\n";
  use += "--screenshots <prefix> : save the views of the events to <prefix>_ev<event>_<view>.png and exit\n";
  use += "                      (still needs an X display, e.g. a VNC session)\n";
  use += "--nevents <n>       : number of events for the screenshots, starting from the first one (-s) (default 1)\n";
  use += "--views <views>     : comma-separated views for the screenshots: 3d, rhophi, rhoz (default all)\n";
  use += "--width <pixels>    : width of the screenshots (default 1920)\n";
  use += "--focus <cluster|track|track+cluster|all>[:<scale>] : frame the most energetic cluster, charged particle,\n";
  use += "                      or both, in the projections of the screenshots, with their box scaled by <scale> filling\n";
  use += "                      the image (default 2 for clusters, 1.2 otherwise), or the full detector (default); the 3D\n";
  use += "                      view is not changed\n";
  use += "--camera <view>=<values> : set the camera of a view in the screenshots (overriding --focus):\n";
  use += "                      rhophi=<x>,<y>,<zoom>  rhoz=<z>,<rho>,<zoom>\n";
  use += "                      3d=<fov>,<center x,y,z>,<eye x,y,z>,<up x,y,z>\n";
  use += "                      (the Camera button of the display prints them for the current view; can be repeated)\n\n";

  if (verbose)
  {
    use += "Default values for options with parameters:\n\n";
    use += "  geomFile: ";
    use += aDisplay.geomFile;
    use += "\n";
    use += "  eventFile: ";
    use += aDisplay.evtFile;
    use += "\n";
    use += "  configFile: ";
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

  // screenshots mode
  std::string screenshotPrefix;
  int screenshotEvents = 1;
  std::vector<std::string> screenshotViews = {"3d", "rhophi", "rhoz"};
  int screenshotWidth = 1920;

  // to find the files distributed with calodisplay (icons, field map, ..)
  PathUtils::setInstallDir(argv[0]);

  // parse the command line
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

      else if (arg == "--screenshots" || arg == "--nevents" || arg == "--views" || arg == "--width" ||
               arg == "--focus" || arg == "--camera")
      {
        if (i == args.size() - 1)
          throw std::runtime_error("missing value after " + std::string(arg) + " option!");
        i++;
        const std::string value(args[i]);
        if (arg == "--screenshots")
          screenshotPrefix = value;
        else if (arg == "--nevents")
          screenshotEvents = std::stoi(value);
        else if (arg == "--width")
          screenshotWidth = std::stoi(value);
        else if (arg == "--focus")
        {
          const size_t colon = value.find(':');
          display.cameraFocus = value.substr(0, colon);
          if (colon != std::string::npos)
            display.cameraFocusScale = std::stod(value.substr(colon + 1));
          if (display.cameraFocus != "cluster" && display.cameraFocus != "track" &&
              display.cameraFocus != "track+cluster" && display.cameraFocus != "all")
            throw std::runtime_error("unknown focus " + display.cameraFocus + " (cluster, track, track+cluster or all)");
        }
        else if (arg == "--camera")
        {
          const size_t eq = value.find('=');
          if (eq == std::string::npos)
            throw std::runtime_error("--camera needs <view>=<values>");
          std::vector<double> values;
          std::stringstream ss(value.substr(eq + 1));
          std::string number;
          while (std::getline(ss, number, ','))
            values.push_back(std::stod(number));
          display.cameraSettings[value.substr(0, eq)] = values;
        }
        else
        {
          screenshotViews.clear();
          std::stringstream ss(value);
          std::string view;
          while (std::getline(ss, view, ','))
            screenshotViews.push_back(view);
        }
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
  if (screenshotPrefix != "")
  {
    // (not a TRint, which would quit at the end of the standard input, e.g.
    // when run from a script)
    TApplication *app = new TApplication("display", &argc, argv);
    display.startDisplay(skipEvents);
    display.saveScreenshots(screenshotPrefix, skipEvents, screenshotEvents, screenshotViews, screenshotWidth);
    app->Terminate(0);
  }
  TRint *app = new TRint("display", &argc, argv);
  // ignore the resizes of the terminal while the display starts: the handler
  // of the ROOT prompt is not async-signal-safe (it calls setenv), and if it
  // interrupts an allocation, e.g. during the long projection of the
  // geometry, the application hangs
  auto resizeHandler = std::signal(SIGWINCH, SIG_IGN);
  display.startDisplay(skipEvents);
  std::signal(SIGWINCH, resizeHandler);
  app->Run();
  return 0;
}
