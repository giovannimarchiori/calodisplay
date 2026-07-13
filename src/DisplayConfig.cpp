/******************************************************************************/
// Simple event display for the ALLEGRO detector with ECAL with inclined modules
// author: Giovanni Marchiori (giovanni.marchiori@cern.ch)
//
// class DisplayConfig: holds the configuration of the event display
//
/******************************************************************************/

#include "DisplayConfig.h"
#include <iostream>
#include <fstream>
#include <TBufferJSON.h>
#include <TString.h>
#include <sstream>

DisplayConfig::DisplayConfig() {
  setBoolConfig("drawGenParticles", true);
  setBoolConfig("drawSimParticles", true);
  setBoolConfig("drawVertexHits", true);
  setBoolConfig("drawMainTrackerHits", true);
  setBoolConfig("drawECalBarrelHits", true);
  setBoolConfig("drawECalBarrelCells", true);
  setBoolConfig("drawECalBarrelMergedCells", false);
  setBoolConfig("drawHCalBarrelHits", true);
  setBoolConfig("drawHCalBarrelCells", true);
  setBoolConfig("drawTopoClusters", true);
  setBoolConfig("drawCaloClusters", false);

  setStringConfig("bField", "GlobalSolenoid");
  setStringConfig("genParticles", "genParticles");
  setStringConfig("simParticles", "SimParticlesSecondaries");
  setStringConfig("vertexInnerBarrelHits", "VertexInnerBarrelCollection");
  setStringConfig("vertexOuterBarrelHits", "VertexOuterBarrelCollection");
  setStringConfig("vertexEndcapHits", "VertexEndcapCollection");
  setStringConfig("mainTrackerHits", "MainTracker_simHits");
  setStringConfig("ecalBarrelHits", "ECalBarrelPositionedHits");
  setStringConfig("ecalBarrelCells", "ECalBarrelPositionedCells");
  setStringConfig("ecalBarrelMergedCells", "ECalParrelPositionedCells2");
  setStringConfig("hcalBarrelHits", "HCalBarrelPositionedHits");
  setStringConfig("hcalBarrelCells", "HCalBarrelPositionedCells");
  setStringConfig("topoClusters", "CaloTopoClusters");
  setStringConfig("topoClusterCells", "CaloTopoClusterCells");
  setStringConfig("caloClusters", "CaloClusters");
  setStringConfig("caloClusterCells", "CaloClusterCells");

  // energy thresholds in GeV
  setFloatConfig("energyThresholdParticles", 1.0);
  setFloatConfig("energyThresholdHits", 0.0);
  setFloatConfig("energyThresholdCells", 0.0);
  setFloatConfig("energyThresholdClusters", 1.0);

  // default bkg color
  setIntConfig("bkgColor3D", 1);

  // files
  setStringConfig("configFile", "");
  setStringConfig("geometryFile", "");
  setStringConfig("eventFile", "");
}


void DisplayConfig::print()
{
  std::cout << std::endl;  
  std::cout << "******************************************************************************" << std::endl;
  std::cout << " Display options: " << std::endl;
  std::cout << "******************************************************************************" << std::endl;
  std::cout << std::endl;

  for (auto& it: configMap) {
    std::cout << it.first << " : " << it.second << std::endl;
  }
  std::cout << std::endl;
}

void DisplayConfig::writeToFile(std::string filename) {
  // to stream DisplayConfig when it uses bool/string:
  // TString json = TBufferJSON::ToJSON(this);
  // to stream the configMap when a map is used:
  TString json = TBufferJSON::ToJSON(&configMap, TBufferJSON::kMapAsObject);
  json = json.ReplaceAll("{", "{\n  ");
  json = json.ReplaceAll("\",", "\",\n ");
  json = json.ReplaceAll("}", "\n}");
  std::ofstream out(filename);
  out << json.Data() << std::endl;
  out.close();
}

void DisplayConfig::readFromFile(std::string filename) {
  std::ifstream in(filename);
  std::stringstream buffer;
  buffer << in.rdbuf();
  std::string json = buffer.str();
  std::unique_ptr<std::map<std::string, std::string>> b = TBufferJSON::FromJSON<std::map<std::string, std::string>>(json);
  // this one keeps only the values in the JSON
  // configMap = *b;
  // this one merges the values in the JSON with the defaults if not specified
  for (const auto& [key, value] : *b) {
    configMap[key] = value;
  }
  in.close();
}
