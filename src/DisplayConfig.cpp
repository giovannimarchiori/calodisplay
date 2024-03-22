/******************************************************************************/
// Simple event display for the LAr calorimeter with inclined modules
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
  setBoolConfig("drawECalBarrelHits", true);
  setBoolConfig("drawECalBarrelCells", true);
  setBoolConfig("drawECalBarrelMergedCells", false);
  setBoolConfig("drawHCalBarrelHits", true);
  setBoolConfig("drawHCalBarrelCells", true);
  setBoolConfig("drawTopoClusters", true);
  setBoolConfig("drawCaloClusters", false);

  setStringConfig("genParticles", "genParticles");
  setStringConfig("simParticles", "SimParticlesSecondaries");
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
}


void DisplayConfig::Print()
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
  
  /*
  std::cout << "drawGenParticles : " << drawGenParticles << std::endl;
  std::cout << "drawSimParticles : " << drawSimParticles << std::endl;
  std::cout << "drawECalBarrelHits : " << drawECalBarrelHits << std::endl;
  std::cout << "drawECalBarrelCells : " << drawECalBarrelCells << std::endl;
  std::cout << "drawECalBarrelMergedCells : " << drawECalBarrelMergedCells << std::endl;
  std::cout << "drawHCalBarrelHits : " << drawHCalBarrelHits << std::endl;
  std::cout << "drawHCalBarrelCells : " << drawHCalBarrelCells << std::endl;
  std::cout << "drawTopoClusters : " << drawTopoClusters << std::endl;
  std::cout << "drawCaloClusters : " << drawCaloClusters << std::endl;

  std::cout << "genParticles : " << genParticles << std::endl;
  std::cout << "simParticles : " << simParticles << std::endl;
  std::cout << "ecalBarrelHits : " << ecalBarrelHits << std::endl;
  std::cout << "ecalBarrelCells : " << ecalBarrelCells << std::endl;
  std::cout << "ecalBarrelMergedCells : " << ecalBarrelMergedCells << std::endl;
  std::cout << "hcalBarrelHits : " << hcalBarrelHits << std::endl;
  std::cout << "hcalBarrelCells : " << hcalBarrelCells << std::endl;
  std::cout << "topoClusters : " << topoClusters << std::endl;
  std::cout << "topoClusterCells : " << topoClusterCells << std::endl;
  std::cout << "caloClusters : " << caloClusters << std::endl;
  std::cout << "caloClusterCells : " << caloClusterCells << std::endl;
  */
}

void DisplayConfig::WriteToFile(std::string filename) {
  // to stream DisplayConfig when it uses bool/string
  // TString json = TBufferJSON::ToJSON(this);
  // to stream the configMap when a map is used
  TString json = TBufferJSON::ToJSON(&configMap, TBufferJSON::kMapAsObject);
  json = json.ReplaceAll("{", "{\n  ");
  json = json.ReplaceAll("\",", "\",\n ");
  json = json.ReplaceAll("}", "\n}");
  std::ofstream out(filename);
  out << json.Data() << std::endl;
  out.close();
}

void DisplayConfig::ReadFromFile(std::string filename) {
  std::ifstream in(filename);
  std::stringstream buffer;
  buffer << in.rdbuf();
  std::string json = buffer.str();
  // std::unique_ptr<DisplayConfig> b = TBufferJSON::FromJSON<DisplayConfig>(json);
  // *this = *b;
  std::unique_ptr<std::map<std::string, std::string>> b = TBufferJSON::FromJSON<std::map<std::string, std::string>>(json);
  configMap = *b;
  in.close();
}
