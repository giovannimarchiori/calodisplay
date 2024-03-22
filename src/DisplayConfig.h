/******************************************************************************/
// Simple event display for the LAr calorimeter with inclined modules
// author: Giovanni Marchiori (giovanni.marchiori@cern.ch)
//
// class DisplayConfig: holds the configuration of the event display
//
/******************************************************************************/

#ifndef DISPLAYCONFIG_H
#define DISPLAYCONFIG_H

#include <string>
#include <map>

class DisplayConfig {
public:
  std::map<std::string, std::string> configMap;

  /*
  bool drawGenParticles;
  bool drawSimParticles;
  bool drawECalBarrelHits;
  bool drawECalBarrelCells;
  bool drawECalBarrelMergedCells;
  bool drawHCalBarrelHits;
  bool drawHCalBarrelCells;
  bool drawTopoClusters;
  bool drawCaloClusters;
  
  std::string genParticles;
  std::string simParticles;
  std::string ecalBarrelHits;
  std::string ecalBarrelCells;
  std::string ecalBarrelMergedCells;
  std::string hcalBarrelHits;
  std::string hcalBarrelCells;
  std::string topoClusters;
  std::string topoClusterCells;
  std::string caloClusters;
  std::string caloClusterCells;
  */
  
  DisplayConfig();
  void Print();
  void ReadFromFile(std::string filename);
  void WriteToFile(std::string filename);

  void setStringConfig(const std::string& key, const std::string& value) {
    configMap[key] = value;
  }

  void setBoolConfig(const std::string& key, bool value) {
    configMap[key] = value ? "true" : "false";
  }

  void setFloatConfig(const std::string& key, float value) {
    configMap[key] = std::to_string(value);
  }

  std::string getStringConfig(const std::string& key) const {
    if (configMap.find(key) != configMap.end()) {
      return configMap.at(key);
    }
    return ""; // or throw exception for missing config
  }

  bool getBoolConfig(const std::string& key) const {
    if (configMap.find(key) != configMap.end()) {
      return configMap.at(key) == "true";
    }
    return false; // or throw exception for missing config
  }

  float getFloatConfig(const std::string& key) const {
    if (configMap.find(key) != configMap.end()) {
      return std::stof(configMap.at(key));
    }
    return -9999999.; // or throw exception for missing config
  }

};

#endif
