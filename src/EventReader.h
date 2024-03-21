/******************************************************************************/
// Simple event display for the LAr calorimeter with inclined modules
// author: Giovanni Marchiori (giovanni.marchiori@cern.ch)
//
// class EventReader: sets up the TTreeReader branches for reading event data
//
/******************************************************************************/

#ifndef EVENTREADER_H
#define EVENTREADER_H

/******************************************************************************/
// dependencies
/******************************************************************************/
#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
using namespace std;

/******************************************************************************/
// main class
/******************************************************************************/

class EventReader
{
public:
  Int_t nEvents = 0; // Number of events in file
  Int_t eventId = 0; // Current event id

  //private:
  TTreeReader* fReader = nullptr;

  // 1a. - primary particles
  bool m_readGenParticles = true;
  std::string m_genParticlesBranch = "genParticles";
  TTreeReaderArray<Int_t> *genParticles_PDG = nullptr;
  TTreeReaderArray<Int_t> *genParticles_generatorStatus = nullptr;
  TTreeReaderArray<Int_t> *genParticles_simulatorStatus = nullptr;
  TTreeReaderArray<Float_t> *genParticles_charge = nullptr;
  TTreeReaderArray<Float_t> *genParticles_time = nullptr;
  TTreeReaderArray<Double_t> *genParticles_mass = nullptr;
  TTreeReaderArray<Double_t> *genParticles_vertex_x = nullptr;
  TTreeReaderArray<Double_t> *genParticles_vertex_y = nullptr;
  TTreeReaderArray<Double_t> *genParticles_vertex_z = nullptr;
  // TTreeReaderArray<Double_t> *genParticles_endpoint_x = nullptr;
  // TTreeReaderArray<Double_t> *genParticles_endpoint_y = nullptr;
  // TTreeReaderArray<Double_t> *genParticles_endpoint_z = nullptr;
  // these variables are double in latest versions
  // TTreeReaderArray<Float_t> *genParticles_momentum_x = nullptr;
  // TTreeReaderArray<Float_t> *genParticles_momentum_y = nullptr;
  // TTreeReaderArray<Float_t> *genParticles_momentum_z = nullptr;
  TTreeReaderArray<Double_t> *genParticles_momentum_x = nullptr;
  TTreeReaderArray<Double_t> *genParticles_momentum_y = nullptr;
  TTreeReaderArray<Double_t> *genParticles_momentum_z = nullptr;

  // 1a. - secondary particles
  bool m_readSimParticles = true;
  std::string m_simParticlesBranch = "SimParticlesSecondaries";
  TTreeReaderArray<Int_t> *SimParticleSecondaries_PDG = nullptr;
  //TTreeReaderArray<Int_t> *SimParticleSecondaries_generatorStatus = nullptr;
  //TTreeReaderArray<Int_t> *SimParticleSecondaries_simulatorStatus = nullptr;
  //TTreeReaderArray<Float_t> *SimParticleSecondaries_charge = nullptr;
  TTreeReaderArray<Float_t> *SimParticleSecondaries_time = nullptr;
  TTreeReaderArray<Double_t> *SimParticleSecondaries_mass = nullptr;
  TTreeReaderArray<Double_t> *SimParticleSecondaries_vertex_x = nullptr;
  TTreeReaderArray<Double_t> *SimParticleSecondaries_vertex_y = nullptr;
  TTreeReaderArray<Double_t> *SimParticleSecondaries_vertex_z = nullptr;
  TTreeReaderArray<Double_t> *SimParticleSecondaries_endpoint_x = nullptr;
  TTreeReaderArray<Double_t> *SimParticleSecondaries_endpoint_y = nullptr;
  TTreeReaderArray<Double_t> *SimParticleSecondaries_endpoint_z = nullptr;
  TTreeReaderArray<Float_t> *SimParticleSecondaries_momentum_x = nullptr;
  TTreeReaderArray<Float_t> *SimParticleSecondaries_momentum_y = nullptr;
  TTreeReaderArray<Float_t> *SimParticleSecondaries_momentum_z = nullptr;

  // the hits in ECal barrel
  bool m_readECalBarrelHits = true;
  std::string m_ecalBarrelHitsBranch = "ECalBarrelPositionedHits";
  TTreeReaderArray<ULong_t> *ECalBarrelPositionedHits_cellID = nullptr;
  TTreeReaderArray<Float_t> *ECalBarrelPositionedHits_energy = nullptr;
  TTreeReaderArray<Float_t> *ECalBarrelPositionedHits_position_x = nullptr;
  TTreeReaderArray<Float_t> *ECalBarrelPositionedHits_position_y = nullptr;
  TTreeReaderArray<Float_t> *ECalBarrelPositionedHits_position_z = nullptr;

  // the cells in ECal barrel
  bool m_readECalBarrelCells = true;
  std::string m_ecalBarrelCellsBranch = "ECalBarrelPositionedCells";
  TTreeReaderArray<ULong_t> *ECalBarrelPositionedCells_cellID = nullptr;
  TTreeReaderArray<Float_t> *ECalBarrelPositionedCells_energy = nullptr;
  TTreeReaderArray<Float_t> *ECalBarrelPositionedCells_position_x = nullptr;
  TTreeReaderArray<Float_t> *ECalBarrelPositionedCells_position_y = nullptr;
  TTreeReaderArray<Float_t> *ECalBarrelPositionedCells_position_z = nullptr;

  // the cells in ECal barrel with coarser merging
  bool m_readECalBarrelCells2 = false;
  std::string m_ecalBarrelCells2Branch = "ECalParrelPositionedCells2";
  TTreeReaderArray<ULong_t> *ECalBarrelPositionedCells2_cellID = nullptr;
  TTreeReaderArray<Float_t> *ECalBarrelPositionedCells2_energy = nullptr;
  TTreeReaderArray<Float_t> *ECalBarrelPositionedCells2_position_x = nullptr;
  TTreeReaderArray<Float_t> *ECalBarrelPositionedCells2_position_y = nullptr;
  TTreeReaderArray<Float_t> *ECalBarrelPositionedCells2_position_z = nullptr;

  // the hits in HCal barrel
  bool m_readHCalBarrelHits = true;
  std::string m_hcalBarrelHitsBranch = "HCalBarrelPositionedHits";
  TTreeReaderArray<ULong_t> *HCalBarrelPositionedHits_cellID = nullptr;
  TTreeReaderArray<Float_t> *HCalBarrelPositionedHits_energy = nullptr;
  TTreeReaderArray<Float_t> *HCalBarrelPositionedHits_position_x = nullptr;
  TTreeReaderArray<Float_t> *HCalBarrelPositionedHits_position_y = nullptr;
  TTreeReaderArray<Float_t> *HCalBarrelPositionedHits_position_z = nullptr;

  // the cells in HCal barrel
  bool m_readHCalBarrelCells = true;
  std::string m_hcalBarrelCellsBranch = "HCalBarrelPositionedCells";
  TTreeReaderArray<ULong_t> *HCalBarrelPositionedCells_cellID = nullptr;
  TTreeReaderArray<Float_t> *HCalBarrelPositionedCells_energy = nullptr;
  TTreeReaderArray<Float_t> *HCalBarrelPositionedCells_position_x = nullptr;
  TTreeReaderArray<Float_t> *HCalBarrelPositionedCells_position_y = nullptr;
  TTreeReaderArray<Float_t> *HCalBarrelPositionedCells_position_z = nullptr;

  // IF WE ARE SURE THAT WE DO NOT PLOT BOTH SW AND TOPOCLUSTERS AT THE SAME TIME
  // we could simplify the code of the EventDisplay class by having a single
  // cluster reader, that points to either SW or topo clusters based on the clusters
  // that we want to draw
  // For the moment keep both (we might want to have SW clusters in ECAL for
  // EM objects and topoclusters for ECAL+HCAL for jets)
  
  // the (corrected) calo topo clusters
  bool m_readTopoClusters = false;
  std::string m_topoClustersBranch = "CaloTopoClusters";
  TTreeReaderArray<Float_t> *CaloTopoClusters_energy = nullptr;
  TTreeReaderArray<Float_t> *CaloTopoClusters_position_x = nullptr;
  TTreeReaderArray<Float_t> *CaloTopoClusters_position_y = nullptr;
  TTreeReaderArray<Float_t> *CaloTopoClusters_position_z = nullptr;
  TTreeReaderArray<UInt_t> *CaloTopoClusters_hits_begin = nullptr;
  TTreeReaderArray<UInt_t> *CaloTopoClusters_hits_end = nullptr;  
  
  // the (positioned) cells in the topo clusters
  bool m_readTopoClusterCells = false;
  std::string m_topoClusterCellsBranch = "CaloTopoClusterCells";
  TTreeReaderArray<ULong_t> *CaloTopoClusterCells_cellID = nullptr;
  TTreeReaderArray<Float_t> *CaloTopoClusterCells_energy = nullptr;
  TTreeReaderArray<Float_t> *CaloTopoClusterCells_position_x = nullptr;
  TTreeReaderArray<Float_t> *CaloTopoClusterCells_position_y = nullptr;
  TTreeReaderArray<Float_t> *CaloTopoClusterCells_position_z = nullptr;

  // the (corrected) sliding window clusters
  bool m_readCaloClusters = false;
  std::string m_caloClustersBranch = "CaloClusters";
  TTreeReaderArray<Float_t> *CaloClusters_energy = nullptr;
  TTreeReaderArray<Float_t> *CaloClusters_position_x = nullptr;
  TTreeReaderArray<Float_t> *CaloClusters_position_y = nullptr;
  TTreeReaderArray<Float_t> *CaloClusters_position_z = nullptr;
  TTreeReaderArray<UInt_t> *CaloClusters_hits_begin = nullptr;
  TTreeReaderArray<UInt_t> *CaloClusters_hits_end = nullptr;

  // the (positioned) cells in the sliding window clusters
  bool m_readCaloClusterCells = false;
  std::string m_caloClusterCellsBranch = "CaloClusterCells";
  TTreeReaderArray<ULong_t> *CaloClusterCells_cellID = nullptr;
  TTreeReaderArray<Float_t> *CaloClusterCells_energy = nullptr;
  TTreeReaderArray<Float_t> *CaloClusterCells_position_x = nullptr;
  TTreeReaderArray<Float_t> *CaloClusterCells_position_y = nullptr;
  TTreeReaderArray<Float_t> *CaloClusterCells_position_z = nullptr;

  
 public:
  EventReader(TFile* f, bool doHCal=false, bool doSW=false, bool doTopo=true, bool drawMergedCells=false);
  ~EventReader();
  void SetBranches();
  void loadEvent(int event);
  bool m_doHCal;
  bool m_doSW;
  bool m_doTopo;
  bool m_readMergedCells;
};

#endif
