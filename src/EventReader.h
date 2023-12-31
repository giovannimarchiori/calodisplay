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
  TTreeReaderArray<Int_t> *genParticles_PDG = nullptr;
  TTreeReaderArray<Int_t> *genParticles_generatorStatus = nullptr;
  TTreeReaderArray<Int_t> *genParticles_simulatorStatus = nullptr;
  TTreeReaderArray<Float_t> *genParticles_charge = nullptr;
  TTreeReaderArray<Float_t> *genParticles_time = nullptr;
  TTreeReaderArray<Double_t> *genParticles_mass = nullptr;
  TTreeReaderArray<Double_t> *genParticles_vertex_x = nullptr;
  TTreeReaderArray<Double_t> *genParticles_vertex_y = nullptr;
  TTreeReaderArray<Double_t> *genParticles_vertex_z = nullptr;
  //TTreeReaderArray<Double_t> *genParticles_endpoint_x = nullptr;
  //TTreeReaderArray<Double_t> *genParticles_endpoint_y = nullptr;
  //TTreeReaderArray<Double_t> *genParticles_endpoint_z = nullptr;
  TTreeReaderArray<Float_t> *genParticles_momentum_x = nullptr;
  TTreeReaderArray<Float_t> *genParticles_momentum_y = nullptr;
  TTreeReaderArray<Float_t> *genParticles_momentum_z = nullptr;

  // 1a. - secondary particles
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
  TTreeReaderArray<ULong_t> *ECalBarrelPositionedHits_cellID = nullptr;
  TTreeReaderArray<Float_t> *ECalBarrelPositionedHits_energy = nullptr;
  TTreeReaderArray<Float_t> *ECalBarrelPositionedHits_position_x = nullptr;
  TTreeReaderArray<Float_t> *ECalBarrelPositionedHits_position_y = nullptr;
  TTreeReaderArray<Float_t> *ECalBarrelPositionedHits_position_z = nullptr;

  // the cells in ECal barrel
  TTreeReaderArray<ULong_t> *ECalBarrelPositionedCells_cellID = nullptr;
  TTreeReaderArray<Float_t> *ECalBarrelPositionedCells_energy = nullptr;
  TTreeReaderArray<Float_t> *ECalBarrelPositionedCells_position_x = nullptr;
  TTreeReaderArray<Float_t> *ECalBarrelPositionedCells_position_y = nullptr;
  TTreeReaderArray<Float_t> *ECalBarrelPositionedCells_position_z = nullptr;

  // the cells in ECal barrel with coarser merging
  TTreeReaderArray<ULong_t> *ECalBarrelPositionedCells2_cellID = nullptr;
  TTreeReaderArray<Float_t> *ECalBarrelPositionedCells2_energy = nullptr;
  TTreeReaderArray<Float_t> *ECalBarrelPositionedCells2_position_x = nullptr;
  TTreeReaderArray<Float_t> *ECalBarrelPositionedCells2_position_y = nullptr;
  TTreeReaderArray<Float_t> *ECalBarrelPositionedCells2_position_z = nullptr;

  // the hits in HCal barrel 
  TTreeReaderArray<ULong_t> *HCalBarrelPositionedHits_cellID = nullptr;
  TTreeReaderArray<Float_t> *HCalBarrelPositionedHits_energy = nullptr;
  TTreeReaderArray<Float_t> *HCalBarrelPositionedHits_position_x = nullptr;
  TTreeReaderArray<Float_t> *HCalBarrelPositionedHits_position_y = nullptr;
  TTreeReaderArray<Float_t> *HCalBarrelPositionedHits_position_z = nullptr;

  // the cells in HCal barrel
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
  
  // the corrected calo topo clusters
  TTreeReaderArray<Float_t> *CorrectedCaloTopoClusters_energy = nullptr;
  TTreeReaderArray<Float_t> *CorrectedCaloTopoClusters_position_x = nullptr;
  TTreeReaderArray<Float_t> *CorrectedCaloTopoClusters_position_y = nullptr;
  TTreeReaderArray<Float_t> *CorrectedCaloTopoClusters_position_z = nullptr;
  TTreeReaderArray<UInt_t> *CorrectedCaloTopoClusters_hits_begin = nullptr;
  TTreeReaderArray<UInt_t> *CorrectedCaloTopoClusters_hits_end = nullptr;  
  
  // the cells in the topo clusters
  TTreeReaderArray<ULong_t> *CaloTopoClusterCells_cellID = nullptr;
  TTreeReaderArray<Float_t> *CaloTopoClusterCells_energy = nullptr;
  TTreeReaderArray<Float_t> *CaloTopoClusterCells_position_x = nullptr;
  TTreeReaderArray<Float_t> *CaloTopoClusterCells_position_y = nullptr;
  TTreeReaderArray<Float_t> *CaloTopoClusterCells_position_z = nullptr;

  // the sliding window clusters (TODO: include upstream/downstream corrections)
  TTreeReaderArray<Float_t> *CaloClusters_energy = nullptr;
  TTreeReaderArray<Float_t> *CaloClusters_position_x = nullptr;
  TTreeReaderArray<Float_t> *CaloClusters_position_y = nullptr;
  TTreeReaderArray<Float_t> *CaloClusters_position_z = nullptr;
  TTreeReaderArray<UInt_t> *CaloClusters_hits_begin = nullptr;
  TTreeReaderArray<UInt_t> *CaloClusters_hits_end = nullptr;

  // the cells in the sliding window clusters
  TTreeReaderArray<ULong_t> *PositionedCaloClusterCells_cellID = nullptr;
  TTreeReaderArray<Float_t> *PositionedCaloClusterCells_energy = nullptr;
  TTreeReaderArray<Float_t> *PositionedCaloClusterCells_position_x = nullptr;
  TTreeReaderArray<Float_t> *PositionedCaloClusterCells_position_y = nullptr;
  TTreeReaderArray<Float_t> *PositionedCaloClusterCells_position_z = nullptr;

  
 public:
  EventReader(TFile* f, bool doHCal=false, bool doSW=false, bool doTopo=true, bool drawMergedCells=false);
  ~EventReader();
  void loadEvent(int event);
  bool m_doHCal;
  bool m_doSW;
  bool m_doTopo;
  bool m_drawMergedCells;
};

#endif
