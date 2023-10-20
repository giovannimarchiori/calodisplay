/******************************************************************************/
// Simple event display for the LAr calorimeter with inclined modules
// author: Giovanni Marchiori (giovanni.marchiori@cern.ch)
//
// class EventReader: sets up the TTreeReader branches for reading event data
//
/******************************************************************************/

#include "EventReader.h"


EventReader::EventReader(TFile* f, bool doHCal)
{
  fReader = new TTreeReader("events", f);
  nEvents = fReader->GetEntries();
  m_doHCal = doHCal;
  
  // primary particles
  genParticles_PDG = new TTreeReaderArray<Int_t>(*fReader, "genParticles.PDG");
  genParticles_generatorStatus = new TTreeReaderArray<Int_t>(*fReader, "genParticles.generatorStatus");
  genParticles_simulatorStatus = new TTreeReaderArray<Int_t>(*fReader, "genParticles.simulatorStatus");
  genParticles_charge = new TTreeReaderArray<Float_t>(*fReader, "genParticles.charge");
  genParticles_time = new TTreeReaderArray<Float_t>(*fReader, "genParticles.time");
  genParticles_mass = new TTreeReaderArray<Double_t>(*fReader, "genParticles.mass");
  genParticles_vertex_x = new TTreeReaderArray<Double_t>(*fReader, "genParticles.vertex.x");
  genParticles_vertex_y = new TTreeReaderArray<Double_t>(*fReader, "genParticles.vertex.y");
  genParticles_vertex_z = new TTreeReaderArray<Double_t>(*fReader, "genParticles.vertex.z");
  //genParticles_endpoint_x = new TTreeReaderArray<Double_t>(*fReader, "genParticles.endpoint.x");
  //genParticles_endpoint_y = new TTreeReaderArray<Double_t>(*fReader, "genParticles.endpoint.y");
  //genParticles_endpoint_z = new TTreeReaderArray<Double_t>(*fReader, "genParticles.endpoint.z");
  genParticles_momentum_x = new TTreeReaderArray<Float_t>(*fReader, "genParticles.momentum.x");
  genParticles_momentum_y = new TTreeReaderArray<Float_t>(*fReader, "genParticles.momentum.y");
  genParticles_momentum_z = new TTreeReaderArray<Float_t>(*fReader, "genParticles.momentum.z");
  
  // secondary particles
  if (fReader->GetTree()->FindBranch("SimParticleSecondaries.PDG")) {
    SimParticleSecondaries_PDG = new TTreeReaderArray<Int_t>(*fReader, "SimParticleSecondaries.PDG");
    //SimParticleSecondaries_generatorStatus = new TTreeReaderArray<Int_t>(*fReader, "SimParticleSecondaries.generatorStatus");
    //SimParticleSecondaries_simulatorStatus = new TTreeReaderArray<Int_t>(*fReader, "SimParticleSecondaries.simulatorStatus");
    //SimParticleSecondaries_charge = new TTreeReaderArray<Float_t>(*fReader, "SimParticleSecondaries.charge");
    SimParticleSecondaries_time = new TTreeReaderArray<Float_t>(*fReader, "SimParticleSecondaries.time");
    SimParticleSecondaries_mass = new TTreeReaderArray<Double_t>(*fReader, "SimParticleSecondaries.mass");
    SimParticleSecondaries_vertex_x = new TTreeReaderArray<Double_t>(*fReader, "SimParticleSecondaries.vertex.x");
    SimParticleSecondaries_vertex_y = new TTreeReaderArray<Double_t>(*fReader, "SimParticleSecondaries.vertex.y");
    SimParticleSecondaries_vertex_z = new TTreeReaderArray<Double_t>(*fReader, "SimParticleSecondaries.vertex.z");
    SimParticleSecondaries_endpoint_x = new TTreeReaderArray<Double_t>(*fReader, "SimParticleSecondaries.endpoint.x");
    SimParticleSecondaries_endpoint_y = new TTreeReaderArray<Double_t>(*fReader, "SimParticleSecondaries.endpoint.y");
    SimParticleSecondaries_endpoint_z = new TTreeReaderArray<Double_t>(*fReader, "SimParticleSecondaries.endpoint.z");
    SimParticleSecondaries_momentum_x = new TTreeReaderArray<Float_t>(*fReader, "SimParticleSecondaries.momentum.x");
    SimParticleSecondaries_momentum_y = new TTreeReaderArray<Float_t>(*fReader, "SimParticleSecondaries.momentum.y");
    SimParticleSecondaries_momentum_z = new TTreeReaderArray<Float_t>(*fReader, "SimParticleSecondaries.momentum.z");
  }
  
  // hits in ECal barrel
  ECalBarrelPositionedHits_cellID     = new TTreeReaderArray<ULong_t>(*fReader, "ECalBarrelPositionedHits.cellID");
  ECalBarrelPositionedHits_energy     = new TTreeReaderArray<Float_t>(*fReader, "ECalBarrelPositionedHits.energy");
  ECalBarrelPositionedHits_position_x = new TTreeReaderArray<Float_t>(*fReader, "ECalBarrelPositionedHits.position.x");
  ECalBarrelPositionedHits_position_y = new TTreeReaderArray<Float_t>(*fReader, "ECalBarrelPositionedHits.position.y");
  ECalBarrelPositionedHits_position_z = new TTreeReaderArray<Float_t>(*fReader, "ECalBarrelPositionedHits.position.z");
  
  // cells in ECal barrel
  ECalBarrelPositionedCells_cellID     = new TTreeReaderArray<ULong_t>(*fReader, "ECalBarrelPositionedCells.cellID");
  ECalBarrelPositionedCells_energy     = new TTreeReaderArray<Float_t>(*fReader, "ECalBarrelPositionedCells.energy");
  ECalBarrelPositionedCells_position_x = new TTreeReaderArray<Float_t>(*fReader, "ECalBarrelPositionedCells.position.x");
  ECalBarrelPositionedCells_position_y = new TTreeReaderArray<Float_t>(*fReader, "ECalBarrelPositionedCells.position.y");
  ECalBarrelPositionedCells_position_z = new TTreeReaderArray<Float_t>(*fReader, "ECalBarrelPositionedCells.position.z");
  
  // cells in ECal barrel with coarser merging
  ECalBarrelPositionedCells2_cellID     = new TTreeReaderArray<ULong_t>(*fReader, "ECalBarrelPositionedCells2.cellID");
  ECalBarrelPositionedCells2_energy     = new TTreeReaderArray<Float_t>(*fReader, "ECalBarrelPositionedCells2.energy");
  ECalBarrelPositionedCells2_position_x = new TTreeReaderArray<Float_t>(*fReader, "ECalBarrelPositionedCells2.position.x");
  ECalBarrelPositionedCells2_position_y = new TTreeReaderArray<Float_t>(*fReader, "ECalBarrelPositionedCells2.position.y");
  ECalBarrelPositionedCells2_position_z = new TTreeReaderArray<Float_t>(*fReader, "ECalBarrelPositionedCells2.position.z");

  if (m_doHCal) {
    // hits in HCal barrel
    HCalBarrelPositionedHits_cellID     = new TTreeReaderArray<ULong_t>(*fReader, "HCalBarrelPositionedHits.cellID");
    HCalBarrelPositionedHits_energy     = new TTreeReaderArray<Float_t>(*fReader, "HCalBarrelPositionedHits.energy");
    HCalBarrelPositionedHits_position_x = new TTreeReaderArray<Float_t>(*fReader, "HCalBarrelPositionedHits.position.x");
    HCalBarrelPositionedHits_position_y = new TTreeReaderArray<Float_t>(*fReader, "HCalBarrelPositionedHits.position.y");
    HCalBarrelPositionedHits_position_z = new TTreeReaderArray<Float_t>(*fReader, "HCalBarrelPositionedHits.position.z");
    
    // cells in HCal barrel
    HCalBarrelPositionedCells_cellID     = new TTreeReaderArray<ULong_t>(*fReader, "HCalBarrelPositionedCells.cellID");
    HCalBarrelPositionedCells_energy     = new TTreeReaderArray<Float_t>(*fReader, "HCalBarrelPositionedCells.energy");
    HCalBarrelPositionedCells_position_x = new TTreeReaderArray<Float_t>(*fReader, "HCalBarrelPositionedCells.position.x");
    HCalBarrelPositionedCells_position_y = new TTreeReaderArray<Float_t>(*fReader, "HCalBarrelPositionedCells.position.y");
    HCalBarrelPositionedCells_position_z = new TTreeReaderArray<Float_t>(*fReader, "HCalBarrelPositionedCells.position.z");
  }
  
  // corrected calo topo clusters
  CorrectedCaloTopoClusters_energy     = new TTreeReaderArray<Float_t>(*fReader, "CorrectedCaloTopoClusters.energy");
  CorrectedCaloTopoClusters_position_x = new TTreeReaderArray<Float_t>(*fReader, "CorrectedCaloTopoClusters.position.x");
  CorrectedCaloTopoClusters_position_y = new TTreeReaderArray<Float_t>(*fReader, "CorrectedCaloTopoClusters.position.y");
  CorrectedCaloTopoClusters_position_z = new TTreeReaderArray<Float_t>(*fReader, "CorrectedCaloTopoClusters.position.z");
  CorrectedCaloTopoClusters_hits_begin = new TTreeReaderArray<UInt_t>(*fReader, "CorrectedCaloTopoClusters.hits_begin");
  CorrectedCaloTopoClusters_hits_end   = new TTreeReaderArray<UInt_t>(*fReader, "CorrectedCaloTopoClusters.hits_end");
  
  // 1f. - cells in the topo clusters
  PositionedCaloTopoClusterCells_cellID     = new TTreeReaderArray<ULong_t>(*fReader, "PositionedCaloTopoClusterCells.cellID");
  PositionedCaloTopoClusterCells_energy     = new TTreeReaderArray<Float_t>(*fReader, "PositionedCaloTopoClusterCells.energy");
  PositionedCaloTopoClusterCells_position_x = new TTreeReaderArray<Float_t>(*fReader, "PositionedCaloTopoClusterCells.position.x");
  PositionedCaloTopoClusterCells_position_y = new TTreeReaderArray<Float_t>(*fReader, "PositionedCaloTopoClusterCells.position.y");
  PositionedCaloTopoClusterCells_position_z = new TTreeReaderArray<Float_t>(*fReader, "PositionedCaloTopoClusterCells.position.z");
  
  //qs_rhoz.clear();
  //qs_rhophi.clear();      
}

EventReader::~EventReader() {
  delete genParticles_PDG;
  delete genParticles_generatorStatus;
  delete genParticles_simulatorStatus;
  delete genParticles_charge;
  delete genParticles_time;
  delete genParticles_mass;
  delete genParticles_vertex_x;
  delete genParticles_vertex_y;
  delete genParticles_vertex_z;
  //delete genParticles_endpoint_x;
  //delete genParticles_endpoint_y;
  //delete genParticles_endpoint_z;
  delete genParticles_momentum_x;
  delete genParticles_momentum_y;
  delete genParticles_momentum_z;
  
  if (SimParticleSecondaries_PDG) {
    delete SimParticleSecondaries_PDG;
    // delete SimParticleSecondaries_generatorStatus;
    // delete SimParticleSecondaries_simulatorStatus;
    // delete SimParticleSecondaries_charge;
    delete SimParticleSecondaries_time;
    delete SimParticleSecondaries_mass;
    delete SimParticleSecondaries_vertex_x;
    delete SimParticleSecondaries_vertex_y;
    delete SimParticleSecondaries_vertex_z;
    delete SimParticleSecondaries_endpoint_x;
    delete SimParticleSecondaries_endpoint_y;
    delete SimParticleSecondaries_endpoint_z;
    delete SimParticleSecondaries_momentum_x;
    delete SimParticleSecondaries_momentum_y;
    delete SimParticleSecondaries_momentum_z;
  }
  delete ECalBarrelPositionedHits_cellID;
  delete ECalBarrelPositionedHits_energy;
  delete ECalBarrelPositionedHits_position_x;
  delete ECalBarrelPositionedHits_position_y;
  delete ECalBarrelPositionedHits_position_z;
  delete ECalBarrelPositionedCells_cellID;
  delete ECalBarrelPositionedCells_energy;
  delete ECalBarrelPositionedCells_position_x;
  delete ECalBarrelPositionedCells_position_y;
  delete ECalBarrelPositionedCells_position_z;
  delete ECalBarrelPositionedCells2_cellID;
  delete ECalBarrelPositionedCells2_energy;
  delete ECalBarrelPositionedCells2_position_x;
  delete ECalBarrelPositionedCells2_position_y;
  delete ECalBarrelPositionedCells2_position_z;
  if (m_doHCal) {
    delete ECalBarrelPositionedHits_cellID;
    delete ECalBarrelPositionedHits_energy;
    delete ECalBarrelPositionedHits_position_x;
    delete ECalBarrelPositionedHits_position_y;
    delete ECalBarrelPositionedHits_position_z;
    delete ECalBarrelPositionedCells_cellID;
    delete ECalBarrelPositionedCells_energy;
    delete ECalBarrelPositionedCells_position_x;
    delete ECalBarrelPositionedCells_position_y;
    delete ECalBarrelPositionedCells_position_z;
  }
  delete CorrectedCaloTopoClusters_energy;
  delete CorrectedCaloTopoClusters_position_x;
  delete CorrectedCaloTopoClusters_position_y;
  delete CorrectedCaloTopoClusters_position_z;
  delete CorrectedCaloTopoClusters_hits_begin;
  delete CorrectedCaloTopoClusters_hits_end;
  delete PositionedCaloTopoClusterCells_cellID;
  delete PositionedCaloTopoClusterCells_energy;
  delete PositionedCaloTopoClusterCells_position_x;
  delete PositionedCaloTopoClusterCells_position_y;
  delete PositionedCaloTopoClusterCells_position_z;
  delete fReader;
}


void EventReader::loadEvent(int event) {
  if (event<0 || event>=nEvents) {
    cout << "Event number " << event << " out of range : [0, " << nEvents-1 << "]" << endl;
  }
  fReader->SetEntry(event);
}
