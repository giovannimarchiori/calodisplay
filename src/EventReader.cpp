/******************************************************************************/
// Simple event display for the LAr calorimeter with inclined modules
// author: Giovanni Marchiori (giovanni.marchiori@cern.ch)
//
// class EventReader: sets up the TTreeReader branches for reading event data
//
/******************************************************************************/

#include "EventReader.h"


EventReader::EventReader(TFile* f, bool doHCal, bool doSW, bool doTopo, bool drawMergedCells)
{
  fReader = new TTreeReader("events", f);
  nEvents = fReader->GetEntries();
  m_doHCal = doHCal;
  m_doSW = doSW;
  m_doTopo = doTopo;
  m_readMergedCells = drawMergedCells;

}

void EventReader::SetBranches()
{
  // primary particles
  if (m_genParticlesBranch=="") m_readGenParticles = false;
  if (m_readGenParticles)
  {
    if (! fReader->GetTree()->FindBranch(Form("%s.PDG", m_genParticlesBranch.c_str())))
    {
      std::cout << "WARNING: branch " << m_genParticlesBranch << ".PDG not found, disabling gen particles" << std::endl;
      m_readGenParticles = false;
      m_genParticlesBranch = "";
    }
    else
    {
      genParticles_PDG = new TTreeReaderArray<Int_t>(*fReader, Form("%s.PDG", m_genParticlesBranch.c_str()));
      genParticles_generatorStatus = new TTreeReaderArray<Int_t>(*fReader, Form("%s.generatorStatus", m_genParticlesBranch.c_str()));
      genParticles_simulatorStatus = new TTreeReaderArray<Int_t>(*fReader, Form("%s.simulatorStatus", m_genParticlesBranch.c_str()));
      genParticles_charge = new TTreeReaderArray<Float_t>(*fReader, Form("%s.charge", m_genParticlesBranch.c_str()));
      genParticles_time = new TTreeReaderArray<Float_t>(*fReader, Form("%s.time", m_genParticlesBranch.c_str()));
      genParticles_mass = new TTreeReaderArray<Double_t>(*fReader, Form("%s.mass", m_genParticlesBranch.c_str()));
      genParticles_vertex_x = new TTreeReaderArray<Double_t>(*fReader, Form("%s.vertex.x", m_genParticlesBranch.c_str()));
      genParticles_vertex_y = new TTreeReaderArray<Double_t>(*fReader, Form("%s.vertex.y", m_genParticlesBranch.c_str()));
      genParticles_vertex_z = new TTreeReaderArray<Double_t>(*fReader, Form("%s.vertex.z", m_genParticlesBranch.c_str()));
      //genParticles_endpoint_x = new TTreeReaderArray<Double_t>(*fReader, Form("%s.endpoint.x", m_genParticlesBranch.c_str()));
      //genParticles_endpoint_y = new TTreeReaderArray<Double_t>(*fReader, Form("%s.endpoint.y", m_genParticlesBranch.c_str()));
      //genParticles_endpoint_z = new TTreeReaderArray<Double_t>(*fReader, Form("%s.endpoint.z", m_genParticlesBranch.c_str()));
      // these variables are double in the latest version of the EDM
      // genParticles_momentum_x = new TTreeReaderArray<Float_t>(*fReader, Form("%s.momentum.x", m_genParticlesBranch.c_str()));
      // genParticles_momentum_y = new TTreeReaderArray<Float_t>(*fReader, Form("%s.momentum.y", m_genParticlesBranch.c_str()));
      // genParticles_momentum_z = new TTreeReaderArray<Float_t>(*fReader, Form("%s.momentum.z", m_genParticlesBranch.c_str()));
      genParticles_momentum_x = new TTreeReaderArray<Double_t>(*fReader, Form("%s.momentum.x", m_genParticlesBranch.c_str()));
      genParticles_momentum_y = new TTreeReaderArray<Double_t>(*fReader, Form("%s.momentum.y", m_genParticlesBranch.c_str()));
      genParticles_momentum_z = new TTreeReaderArray<Double_t>(*fReader, Form("%s.momentum.z", m_genParticlesBranch.c_str()));
    }
  }
  // secondary particles
  if (m_simParticlesBranch=="") m_readSimParticles = false;
  if (m_readSimParticles) {
    if (! fReader->GetTree()->FindBranch(Form("%s.PDG",m_simParticlesBranch.c_str()))) {
      std::cout << "WARNING: branch " << m_simParticlesBranch << ".PDG not found, disabling sim particles" << std::endl;
      m_readSimParticles = false;
      m_simParticlesBranch = "";
    }
    else
    {
      SimParticleSecondaries_PDG = new TTreeReaderArray<Int_t>(*fReader, Form("%s.PDG", m_simParticlesBranch.c_str()));
      //SimParticleSecondaries_generatorStatus = new TTreeReaderArray<Int_t>(*fReader, Form("%s.generatorStatus", m_simParticlesBranch.c_str()));
      //SimParticleSecondaries_simulatorStatus = new TTreeReaderArray<Int_t>(*fReader, Form("%s.simulatorStatus", m_simParticlesBranch.c_str()));
      //SimParticleSecondaries_charge = new TTreeReaderArray<Float_t>(*fReader, Form("%s.charge", m_simParticlesBranch.c_str()));
      SimParticleSecondaries_time = new TTreeReaderArray<Float_t>(*fReader, Form("%s.time", m_simParticlesBranch.c_str()));
      SimParticleSecondaries_mass = new TTreeReaderArray<Double_t>(*fReader, Form("%s.mass", m_simParticlesBranch.c_str()));
      SimParticleSecondaries_vertex_x = new TTreeReaderArray<Double_t>(*fReader, Form("%s.vertex.x", m_simParticlesBranch.c_str()));
      SimParticleSecondaries_vertex_y = new TTreeReaderArray<Double_t>(*fReader, Form("%s.vertex.y", m_simParticlesBranch.c_str()));
      SimParticleSecondaries_vertex_z = new TTreeReaderArray<Double_t>(*fReader, Form("%s.vertex.z", m_simParticlesBranch.c_str()));
      SimParticleSecondaries_endpoint_x = new TTreeReaderArray<Double_t>(*fReader, Form("%s.endpoint.x", m_simParticlesBranch.c_str()));
      SimParticleSecondaries_endpoint_y = new TTreeReaderArray<Double_t>(*fReader, Form("%s.endpoint.y", m_simParticlesBranch.c_str()));
      SimParticleSecondaries_endpoint_z = new TTreeReaderArray<Double_t>(*fReader, Form("%s.endpoint.z", m_simParticlesBranch.c_str()));
      SimParticleSecondaries_momentum_x = new TTreeReaderArray<Float_t>(*fReader, Form("%s.momentum.x", m_simParticlesBranch.c_str()));
      SimParticleSecondaries_momentum_y = new TTreeReaderArray<Float_t>(*fReader, Form("%s.momentum.y", m_simParticlesBranch.c_str()));
      SimParticleSecondaries_momentum_z = new TTreeReaderArray<Float_t>(*fReader, Form("%s.momentum.z", m_simParticlesBranch.c_str()));
    }
  }
  
  // hits in ECal barrel
  if (m_ecalBarrelHitsBranch=="") m_readECalBarrelHits = false;
  if (m_readECalBarrelHits)
  {
    if (! fReader->GetTree()->FindBranch(Form("%s.cellID", m_ecalBarrelHitsBranch.c_str())))
    {
      std::cout << "WARNING: branch " << m_ecalBarrelHitsBranch << ".cellID not found, disabling ecal barrel hits" << std::endl;
      m_readECalBarrelHits = false;
      m_ecalBarrelHitsBranch = "";
    }
    else
    {
      ECalBarrelPositionedHits_cellID     = new TTreeReaderArray<ULong_t>(*fReader, Form("%s.cellID", m_ecalBarrelHitsBranch.c_str()));
      ECalBarrelPositionedHits_energy     = new TTreeReaderArray<Float_t>(*fReader, Form("%s.energy", m_ecalBarrelHitsBranch.c_str()));
      ECalBarrelPositionedHits_position_x = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.x", m_ecalBarrelHitsBranch.c_str()));
      ECalBarrelPositionedHits_position_y = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.y", m_ecalBarrelHitsBranch.c_str()));
      ECalBarrelPositionedHits_position_z = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.z", m_ecalBarrelHitsBranch.c_str()));
    }
  }
  
  // cells in ECal barrel
  if (m_ecalBarrelCellsBranch=="") m_readECalBarrelCells = false;
  if (m_readECalBarrelCells)
  {
    if (! fReader->GetTree()->FindBranch(Form("%s.cellID", m_ecalBarrelCellsBranch.c_str())))
    {
      std::cout << "WARNING: branch " << m_ecalBarrelCellsBranch << ".cellID not found, disabling ecal barrel cells" << std::endl;
      m_readECalBarrelCells = false;
      m_ecalBarrelCellsBranch = "";
    }
    else
    {
      ECalBarrelPositionedCells_cellID     = new TTreeReaderArray<ULong_t>(*fReader, Form("%s.cellID", m_ecalBarrelCellsBranch.c_str()));
      ECalBarrelPositionedCells_energy     = new TTreeReaderArray<Float_t>(*fReader, Form("%s.energy", m_ecalBarrelCellsBranch.c_str()));
      ECalBarrelPositionedCells_position_x = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.x", m_ecalBarrelCellsBranch.c_str()));
      ECalBarrelPositionedCells_position_y = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.y", m_ecalBarrelCellsBranch.c_str()));
      ECalBarrelPositionedCells_position_z = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.z", m_ecalBarrelCellsBranch.c_str()));
    }
  }
  
  // cells in ECal barrel with coarser merging
  if (m_readMergedCells) {
    ECalBarrelPositionedCells2_cellID     = new TTreeReaderArray<ULong_t>(*fReader, "ECalBarrelPositionedCells2.cellID");
    ECalBarrelPositionedCells2_energy     = new TTreeReaderArray<Float_t>(*fReader, "ECalBarrelPositionedCells2.energy");
    ECalBarrelPositionedCells2_position_x = new TTreeReaderArray<Float_t>(*fReader, "ECalBarrelPositionedCells2.position.x");
    ECalBarrelPositionedCells2_position_y = new TTreeReaderArray<Float_t>(*fReader, "ECalBarrelPositionedCells2.position.y");
    ECalBarrelPositionedCells2_position_z = new TTreeReaderArray<Float_t>(*fReader, "ECalBarrelPositionedCells2.position.z");
  }

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

  if (m_doSW) {
    // SW clusters
    if (m_caloClustersBranch=="") m_readCaloClusters = false;
    if (m_readCaloClusters) {
      if (! fReader->GetTree()->FindBranch(Form("%s.energy", m_caloClustersBranch.c_str()))) {
	std::cout << "WARNING: branch " << m_caloClustersBranch << ".energy not found, disabling calo clusters" << std::endl;
	m_readCaloClusters= false;
	m_caloClustersBranch = "";
      }
      else {
	CaloClusters_energy     = new TTreeReaderArray<Float_t>(*fReader, Form("%s.energy", m_caloClustersBranch.c_str()));
	CaloClusters_position_x = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.x", m_caloClustersBranch.c_str()));
	CaloClusters_position_y = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.y", m_caloClustersBranch.c_str()));
	CaloClusters_position_z = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.z", m_caloClustersBranch.c_str()));
	CaloClusters_hits_begin = new TTreeReaderArray<UInt_t> (*fReader, Form("%s.hits_begin", m_caloClustersBranch.c_str()));
	CaloClusters_hits_end   = new TTreeReaderArray<UInt_t> (*fReader, Form("%s.hits_end", m_caloClustersBranch.c_str()));
      }
    }
    if (m_caloClusterCellsBranch=="") m_readCaloClusterCells = false;
    if (m_readCaloClusterCells) {
      if (! fReader->GetTree()->FindBranch(Form("%s.cellID", m_caloClusterCellsBranch.c_str()))) {
	std::cout << "WARNING: branch " << m_caloClusterCellsBranch << ".cellID not found, disabling calo cluster cells" << std::endl;
	m_readCaloClusterCells= false;
	m_caloClusterCellsBranch = "";
      }
      else {
	// cells in the SW clusters
	CaloClusterCells_cellID     = new TTreeReaderArray<ULong_t>(*fReader, Form("%s.cellID", m_caloClusterCellsBranch.c_str()));
	CaloClusterCells_energy     = new TTreeReaderArray<Float_t>(*fReader, Form("%s.energy", m_caloClusterCellsBranch.c_str()));
	CaloClusterCells_position_x = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.x", m_caloClusterCellsBranch.c_str()));
	CaloClusterCells_position_y = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.y", m_caloClusterCellsBranch.c_str()));
	CaloClusterCells_position_z = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.z", m_caloClusterCellsBranch.c_str()));
      }
    }
  }

  if (m_doTopo) {
    // topo clusters
    if (m_topoClustersBranch=="") m_readTopoClusters = false;
    if (m_readTopoClusters) {
      if (! fReader->GetTree()->FindBranch(Form("%s.energy", m_topoClustersBranch.c_str()))) {
	std::cout << "WARNING: branch " << m_topoClustersBranch << ".energy not found, disabling topo clusters" << std::endl;
	m_readTopoClusters= false;
	m_topoClustersBranch = "";
      }
      else {
	CaloTopoClusters_energy     = new TTreeReaderArray<Float_t>(*fReader, Form("%s.energy", m_topoClustersBranch.c_str()));
	CaloTopoClusters_position_x = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.x", m_topoClustersBranch.c_str()));
	CaloTopoClusters_position_y = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.y", m_topoClustersBranch.c_str()));
	CaloTopoClusters_position_z = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.z", m_topoClustersBranch.c_str()));
	CaloTopoClusters_hits_begin = new TTreeReaderArray<UInt_t> (*fReader, Form("%s.hits_begin", m_topoClustersBranch.c_str()));
	CaloTopoClusters_hits_end   = new TTreeReaderArray<UInt_t> (*fReader, Form("%s.hits_end", m_topoClustersBranch.c_str()));
      }
    }
    
    // cells in the topo clusters
    if (m_topoClusterCellsBranch=="") m_readTopoClusterCells = false;
    if (m_readTopoClusterCells) {
      if (! fReader->GetTree()->FindBranch(Form("%s.cellID", m_topoClusterCellsBranch.c_str()))) {
	std::cout << "WARNING: branch " << m_topoClusterCellsBranch << ".cellID not found, disabling topo cluster cells" << std::endl;
	m_readTopoClusterCells= false;
	m_topoClusterCellsBranch = "";
      }
      else {
	CaloTopoClusterCells_cellID     = new TTreeReaderArray<ULong_t>(*fReader, Form("%s.cellID", m_topoClusterCellsBranch.c_str()));
	CaloTopoClusterCells_energy     = new TTreeReaderArray<Float_t>(*fReader, Form("%s.energy", m_topoClusterCellsBranch.c_str()));
	CaloTopoClusterCells_position_x = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.x", m_topoClusterCellsBranch.c_str()));
	CaloTopoClusterCells_position_y = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.y", m_topoClusterCellsBranch.c_str()));
	CaloTopoClusterCells_position_z = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.z", m_topoClusterCellsBranch.c_str()));
      }
    }
  }
}

EventReader::~EventReader() {
  if (m_readGenParticles) {
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
  }
  
  if (m_readSimParticles) {
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

  if (m_readECalBarrelHits) {
    delete ECalBarrelPositionedHits_cellID;
    delete ECalBarrelPositionedHits_energy;
    delete ECalBarrelPositionedHits_position_x;
    delete ECalBarrelPositionedHits_position_y;
    delete ECalBarrelPositionedHits_position_z;
  }

  if (m_readECalBarrelCells) {
    delete ECalBarrelPositionedCells_cellID;
    delete ECalBarrelPositionedCells_energy;
    delete ECalBarrelPositionedCells_position_x;
    delete ECalBarrelPositionedCells_position_y;
    delete ECalBarrelPositionedCells_position_z;
  }

  if (m_readMergedCells) {
    delete ECalBarrelPositionedCells2_cellID;
    delete ECalBarrelPositionedCells2_energy;
    delete ECalBarrelPositionedCells2_position_x;
    delete ECalBarrelPositionedCells2_position_y;
    delete ECalBarrelPositionedCells2_position_z;
  }
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
  if (m_doSW) {
    if (m_readCaloClusters) {
      delete CaloClusters_energy;
      delete CaloClusters_position_x;
      delete CaloClusters_position_y;
      delete CaloClusters_position_z;
      delete CaloClusters_hits_begin;
      delete CaloClusters_hits_end;
    }
    if (m_readCaloClusterCells) {
      delete CaloClusterCells_cellID;
      delete CaloClusterCells_energy;
      delete CaloClusterCells_position_x;
      delete CaloClusterCells_position_y;
      delete CaloClusterCells_position_z;
    }
  }
  if (m_doTopo) {
    if (m_readTopoClusters) {
      delete CaloTopoClusters_energy;
      delete CaloTopoClusters_position_x;
      delete CaloTopoClusters_position_y;
      delete CaloTopoClusters_position_z;
      delete CaloTopoClusters_hits_begin;
      delete CaloTopoClusters_hits_end;
    }
    if (m_readCaloClusterCells) {
      delete CaloTopoClusterCells_cellID;
      delete CaloTopoClusterCells_energy;
      delete CaloTopoClusterCells_position_x;
      delete CaloTopoClusterCells_position_y;
      delete CaloTopoClusterCells_position_z;
    }
  }
  delete fReader;
}


void EventReader::loadEvent(int event) {
  if (event<0 || event>=nEvents) {
    cout << "Event number " << event << " out of range : [0, " << nEvents-1 << "]" << endl;
  }
  fReader->SetEntry(event);
}
