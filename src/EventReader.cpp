/******************************************************************************/
// Simple event display for the LAr calorimeter with inclined modules
// author: Giovanni Marchiori (giovanni.marchiori@cern.ch)
//
// class EventReader: sets up the TTreeReader branches for reading event data
//
/******************************************************************************/

#include "EventReader.h"
#include "Globals.h"

EventReader::EventReader(TFile* f, bool doHCal)
{
  fReader = new TTreeReader("events", f);
  nEvents = fReader->GetEntries();
  m_doHCal = doHCal;
}

void EventReader::SetBranches()
{
  // primary particles
  if (displayConfig.getBoolConfig("drawGenParticles"))
  {
    std::string branchName = displayConfig.getStringConfig("genParticles");
    const char* branch = branchName.c_str();
    if (branchName == "") {
      std::cout << "WARNING: genParticles not set, disabling gen particles" << std::endl;
      displayConfig.setBoolConfig("drawGenParticles", false);
    }
    else {
      if (! fReader->GetTree()->FindBranch(Form("%s.PDG", branch)))
      {
        std::cout << "WARNING: branch " << branch << ".PDG not found, disabling gen particles" << std::endl;
        displayConfig.setBoolConfig("drawGenParticles", false);
        displayConfig.setStringConfig("genParticles", "");
      }
      else
      {
        genParticles_PDG = new TTreeReaderArray<Int_t>(*fReader, Form("%s.PDG", branch));
        genParticles_generatorStatus = new TTreeReaderArray<Int_t>(*fReader, Form("%s.generatorStatus", branch));
        genParticles_simulatorStatus = new TTreeReaderArray<Int_t>(*fReader, Form("%s.simulatorStatus", branch));
        genParticles_charge = new TTreeReaderArray<Float_t>(*fReader, Form("%s.charge", branch));
        genParticles_time = new TTreeReaderArray<Float_t>(*fReader, Form("%s.time", branch));
        genParticles_mass = new TTreeReaderArray<Double_t>(*fReader, Form("%s.mass", branch));
        genParticles_vertex_x = new TTreeReaderArray<Double_t>(*fReader, Form("%s.vertex.x", branch));
        genParticles_vertex_y = new TTreeReaderArray<Double_t>(*fReader, Form("%s.vertex.y", branch));
        genParticles_vertex_z = new TTreeReaderArray<Double_t>(*fReader, Form("%s.vertex.z", branch));
        //genParticles_endpoint_x = new TTreeReaderArray<Double_t>(*fReader, Form("%s.endpoint.x", branch));
        //genParticles_endpoint_y = new TTreeReaderArray<Double_t>(*fReader, Form("%s.endpoint.y", branch));
        //genParticles_endpoint_z = new TTreeReaderArray<Double_t>(*fReader, Form("%s.endpoint.z", branch));
        // these variables are double in the latest version of the EDM
        // genParticles_momentum_x = new TTreeReaderArray<Float_t>(*fReader, Form("%s.momentum.x", branch));
        // genParticles_momentum_y = new TTreeReaderArray<Float_t>(*fReader, Form("%s.momentum.y", branch));
        // genParticles_momentum_z = new TTreeReaderArray<Float_t>(*fReader, Form("%s.momentum.z", branch));
        genParticles_momentum_x = new TTreeReaderArray<Double_t>(*fReader, Form("%s.momentum.x", branch));
        genParticles_momentum_y = new TTreeReaderArray<Double_t>(*fReader, Form("%s.momentum.y", branch));
        genParticles_momentum_z = new TTreeReaderArray<Double_t>(*fReader, Form("%s.momentum.z", branch));
      }
    }
  }
  
  // secondary particles
  if (displayConfig.getBoolConfig("drawSimParticles"))
  {
    std::string branchName = displayConfig.getStringConfig("simParticles");
    const char* branch = branchName.c_str();
    if (branchName == "") {
      std::cout << "WARNING: simParticles not set, disabling sim particles" << std::endl;
      displayConfig.setBoolConfig("drawSimParticles", false);
    }
    else {
      if (! fReader->GetTree()->FindBranch(Form("%s.PDG", branch))) {
        std::cout << "WARNING: branch " << branch << ".PDG not found, disabling sim particles" << std::endl;
        displayConfig.setBoolConfig("drawSimParticles", false);
        displayConfig.setStringConfig("simParticles", "");
      }
      else
      {
        SimParticleSecondaries_PDG = new TTreeReaderArray<Int_t>(*fReader, Form("%s.PDG", branch));
        //SimParticleSecondaries_generatorStatus = new TTreeReaderArray<Int_t>(*fReader, Form("%s.generatorStatus", branch));
        //SimParticleSecondaries_simulatorStatus = new TTreeReaderArray<Int_t>(*fReader, Form("%s.simulatorStatus", branch));
        //SimParticleSecondaries_charge = new TTreeReaderArray<Float_t>(*fReader, Form("%s.charge", branch));
        SimParticleSecondaries_time = new TTreeReaderArray<Float_t>(*fReader, Form("%s.time", branch));
        SimParticleSecondaries_mass = new TTreeReaderArray<Double_t>(*fReader, Form("%s.mass", branch));
        SimParticleSecondaries_vertex_x = new TTreeReaderArray<Double_t>(*fReader, Form("%s.vertex.x", branch));
        SimParticleSecondaries_vertex_y = new TTreeReaderArray<Double_t>(*fReader, Form("%s.vertex.y", branch));
        SimParticleSecondaries_vertex_z = new TTreeReaderArray<Double_t>(*fReader, Form("%s.vertex.z", branch));
        SimParticleSecondaries_endpoint_x = new TTreeReaderArray<Double_t>(*fReader, Form("%s.endpoint.x", branch));
        SimParticleSecondaries_endpoint_y = new TTreeReaderArray<Double_t>(*fReader, Form("%s.endpoint.y", branch));
        SimParticleSecondaries_endpoint_z = new TTreeReaderArray<Double_t>(*fReader, Form("%s.endpoint.z", branch));
        SimParticleSecondaries_momentum_x = new TTreeReaderArray<Float_t>(*fReader, Form("%s.momentum.x", branch));
        SimParticleSecondaries_momentum_y = new TTreeReaderArray<Float_t>(*fReader, Form("%s.momentum.y", branch));
        SimParticleSecondaries_momentum_z = new TTreeReaderArray<Float_t>(*fReader, Form("%s.momentum.z", branch));
      }
    }
  }
  
  // hits in ECal barrel
  if (displayConfig.getBoolConfig("drawECalBarrelHits"))
  {
    std::string branchName = displayConfig.getStringConfig("ecalBarrelHits");
    const char* branch = branchName.c_str();
    if (branchName == "") {
      std::cout << "WARNING: ecalBarrelHits not set, disabling ecal barrel hits" << std::endl;
      displayConfig.setBoolConfig("drawECalBarrelHits", false);
    }
    else {
      if (! fReader->GetTree()->FindBranch(Form("%s.cellID", branch)))
      {
        std::cout << "WARNING: branch " << branch << ".cellID not found, disabling ecal barrel hits" << std::endl;
        displayConfig.setBoolConfig("drawECalBarrelHits", false);
        displayConfig.setStringConfig("ecalBarrelHits", "");
      }
      else
      {
        ECalBarrelPositionedHits_cellID     = new TTreeReaderArray<ULong_t>(*fReader, Form("%s.cellID", branch));
        ECalBarrelPositionedHits_energy     = new TTreeReaderArray<Float_t>(*fReader, Form("%s.energy", branch));
        ECalBarrelPositionedHits_position_x = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.x", branch));
        ECalBarrelPositionedHits_position_y = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.y", branch));
        ECalBarrelPositionedHits_position_z = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.z", branch));
      }
    }
  }
  
  // cells in ECal barrel
  if (displayConfig.getBoolConfig("drawECalBarrelCells"))
  {
    std::string branchName = displayConfig.getStringConfig("ecalBarrelCells");
    const char* branch = branchName.c_str();
    if (branchName == "") {
      std::cout << "WARNING: ecalBarrelCells not set, disabling ecal barrel cells" << std::endl;
      displayConfig.setBoolConfig("drawECalBarrelCells", false);
    }
    else {
      if (! fReader->GetTree()->FindBranch(Form("%s.cellID", branch)))
      {
        std::cout << "WARNING: branch " << branch << ".cellID not found, disabling ecal barrel cells" << std::endl;
        displayConfig.setBoolConfig("drawECalBarrelCells", false);
        displayConfig.setStringConfig("ecalBarrelCells", "");
      }
      else
      {
        ECalBarrelPositionedCells_cellID     = new TTreeReaderArray<ULong_t>(*fReader, Form("%s.cellID", branch));
        ECalBarrelPositionedCells_energy     = new TTreeReaderArray<Float_t>(*fReader, Form("%s.energy", branch));
        ECalBarrelPositionedCells_position_x = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.x", branch));
        ECalBarrelPositionedCells_position_y = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.y", branch));
        ECalBarrelPositionedCells_position_z = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.z", branch));
      }
    }
  }
  
  // cells in ECal barrel with coarser merging
  if (displayConfig.getBoolConfig("drawECalBarrelMergedCells"))
  {
    std::string branchName = displayConfig.getStringConfig("ecalBarrelMergedCells");
    const char* branch = branchName.c_str();
    if (branchName == "") {
      std::cout << "WARNING: ecalBarrelMergedCells not set, disabling ecal barrel merged cells" << std::endl;
      displayConfig.setBoolConfig("drawECalBarrelMergedCells", false);
    }
    else {
      if (! fReader->GetTree()->FindBranch(Form("%s.cellID", branch)))
      {
        std::cout << "WARNING: branch " << branch << ".cellID not found, disabling ecal barrel cells" << std::endl;
        displayConfig.setBoolConfig("drawECalBarrelMergedCells", false);
        displayConfig.setStringConfig("ecalBarrelMergedCells", "");
      }
      else
      {
        ECalBarrelPositionedCells2_cellID     = new TTreeReaderArray<ULong_t>(*fReader, Form("%s.cellID", branch));    
        ECalBarrelPositionedCells2_energy     = new TTreeReaderArray<Float_t>(*fReader, Form("%s.energy", branch));    
        ECalBarrelPositionedCells2_position_x = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.x", branch));
        ECalBarrelPositionedCells2_position_y = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.y", branch));
        ECalBarrelPositionedCells2_position_z = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.z", branch));
      }
    }
  }
  else {
    displayConfig.setStringConfig("ecalBarrelMergedCells", "");
  }

  if (m_doHCal) {
    // hits in HCal barrel
    if (displayConfig.getBoolConfig("drawHcalBarrelHits"))
    {
      std::string branchName = displayConfig.getStringConfig("hcalBarrelHits");
      const char* branch = branchName.c_str();
      if (branchName == "") {
        std::cout << "WARNING: hcalBarrelHits not set, disabling hcal barrel hits" << std::endl;
        displayConfig.setBoolConfig("drawHcalBarrelHits", false);
      }
      else {
        if (! fReader->GetTree()->FindBranch(Form("%s.cellID", branch)))
        {
          std::cout << "WARNING: branch " << branch << ".cellID not found, disabling hcal barrel hits" << std::endl;
          displayConfig.setBoolConfig("drawHcalBarrelHits", false);
          displayConfig.setStringConfig("hcalBarrelHits", "");
        }
        else
        {
          HCalBarrelPositionedHits_cellID     = new TTreeReaderArray<ULong_t>(*fReader, Form("%s.cellID", branch));    
          HCalBarrelPositionedHits_energy     = new TTreeReaderArray<Float_t>(*fReader, Form("%s.energy", branch));    
          HCalBarrelPositionedHits_position_x = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.x", branch));
          HCalBarrelPositionedHits_position_y = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.y", branch));
          HCalBarrelPositionedHits_position_z = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.z", branch));
        }
      }
    }
    
    // cells in HCal barrel
    if (displayConfig.getBoolConfig("drawHcalBarrelCells"))
    {
      std::string branchName = displayConfig.getStringConfig("hcalBarrelCells");
      const char* branch = branchName.c_str();
      if (branchName == "") {
       std::cout << "WARNING: hcalBarrelCells not set, disabling hcal barrel cells" << std::endl;
        displayConfig.setBoolConfig("drawHcalBarrelCells", false);
      }
      else {
        if (! fReader->GetTree()->FindBranch(Form("%s.cellID", branch)))
        {
          std::cout << "WARNING: branch " << branch << ".cellID not found, disabling hcal barrel cells" << std::endl;
          displayConfig.setBoolConfig("drawHcalBarrelCells", false);
          displayConfig.setStringConfig("hcalBarrelCells", "");
        }
        else
        {
          HCalBarrelPositionedCells_cellID     = new TTreeReaderArray<ULong_t>(*fReader, Form("%s.cellID", branch));
          HCalBarrelPositionedCells_energy     = new TTreeReaderArray<Float_t>(*fReader, Form("%s.energy", branch));
          HCalBarrelPositionedCells_position_x = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.x", branch));
          HCalBarrelPositionedCells_position_y = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.y", branch));
          HCalBarrelPositionedCells_position_z = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.z", branch));
        }
      }
    }
  }
  else {
    displayConfig.setBoolConfig("drawHCalBarrelHits", false);
    displayConfig.setStringConfig("hcalBarrelHits", "");
    displayConfig.setBoolConfig("drawHCalBarrelCells", false);
    displayConfig.setStringConfig("hcalBarrelCells", "");
  }
  
  // SW clusters
  if (displayConfig.getBoolConfig("drawCaloClusters")) {
    std::string branchName = displayConfig.getStringConfig("caloClusters");
    const char* branch = branchName.c_str();
    if (branchName == "") {
      std::cout << "WARNING: caloClusters not set, disabling calo clusters" << std::endl;
      displayConfig.setBoolConfig("drawCaloClusters", false);
    }
    else {
      if (! fReader->GetTree()->FindBranch(Form("%s.energy", branch)))
      {
	std::cout << "WARNING: branch " << branch << ".energy not found, disabling calo clusters" << std::endl;
	displayConfig.setBoolConfig("drawCaloClusters", false);
        displayConfig.setStringConfig("caloClusters", "");
      }
      else {
        CaloClusters_energy     = new TTreeReaderArray<Float_t>(*fReader, Form("%s.energy", branch));
        CaloClusters_position_x = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.x", branch));
        CaloClusters_position_y = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.y", branch));
        CaloClusters_position_z = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.z", branch));
        CaloClusters_hits_begin = new TTreeReaderArray<UInt_t> (*fReader, Form("%s.hits_begin", branch));
        CaloClusters_hits_end   = new TTreeReaderArray<UInt_t> (*fReader, Form("%s.hits_end", branch));
      }
    }
  }
  // cells in the SW clusters
  // need to modify EventDisplay.cpp FillClusters and DrawClusters if I don't want to use cell info
  //  if (displayConfig.getBoolConfig("drawCaloClusters") && displayConfig.getBoolConfig("drawCaloClusterCells")) {
  if (displayConfig.getBoolConfig("drawCaloClusters")) {
    std::string branchName = displayConfig.getStringConfig("caloClusterCells");
    if (branchName == "") {
      // try to use default collection
      branchName = displayConfig.getStringConfig("caloClusters");
      branchName.pop_back();
      branchName += "Cells";
      std::cout << "WARNING: caloClusterCells not set, trying with " << branchName << std::endl;
      displayConfig.setStringConfig("caloClusterCells", branchName);
    }
    const char* branch = branchName.c_str();
    if (! fReader->GetTree()->FindBranch(Form("%s.cellID", branch))) {
      std::cout << "WARNING: branch " << branch << ".cellID not found, disabling calo cluster cells" << std::endl;
      displayConfig.setBoolConfig("drawCaloClusterCells", false);
      displayConfig.setStringConfig("caloClusterCells", "");
      // also disable the clusters until I fix EventDisplay.cpp
      displayConfig.setBoolConfig("drawCaloClusters", false);
      displayConfig.setStringConfig("caloClusters", "");
    }
    else {
      CaloClusterCells_cellID     = new TTreeReaderArray<ULong_t>(*fReader, Form("%s.cellID", branch));
      CaloClusterCells_energy     = new TTreeReaderArray<Float_t>(*fReader, Form("%s.energy", branch));
      CaloClusterCells_position_x = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.x", branch));
      CaloClusterCells_position_y = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.y", branch));
      CaloClusterCells_position_z = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.z", branch));
    }
  }


  // topo clusters
  if (displayConfig.getBoolConfig("drawTopoClusters")) {
    std::string branchName = displayConfig.getStringConfig("topoClusters");
    const char* branch = branchName.c_str();
    if (branchName == "") {
      std::cout << "WARNING: topoClusters not set, disabling topo clusters" << std::endl;
      displayConfig.setBoolConfig("drawTopoClusters", false);
    }
    else {
      if (! fReader->GetTree()->FindBranch(Form("%s.energy", branch)))
      {
	std::cout << "WARNING: branch " << branch << ".energy not found, disabling topo clusters" << std::endl;
	displayConfig.setBoolConfig("drawTopoClusters", false);
        displayConfig.setStringConfig("topoClusters", "");
      }
      else {
        CaloTopoClusters_energy     = new TTreeReaderArray<Float_t>(*fReader, Form("%s.energy", branch));
        CaloTopoClusters_position_x = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.x", branch));
        CaloTopoClusters_position_y = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.y", branch));
        CaloTopoClusters_position_z = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.z", branch));
        CaloTopoClusters_hits_begin = new TTreeReaderArray<UInt_t> (*fReader, Form("%s.hits_begin", branch));
        CaloTopoClusters_hits_end   = new TTreeReaderArray<UInt_t> (*fReader, Form("%s.hits_end", branch));
      }
    }
  }
  
  // cells in the topo clusters
  // need to modify EventDisplay.cpp FillClusters and DrawClusters if I don't want to use cell info
  // if (displayConfig.getBoolConfig("drawTopoClusters") && displayConfig.getBoolConfig("drawTopoClusterCells")) {
  if (displayConfig.getBoolConfig("drawTopoClusters")) {
    std::string branchName = displayConfig.getStringConfig("topoClusterCells");
    if (branchName == "") {
      // try to use default collection
      branchName = displayConfig.getStringConfig("topoClusters");
      branchName.pop_back();
      branchName += "Cells";
      std::cout << "WARNING: topoClusterCells not set, trying with " << branchName << std::endl;
      displayConfig.setStringConfig("topoClusterCells", branchName);
    }
    const char* branch = branchName.c_str();
    if (! fReader->GetTree()->FindBranch(Form("%s.cellID", branch))) {
      std::cout << "WARNING: branch " << branch << ".cellID not found, disabling topo cluster cells" << std::endl;
      displayConfig.setBoolConfig("drawTopoClusterCells", false);
      displayConfig.setStringConfig("topoClusterCells", "");
      // also disable the clusters until I fix EventDisplay.cpp
      displayConfig.setBoolConfig("drawTopoClusters", false);
      displayConfig.setStringConfig("topoClusters", "");
    }
    else {
      CaloTopoClusterCells_cellID     = new TTreeReaderArray<ULong_t>(*fReader, Form("%s.cellID", branch));
      CaloTopoClusterCells_energy     = new TTreeReaderArray<Float_t>(*fReader, Form("%s.energy", branch));
      CaloTopoClusterCells_position_x = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.x", branch));
      CaloTopoClusterCells_position_y = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.y", branch));
      CaloTopoClusterCells_position_z = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.z", branch));
    }
  }
}

EventReader::~EventReader() {
  if (genParticles_PDG) {
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
  if (ECalBarrelPositionedHits_cellID) {
    delete ECalBarrelPositionedHits_cellID;
    delete ECalBarrelPositionedHits_energy;
    delete ECalBarrelPositionedHits_position_x;
    delete ECalBarrelPositionedHits_position_y;
    delete ECalBarrelPositionedHits_position_z;
  }
  if (ECalBarrelPositionedCells_cellID) {
    delete ECalBarrelPositionedCells_cellID;
    delete ECalBarrelPositionedCells_energy;
    delete ECalBarrelPositionedCells_position_x;
    delete ECalBarrelPositionedCells_position_y;
    delete ECalBarrelPositionedCells_position_z;
  }
  if (ECalBarrelPositionedCells2_cellID) {
    delete ECalBarrelPositionedCells2_cellID;
    delete ECalBarrelPositionedCells2_energy;
    delete ECalBarrelPositionedCells2_position_x;
    delete ECalBarrelPositionedCells2_position_y;
    delete ECalBarrelPositionedCells2_position_z;
  }
  if (HCalBarrelPositionedHits_cellID) {
    delete HCalBarrelPositionedHits_cellID;
    delete HCalBarrelPositionedHits_energy;
    delete HCalBarrelPositionedHits_position_x;
    delete HCalBarrelPositionedHits_position_y;
    delete HCalBarrelPositionedHits_position_z;
  }
  if (HCalBarrelPositionedCells_cellID) {
    delete HCalBarrelPositionedCells_cellID;
    delete HCalBarrelPositionedCells_energy;
    delete HCalBarrelPositionedCells_position_x;
    delete HCalBarrelPositionedCells_position_y;
    delete HCalBarrelPositionedCells_position_z;
  }
  if (CaloClusters_energy) {
    delete CaloClusters_energy;
    delete CaloClusters_position_x;
    delete CaloClusters_position_y;
    delete CaloClusters_position_z;
    delete CaloClusters_hits_begin;
    delete CaloClusters_hits_end;
  }
  if (CaloClusterCells_cellID) {
    delete CaloClusterCells_cellID;
    delete CaloClusterCells_energy;
    delete CaloClusterCells_position_x;
    delete CaloClusterCells_position_y;
    delete CaloClusterCells_position_z;
  }
  if (CaloTopoClusters_energy) {
    delete CaloTopoClusters_energy;
    delete CaloTopoClusters_position_x;
    delete CaloTopoClusters_position_y;
    delete CaloTopoClusters_position_z;
    delete CaloTopoClusters_hits_begin;
    delete CaloTopoClusters_hits_end;
  }
  if (CaloTopoClusterCells_cellID) {
    delete CaloTopoClusterCells_cellID;
    delete CaloTopoClusterCells_energy;
    delete CaloTopoClusterCells_position_x;
    delete CaloTopoClusterCells_position_y;
    delete CaloTopoClusterCells_position_z;
  }
  delete fReader;
}


void EventReader::loadEvent(int event) {
  if (event<0 || event>=nEvents) {
    cout << "Event number " << event << " out of range : [0, " << nEvents-1 << "]" << endl;
  }
  fReader->SetEntry(event);
}
