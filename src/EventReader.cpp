/******************************************************************************/
// Simple event display for the ALLEGRO detector with ECAL with inclined modules
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

bool EventReader::checkCollection(const std::string& drawFlag,
				  const std::string& flagName,
				  const std::string& requiredBranch,
				  std::string& collectionName)
{
  if (!displayConfig.getBoolConfig(drawFlag))
    return false;

  collectionName = displayConfig.getStringConfig(flagName);
  if (collectionName == "") {
    std::cout << "WARNING: " << flagName <<  " not set, disabling " << flagName << std::endl;
    displayConfig.setBoolConfig(drawFlag, false);
    return false;
  }
  else {
    if (! fReader->GetTree()->FindBranch(Form("%s.%s", collectionName.c_str(), requiredBranch.c_str()))) {
      std::cout << "WARNING: branch " << collectionName << "." << requiredBranch << " not found, disabling " << flagName << std::endl;
      displayConfig.setBoolConfig(drawFlag, false);
      displayConfig.setStringConfig(flagName, "");
      return false;
    }
    else
      return true;
  }
}

bool EventReader::checkCollections(const std::string& drawFlag,
				   const std::string& flagName1,
				   const std::string& flagName2,
				   const std::string& requiredBranch,
				   std::string& collectionName1,
				   std::string& collectionName2)
{
  if (!displayConfig.getBoolConfig(drawFlag))
    return false;

  collectionName1 = displayConfig.getStringConfig(flagName1);
  collectionName2 = displayConfig.getStringConfig(flagName2);
  if (collectionName1 == "" or collectionName2 == "") {
    std::cout << "WARNING: at least one of " << flagName1 << " or " << flagName2 << " not set, disabling " << flagName1 << " and " << flagName2 << std::endl;
    displayConfig.setBoolConfig(drawFlag, false);
    displayConfig.setStringConfig(flagName1, "");
    displayConfig.setStringConfig(flagName2, "");
    return false;
  }
  else {
    if (! fReader->GetTree()->FindBranch(Form("%s.%s", collectionName1.c_str(), requiredBranch.c_str())) or
	! fReader->GetTree()->FindBranch(Form("%s.%s", collectionName2.c_str(), requiredBranch.c_str())) ) {
      std::cout << "WARNING: branch " << collectionName1 << "." << requiredBranch << " or "
		<< collectionName2 << "." << requiredBranch << " not found, disabling " << flagName1 << " and " << flagName2 << std::endl;
      displayConfig.setBoolConfig(drawFlag, false);
      displayConfig.setStringConfig(flagName1, "");
      displayConfig.setStringConfig(flagName2, "");
      return false;
    }
    else
      return true;
  }
}

void EventReader::setBranches()
{
  std::string collection, collection2;

  // primary particles
  if (checkCollection("drawGenParticles", "genParticles", "PDG", collection))
  {
    genParticles_PDG = makeArray<Int_t>(collection, "PDG");
    genParticles_generatorStatus = makeArray<Int_t>(collection, "generatorStatus");
    genParticles_simulatorStatus = makeArray<Int_t>(collection, "simulatorStatus");
    genParticles_charge = makeArray<Float_t>(collection, "charge");
    genParticles_time = makeArray<Float_t>(collection, "time");
    genParticles_mass = makeArray<Double_t>(collection, "mass");
    genParticles_vertex_x = makeArray<Double_t>(collection, "vertex.x");
    genParticles_vertex_y = makeArray<Double_t>(collection, "vertex.y");
    genParticles_vertex_z = makeArray<Double_t>(collection, "vertex.z");
    genParticles_endpoint_x = makeArray<Double_t>(collection, "endpoint.x");
    genParticles_endpoint_y = makeArray<Double_t>(collection, "endpoint.y");
    genParticles_endpoint_z = makeArray<Double_t>(collection, "endpoint.z");
    genParticles_momentum_x = makeArray<Double_t>(collection, "momentum.x");
    genParticles_momentum_y = makeArray<Double_t>(collection, "momentum.y");
    genParticles_momentum_z = makeArray<Double_t>(collection, "momentum.z");
    genParticles_parents_begin = makeArray<UInt_t>(collection, "parents_begin");
    genParticles_parents_end = makeArray<UInt_t>(collection, "parents_end");
    genParticles_parents_index = makeArray<Int_t>("_" + collection + "_parents", "index");
    genParticles_daughters_begin = makeArray<UInt_t>(collection, "daughters_begin");
    genParticles_daughters_end = makeArray<UInt_t>(collection, "daughters_end");
    genParticles_daughters_index = makeArray<Int_t>("_" + collection + "_daughters", "index");
  }
  
  // secondary particles
  if (checkCollection("drawSimParticles", "simParticles", "PDG", collection))
  {
    SimParticleSecondaries_PDG = makeArray<Int_t>(collection, "PDG");
    //SimParticleSecondaries_generatorStatus = makeArray<Int_t>(collection, "generatorStatus");
    //SimParticleSecondaries_simulatorStatus = makeArray<Int_t>(collection, "simulatorStatus");
    //SimParticleSecondaries_charge = makeArray<Float_t>(collection, "charge");
    SimParticleSecondaries_time = makeArray<Float_t>(collection, "time");
    SimParticleSecondaries_mass = makeArray<Double_t>(collection, "mass");
    SimParticleSecondaries_vertex_x = makeArray<Double_t>(collection, "vertex.x");
    SimParticleSecondaries_vertex_y = makeArray<Double_t>(collection, "vertex.y");
    SimParticleSecondaries_vertex_z = makeArray<Double_t>(collection, "vertex.z");
    SimParticleSecondaries_endpoint_x = makeArray<Double_t>(collection, "endpoint.x");
    SimParticleSecondaries_endpoint_y = makeArray<Double_t>(collection, "endpoint.y");
    SimParticleSecondaries_endpoint_z = makeArray<Double_t>(collection, "endpoint.z");
    SimParticleSecondaries_momentum_x = makeArray<Double_t>(collection, "momentum.x");
    SimParticleSecondaries_momentum_y = makeArray<Double_t>(collection, "momentum.y");
    SimParticleSecondaries_momentum_z = makeArray<Double_t>(collection, "momentum.z");
  }

  // reconstructed tracks
  if (checkCollection("drawTracks", "tracks", "trackStates_begin", collection))
  {
    std::string trackStateColl = "_" + collection + "_trackStates";
    Tracks_trackStates_begin = makeArray<UInt_t>(collection, "trackStates_begin");
    Tracks_trackStates_end = makeArray<UInt_t>(collection, "trackStates_end");
    Tracks_subdetectorHitNumbers_begin = makeArray<UInt_t>(collection, "subdetectorHitNumbers_begin");
    Tracks_subdetectorHitNumbers_end = makeArray<UInt_t>(collection, "subdetectorHitNumbers_end");
    _Tracks_trackStates_location = makeArray<Int_t>(trackStateColl, "location");
    _Tracks_trackStates_omega = makeArray<Float_t>(trackStateColl, "omega");
    _Tracks_trackStates_phi = makeArray<Float_t>(trackStateColl, "phi");
    _Tracks_trackStates_tanLambda = makeArray<Float_t>(trackStateColl, "tanLambda");
    _Tracks_trackStates_referencePoint_x = makeArray<Float_t>(trackStateColl, "referencePoint.x");
    _Tracks_trackStates_referencePoint_y = makeArray<Float_t>(trackStateColl, "referencePoint.y");
    _Tracks_trackStates_referencePoint_z = makeArray<Float_t>(trackStateColl, "referencePoint.z");
    _Tracks_subdetectorHitNumbers = makeArray<Int_t>("", "_" + collection + "_subdetectorHitNumbers");
  }

  // hits in vertex detector
  if (checkCollections("drawVertexHits", "vertexBarrelHits", "vertexEndcapHits", "cellID", collection, collection2))
  {
    VertexBarrelHits_cellID     = makeArray<ULong_t>(collection, "cellID");
    VertexBarrelHits_energy     = makeArray<Float_t>(collection, "eDep");
    VertexBarrelHits_position_x = makeArray<Double_t>(collection, "position.x");
    VertexBarrelHits_position_y = makeArray<Double_t>(collection, "position.y");
    VertexBarrelHits_position_z = makeArray<Double_t>(collection, "position.z");

    VertexEndcapHits_cellID     = makeArray<ULong_t>(collection2, "cellID");
    VertexEndcapHits_energy     = makeArray<Float_t>(collection2, "eDep");
    VertexEndcapHits_position_x = makeArray<Double_t>(collection2, "position.x");
    VertexEndcapHits_position_y = makeArray<Double_t>(collection2, "position.y");
    VertexEndcapHits_position_z = makeArray<Double_t>(collection2, "position.z");
  }

  // hits in DCH or STT
  if (checkCollection("drawMainTrackerHits", "mainTrackerHits", "cellID", collection))
  {
    MainTrackerHits_cellID     = makeArray<ULong_t>(collection, "cellID");
    MainTrackerHits_energy     = makeArray<Float_t>(collection, "eDep");
    MainTrackerHits_position_x = makeArray<Double_t>(collection, "position.x");
    MainTrackerHits_position_y = makeArray<Double_t>(collection, "position.y");
    MainTrackerHits_position_z = makeArray<Double_t>(collection, "position.z");
  }

  // hits in si wrapper
  if (checkCollections("drawSiWrapperHits", "siWrapperBarrelHits", "siWrapperEndcapHits", "cellID", collection, collection2))
  {
    SiWrapperBarrelHits_cellID     = makeArray<ULong_t>(collection, "cellID");
    SiWrapperBarrelHits_energy     = makeArray<Float_t>(collection, "eDep");
    SiWrapperBarrelHits_position_x = makeArray<Double_t>(collection, "position.x");
    SiWrapperBarrelHits_position_y = makeArray<Double_t>(collection, "position.y");
    SiWrapperBarrelHits_position_z = makeArray<Double_t>(collection, "position.z");

    SiWrapperEndcapHits_cellID     = makeArray<ULong_t>(collection2, "cellID");
    SiWrapperEndcapHits_energy     = makeArray<Float_t>(collection2, "eDep");
    SiWrapperEndcapHits_position_x = makeArray<Double_t>(collection2, "position.x");
    SiWrapperEndcapHits_position_y = makeArray<Double_t>(collection2, "position.y");
    SiWrapperEndcapHits_position_z = makeArray<Double_t>(collection2, "position.z");
  }

  // digis in vertex detector
  if (checkCollections("drawVertexDigis", "vertexBarrelDigis", "vertexEndcapDigis", "cellID", collection, collection2))
  {
    VertexBarrelDigis_cellID     = makeArray<ULong_t>(collection, "cellID");
    VertexBarrelDigis_energy     = makeArray<Float_t>(collection, "eDep");
    VertexBarrelDigis_position_x = makeArray<Double_t>(collection, "position.x");
    VertexBarrelDigis_position_y = makeArray<Double_t>(collection, "position.y");
    VertexBarrelDigis_position_z = makeArray<Double_t>(collection, "position.z");

    VertexEndcapDigis_cellID     = makeArray<ULong_t>(collection2, "cellID");
    VertexEndcapDigis_energy     = makeArray<Float_t>(collection2, "eDep");
    VertexEndcapDigis_position_x = makeArray<Double_t>(collection2, "position.x");
    VertexEndcapDigis_position_y = makeArray<Double_t>(collection2, "position.y");
    VertexEndcapDigis_position_z = makeArray<Double_t>(collection2, "position.z");
  }

  // digis in DCH / STT
  if (checkCollection("drawMainTrackerDigis", "mainTrackerDigis", "cellID", collection))
  {
    MainTrackerDigis_cellID     = makeArray<ULong_t>(collection, "cellID");
    MainTrackerDigis_energy     = makeArray<Float_t>(collection, "eDep");
    MainTrackerDigis_position_x = makeArray<Double_t>(collection, "position.x");
    MainTrackerDigis_position_y = makeArray<Double_t>(collection, "position.y");
    MainTrackerDigis_position_z = makeArray<Double_t>(collection, "position.z");
  }

  // digis in si wrapper
  if (checkCollections("drawSiWrapperDigis", "siWrapperBarrelDigis", "siWrapperEndcapDigis", "cellID", collection, collection2))
  {
    SiWrapperBarrelDigis_cellID     = makeArray<ULong_t>(collection, "cellID");
    SiWrapperBarrelDigis_energy     = makeArray<Float_t>(collection, "eDep");
    SiWrapperBarrelDigis_position_x = makeArray<Double_t>(collection, "position.x");
    SiWrapperBarrelDigis_position_y = makeArray<Double_t>(collection, "position.y");
    SiWrapperBarrelDigis_position_z = makeArray<Double_t>(collection, "position.z");

    SiWrapperEndcapDigis_cellID     = makeArray<ULong_t>(collection2, "cellID");
    SiWrapperEndcapDigis_energy     = makeArray<Float_t>(collection2, "eDep");
    SiWrapperEndcapDigis_position_x = makeArray<Double_t>(collection2, "position.x");
    SiWrapperEndcapDigis_position_y = makeArray<Double_t>(collection2, "position.y");
    SiWrapperEndcapDigis_position_z = makeArray<Double_t>(collection2, "position.z");
  }

  // hits in ECal barrel
  if (checkCollection("drawECalBarrelHits", "ecalBarrelHits", "energy", collection))
  {
    ECalBarrelHits_energy     = makeArray<Float_t>(collection, "energy");
    ECalBarrelHits_position_x = makeArray<Float_t>(collection, "stepPosition.x");
    ECalBarrelHits_position_y = makeArray<Float_t>(collection, "stepPosition.y");
    ECalBarrelHits_position_z = makeArray<Float_t>(collection, "stepPosition.z");
  }

  // cells in ECal barrel
  if (checkCollection("drawECalBarrelCells", "ecalBarrelCells", "energy", collection))
  {
    ECalBarrelCells_cellID     = makeArray<ULong_t>(collection, "cellID");
    ECalBarrelCells_energy     = makeArray<Float_t>(collection, "energy");
    ECalBarrelCells_position_x = makeArray<Float_t>(collection, "position.x");
    ECalBarrelCells_position_y = makeArray<Float_t>(collection, "position.y");
    ECalBarrelCells_position_z = makeArray<Float_t>(collection, "position.z");
  }
  
  // cells in ECal barrel with coarser merging
  if (checkCollection("drawECalBarrelMergedCells", "ecalBarrelMergedCells", "energy", collection))
  {
    // ECalBarrelCells2_cellID     = makeArray<ULong_t>(collection, "cellID");    
    ECalBarrelCells2_energy     = makeArray<Float_t>(collection, "energy");    
    ECalBarrelCells2_position_x = makeArray<Float_t>(collection, "position.x");
    ECalBarrelCells2_position_y = makeArray<Float_t>(collection, "position.y");
    ECalBarrelCells2_position_z = makeArray<Float_t>(collection, "position.z");
  }

  // hits in ECal endcap
  if (checkCollection("drawECalEndcapHits", "ecalEndcapHits", "energy", collection))
  {
    ECalEndcapHits_energy     = makeArray<Float_t>(collection, "energy");
    ECalEndcapHits_position_x = makeArray<Float_t>(collection, "stepPosition.x");
    ECalEndcapHits_position_y = makeArray<Float_t>(collection, "stepPosition.y");
    ECalEndcapHits_position_z = makeArray<Float_t>(collection, "stepPosition.z");
  }

  // cells in ECal endcap
  if (checkCollection("drawECalEndcapCells", "ecalEndcapCells", "energy", collection))
  {
    ECalEndcapCells_cellID     = makeArray<ULong_t>(collection, "cellID");
    ECalEndcapCells_energy     = makeArray<Float_t>(collection, "energy");
    ECalEndcapCells_position_x = makeArray<Float_t>(collection, "position.x");
    ECalEndcapCells_position_y = makeArray<Float_t>(collection, "position.y");
    ECalEndcapCells_position_z = makeArray<Float_t>(collection, "position.z");
  }
  
  if (m_doHCal) {
    // hits in HCal barrel
    if (checkCollection("drawHCalBarrelHits", "hcalBarrelHits", "energy", collection))
    {
      HCalBarrelHits_energy     = makeArray<Float_t>(collection, "energy");    
      HCalBarrelHits_position_x = makeArray<Float_t>(collection, "stepPosition.x");
      HCalBarrelHits_position_y = makeArray<Float_t>(collection, "stepPosition.y");
      HCalBarrelHits_position_z = makeArray<Float_t>(collection, "stepPosition.z");
    }
    
    // cells in HCal barrel
    if (checkCollection("drawHCalBarrelCells", "hcalBarrelCells", "energy", collection))
    {
      HCalBarrelCells_cellID     = makeArray<ULong_t>(collection, "cellID");
      HCalBarrelCells_energy     = makeArray<Float_t>(collection, "energy");
      HCalBarrelCells_position_x = makeArray<Float_t>(collection, "position.x");
      HCalBarrelCells_position_y = makeArray<Float_t>(collection, "position.y");
      HCalBarrelCells_position_z = makeArray<Float_t>(collection, "position.z");
    }

    // hits in HCal endcap
    if (checkCollection("drawHCalEndcapHits", "hcalEndcapHits", "energy", collection))
    {
      HCalEndcapHits_energy     = makeArray<Float_t>(collection, "energy");    
      HCalEndcapHits_position_x = makeArray<Float_t>(collection, "stepPosition.x");
      HCalEndcapHits_position_y = makeArray<Float_t>(collection, "stepPosition.y");
      HCalEndcapHits_position_z = makeArray<Float_t>(collection, "stepPosition.z");
    }
    
    // cells in HCal endcap
    if (checkCollection("drawHCalEndcapCells", "hcalEndcapCells", "energy", collection))
    {
      HCalEndcapCells_cellID     = makeArray<ULong_t>(collection, "cellID");
      HCalEndcapCells_energy     = makeArray<Float_t>(collection, "energy");
      HCalEndcapCells_position_x = makeArray<Float_t>(collection, "position.x");
      HCalEndcapCells_position_y = makeArray<Float_t>(collection, "position.y");
      HCalEndcapCells_position_z = makeArray<Float_t>(collection, "position.z");
    }
  }
  else {
    std::cout << "doHCal is false, setting drawHCal*Hits and drawHCal*Cells to false" << std::endl;
    displayConfig.setBoolConfig("drawHCalBarrelHits", false);
    displayConfig.setStringConfig("hcalBarrelHits", "");
    displayConfig.setBoolConfig("drawHCalBarrelCells", false);
    displayConfig.setStringConfig("hcalBarrelCells", "");
    displayConfig.setBoolConfig("drawHCalEndcapHits", false);
    displayConfig.setStringConfig("hcalEndcapHits", "");
    displayConfig.setBoolConfig("drawHCalEndcapCells", false);
    displayConfig.setStringConfig("hcalEndcapCells", "");
  }

  // hits in muon tagger
  if (checkCollections("drawMuonHits", "muonBarrelHits", "muonEndcapHits", "energy", collection, collection2))
  {
    MuonBarrelHits_energy     = makeArray<Float_t>(collection, "energy");
    MuonBarrelHits_position_x = makeArray<Float_t>(collection, "stepPosition.x");
    MuonBarrelHits_position_y = makeArray<Float_t>(collection, "stepPosition.y");
    MuonBarrelHits_position_z = makeArray<Float_t>(collection, "stepPosition.z");

    MuonEndcapHits_energy     = makeArray<Float_t>(collection2, "energy");
    MuonEndcapHits_position_x = makeArray<Float_t>(collection2, "stepPosition.x");
    MuonEndcapHits_position_y = makeArray<Float_t>(collection2, "stepPosition.y");
    MuonEndcapHits_position_z = makeArray<Float_t>(collection2, "stepPosition.z");
  }

  // cells in muon tagger
  if (checkCollections("drawMuonCells", "muonBarrelCells", "muonEndcapCells", "energy", collection, collection2))
  {
    MuonBarrelCells_cellID     = makeArray<ULong_t>(collection, "cellID");
    MuonBarrelCells_energy     = makeArray<Float_t>(collection, "energy");
    MuonBarrelCells_position_x = makeArray<Float_t>(collection, "position.x");
    MuonBarrelCells_position_y = makeArray<Float_t>(collection, "position.y");
    MuonBarrelCells_position_z = makeArray<Float_t>(collection, "position.z");

    MuonEndcapCells_cellID     = makeArray<ULong_t>(collection2, "cellID");
    MuonEndcapCells_energy     = makeArray<Float_t>(collection2, "energy");
    MuonEndcapCells_position_x = makeArray<Float_t>(collection2, "position.x");
    MuonEndcapCells_position_y = makeArray<Float_t>(collection2, "position.y");
    MuonEndcapCells_position_z = makeArray<Float_t>(collection2, "position.z");
  }


  // SW clusters
  if (checkCollection("drawCaloClusters", "caloClusters", "energy", collection))
  {
    CaloClusters_energy     = makeArray<Float_t>(collection, "energy");
    CaloClusters_position_x = makeArray<Float_t>(collection, "position.x");
    CaloClusters_position_y = makeArray<Float_t>(collection, "position.y");
    CaloClusters_position_z = makeArray<Float_t>(collection, "position.z");
    CaloClusters_theta      = makeArray<Float_t>(collection, "iTheta");
    CaloClusters_phi        = makeArray<Float_t>(collection, "phi");
    CaloClusters_hits_begin = makeArray<UInt_t> (collection, "hits_begin");
    CaloClusters_hits_end   = makeArray<UInt_t> (collection, "hits_end");
  }

  // cells in the SW clusters
  // need to modify EventDisplay.cpp FillClusters and DrawClusters if I don't want to use cell info
  //  if (displayConfig.getBoolConfig("drawCaloClusters") && displayConfig.getBoolConfig("drawCaloClusterCells")) {
  if (displayConfig.getBoolConfig("drawCaloClusters")) {
    std::string collection = displayConfig.getStringConfig("caloClusterCells");
    if (collection == "") {
      // try to use default collection
      collection = displayConfig.getStringConfig("caloClusters");
      collection.pop_back();
      collection += "Cells";
      std::cout << "WARNING: caloClusterCells not set, trying with " << collection << std::endl;
      displayConfig.setStringConfig("caloClusterCells", collection);
    }
    if (! fReader->GetTree()->FindBranch(Form("%s.cellID", collection.c_str()))) {
      std::cout << "WARNING: branch " << collection << ".cellID not found, disabling calo cluster cells" << std::endl;
      displayConfig.setBoolConfig("drawCaloClusterCells", false);
      displayConfig.setStringConfig("caloClusterCells", "");
      // also disable the clusters until I fix EventDisplay.cpp
      displayConfig.setBoolConfig("drawCaloClusters", false);
      displayConfig.setStringConfig("caloClusters", "");
    }
    else {
      CaloClusterCells_cellID     = makeArray<ULong_t>(collection, "cellID");
      CaloClusterCells_energy     = makeArray<Float_t>(collection, "energy");
      CaloClusterCells_position_x = makeArray<Float_t>(collection, "position.x");
      CaloClusterCells_position_y = makeArray<Float_t>(collection, "position.y");
      CaloClusterCells_position_z = makeArray<Float_t>(collection, "position.z");
    }
  }

  // topo clusters
  if (checkCollection("drawTopoClusters", "topoClusters", "energy", collection))
  {
    CaloTopoClusters_energy     = makeArray<Float_t>(collection, "energy");
    CaloTopoClusters_position_x = makeArray<Float_t>(collection, "position.x");
    CaloTopoClusters_position_y = makeArray<Float_t>(collection, "position.y");
    CaloTopoClusters_position_z = makeArray<Float_t>(collection, "position.z");
    CaloTopoClusters_theta      = makeArray<Float_t>(collection, "iTheta");
    CaloTopoClusters_phi        = makeArray<Float_t>(collection, "phi");
    CaloTopoClusters_hits_begin = makeArray<UInt_t> (collection, "hits_begin");
    CaloTopoClusters_hits_end   = makeArray<UInt_t> (collection, "hits_end");
  }
  
  // cells in the topo clusters
  // need to modify EventDisplay.cpp FillClusters and DrawClusters if I don't want to use cell info
  // if (displayConfig.getBoolConfig("drawTopoClusters") && displayConfig.getBoolConfig("drawTopoClusterCells")) {
  if (displayConfig.getBoolConfig("drawTopoClusters")) {
    std::string collection = displayConfig.getStringConfig("topoClusterCells");
    if (collection == "") {
      // try to use default collection
      collection = displayConfig.getStringConfig("topoClusters");
      collection.pop_back();
      collection += "Cells";
      std::cout << "WARNING: topoClusterCells not set, trying with " << collection << std::endl;
      displayConfig.setStringConfig("topoClusterCells", collection);
    }
    if (! fReader->GetTree()->FindBranch(Form("%s.cellID", collection.c_str()))) {
      std::cout << "WARNING: branch " << collection << ".cellID not found, disabling topo cluster cells" << std::endl;
      displayConfig.setBoolConfig("drawTopoClusterCells", false);
      displayConfig.setStringConfig("topoClusterCells", "");
      // also disable the clusters until I fix EventDisplay.cpp
      displayConfig.setBoolConfig("drawTopoClusters", false);
      displayConfig.setStringConfig("topoClusters", "");
    }
    else {
      CaloTopoClusterCells_cellID     = makeArray<ULong_t>(collection, "cellID");
      CaloTopoClusterCells_energy     = makeArray<Float_t>(collection, "energy");
      CaloTopoClusterCells_position_x = makeArray<Float_t>(collection, "position.x");
      CaloTopoClusterCells_position_y = makeArray<Float_t>(collection, "position.y");
      CaloTopoClusterCells_position_z = makeArray<Float_t>(collection, "position.z");
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
    delete genParticles_endpoint_x;
    delete genParticles_endpoint_y;
    delete genParticles_endpoint_z;
    delete genParticles_momentum_x;
    delete genParticles_momentum_y;
    delete genParticles_momentum_z;
    delete genParticles_parents_begin;
    delete genParticles_parents_end;
    delete genParticles_parents_index;
    delete genParticles_daughters_begin;
    delete genParticles_daughters_end;
    delete genParticles_daughters_index;
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
  if (Tracks_trackStates_begin) {
    delete Tracks_trackStates_begin;
    delete Tracks_trackStates_end;
    delete Tracks_subdetectorHitNumbers_begin;
    delete Tracks_subdetectorHitNumbers_end;
    delete _Tracks_trackStates_location;
    delete _Tracks_trackStates_omega;
    delete _Tracks_trackStates_phi;
    delete _Tracks_trackStates_tanLambda;
    delete _Tracks_trackStates_referencePoint_x;
    delete _Tracks_trackStates_referencePoint_y;
    delete _Tracks_trackStates_referencePoint_z;
    delete _Tracks_subdetectorHitNumbers;
  }
  if (VertexBarrelHits_cellID) {
    delete VertexBarrelHits_cellID;
    delete VertexBarrelHits_energy;
    delete VertexBarrelHits_position_x;
    delete VertexBarrelHits_position_y;
    delete VertexBarrelHits_position_z;
  }
  if (VertexEndcapHits_cellID) {
    delete VertexEndcapHits_cellID;
    delete VertexEndcapHits_energy;
    delete VertexEndcapHits_position_x;
    delete VertexEndcapHits_position_y;
    delete VertexEndcapHits_position_z;
  }
  if (MainTrackerHits_cellID) {
    delete MainTrackerHits_cellID;
    delete MainTrackerHits_energy;
    delete MainTrackerHits_position_x;
    delete MainTrackerHits_position_y;
    delete MainTrackerHits_position_z;
  }
  if (SiWrapperBarrelHits_cellID) {
    delete SiWrapperBarrelHits_cellID;
    delete SiWrapperBarrelHits_energy;
    delete SiWrapperBarrelHits_position_x;
    delete SiWrapperBarrelHits_position_y;
    delete SiWrapperBarrelHits_position_z;
  }
  if (SiWrapperEndcapHits_cellID) {
    delete SiWrapperEndcapHits_cellID;
    delete SiWrapperEndcapHits_energy;
    delete SiWrapperEndcapHits_position_x;
    delete SiWrapperEndcapHits_position_y;
    delete SiWrapperEndcapHits_position_z;
  }
  if (VertexBarrelDigis_cellID) {
    delete VertexBarrelDigis_cellID;
    delete VertexBarrelDigis_energy;
    delete VertexBarrelDigis_position_x;
    delete VertexBarrelDigis_position_y;
    delete VertexBarrelDigis_position_z;
  }
  if (VertexEndcapDigis_cellID) {
    delete VertexEndcapDigis_cellID;
    delete VertexEndcapDigis_energy;
    delete VertexEndcapDigis_position_x;
    delete VertexEndcapDigis_position_y;
    delete VertexEndcapDigis_position_z;
  }
  if (MainTrackerDigis_cellID) {
    delete MainTrackerDigis_cellID;
    delete MainTrackerDigis_energy;
    delete MainTrackerDigis_position_x;
    delete MainTrackerDigis_position_y;
    delete MainTrackerDigis_position_z;
  }
  if (SiWrapperBarrelDigis_cellID) {
    delete SiWrapperBarrelDigis_cellID;
    delete SiWrapperBarrelDigis_energy;
    delete SiWrapperBarrelDigis_position_x;
    delete SiWrapperBarrelDigis_position_y;
    delete SiWrapperBarrelDigis_position_z;
  }
  if (SiWrapperEndcapDigis_cellID) {
    delete SiWrapperEndcapDigis_cellID;
    delete SiWrapperEndcapDigis_energy;
    delete SiWrapperEndcapDigis_position_x;
    delete SiWrapperEndcapDigis_position_y;
    delete SiWrapperEndcapDigis_position_z;
  }
  if (ECalBarrelHits_energy) {
    delete ECalBarrelHits_energy;
    delete ECalBarrelHits_position_x;
    delete ECalBarrelHits_position_y;
    delete ECalBarrelHits_position_z;
  }
  if (ECalBarrelCells_cellID) {
    delete ECalBarrelCells_cellID;
    delete ECalBarrelCells_energy;
    delete ECalBarrelCells_position_x;
    delete ECalBarrelCells_position_y;
    delete ECalBarrelCells_position_z;
  }
  if (ECalBarrelCells2_cellID) {
    delete ECalBarrelCells2_cellID;
    delete ECalBarrelCells2_energy;
    delete ECalBarrelCells2_position_x;
    delete ECalBarrelCells2_position_y;
    delete ECalBarrelCells2_position_z;
  }
  if (ECalEndcapHits_energy) {
    delete ECalEndcapHits_energy;
    delete ECalEndcapHits_position_x;
    delete ECalEndcapHits_position_y;
    delete ECalEndcapHits_position_z;
  }
  if (ECalEndcapCells_cellID) {
    delete ECalEndcapCells_cellID;
    delete ECalEndcapCells_energy;
    delete ECalEndcapCells_position_x;
    delete ECalEndcapCells_position_y;
    delete ECalEndcapCells_position_z;
  }
  if (HCalBarrelHits_energy) {
    delete HCalBarrelHits_energy;
    delete HCalBarrelHits_position_x;
    delete HCalBarrelHits_position_y;
    delete HCalBarrelHits_position_z;
  }
  if (HCalBarrelCells_cellID) {
    delete HCalBarrelCells_cellID;
    delete HCalBarrelCells_energy;
    delete HCalBarrelCells_position_x;
    delete HCalBarrelCells_position_y;
    delete HCalBarrelCells_position_z;
  }
  if (HCalEndcapHits_energy) {
    delete HCalEndcapHits_energy;
    delete HCalEndcapHits_position_x;
    delete HCalEndcapHits_position_y;
    delete HCalEndcapHits_position_z;
  }
  if (HCalEndcapCells_cellID) {
    delete HCalEndcapCells_cellID;
    delete HCalEndcapCells_energy;
    delete HCalEndcapCells_position_x;
    delete HCalEndcapCells_position_y;
    delete HCalEndcapCells_position_z;
  }
  if (MuonBarrelHits_cellID) {
    delete MuonBarrelHits_cellID;
    delete MuonBarrelHits_energy;
    delete MuonBarrelHits_position_x;
    delete MuonBarrelHits_position_y;
    delete MuonBarrelHits_position_z;
  }
  if (MuonEndcapHits_cellID) {
    delete MuonEndcapHits_cellID;
    delete MuonEndcapHits_energy;
    delete MuonEndcapHits_position_x;
    delete MuonEndcapHits_position_y;
    delete MuonEndcapHits_position_z;
  }
  if (MuonBarrelCells_cellID) {
    delete MuonBarrelCells_cellID;
    delete MuonBarrelCells_energy;
    delete MuonBarrelCells_position_x;
    delete MuonBarrelCells_position_y;
    delete MuonBarrelCells_position_z;
  }
  if (MuonEndcapCells_cellID) {
    delete MuonEndcapCells_cellID;
    delete MuonEndcapCells_energy;
    delete MuonEndcapCells_position_x;
    delete MuonEndcapCells_position_y;
    delete MuonEndcapCells_position_z;
  }
  if (CaloClusters_energy) {
    delete CaloClusters_energy;
    delete CaloClusters_position_x;
    delete CaloClusters_position_y;
    delete CaloClusters_position_z;
    delete CaloClusters_theta;
    delete CaloClusters_phi;
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
    delete CaloTopoClusters_theta;
    delete CaloTopoClusters_phi;
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
