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
        genParticles_endpoint_x = new TTreeReaderArray<Double_t>(*fReader, Form("%s.endpoint.x", branch));
        genParticles_endpoint_y = new TTreeReaderArray<Double_t>(*fReader, Form("%s.endpoint.y", branch));
        genParticles_endpoint_z = new TTreeReaderArray<Double_t>(*fReader, Form("%s.endpoint.z", branch));
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
        SimParticleSecondaries_momentum_x = new TTreeReaderArray<Double_t>(*fReader, Form("%s.momentum.x", branch));
        SimParticleSecondaries_momentum_y = new TTreeReaderArray<Double_t>(*fReader, Form("%s.momentum.y", branch));
        SimParticleSecondaries_momentum_z = new TTreeReaderArray<Double_t>(*fReader, Form("%s.momentum.z", branch));
      }
    }
  }

  // reconstructed tracks
  if (displayConfig.getBoolConfig("drawTracks"))
  {
    std::string branchName = displayConfig.getStringConfig("tracks");
    const char* branch = branchName.c_str();
    if (branchName == "") {
      std::cout << "WARNING: tracks not set, disabling tracks" << std::endl;
      displayConfig.setBoolConfig("drawTracks", false);
    }
    else {
      if (! fReader->GetTree()->FindBranch(Form("%s.trackStates_begin", branch)))
      {
        std::cout << "WARNING: branch " << branch << ".trackStates_begin not found, disabling tracks" << std::endl;
        displayConfig.setBoolConfig("drawTracks", false);
        displayConfig.setStringConfig("tracks", "");
      }
      else
      {
	TracksFromGenParticles_trackStates_begin = new TTreeReaderArray<UInt_t>(*fReader, Form("%s.trackStates_begin", branch));
	TracksFromGenParticles_trackStates_end = new TTreeReaderArray<UInt_t>(*fReader, Form("%s.trackStates_end", branch));
	TracksFromGenParticles_subdetectorHitNumbers_begin = new TTreeReaderArray<UInt_t>(*fReader, Form("%s.subdetectorHitNumbers_begin", branch));
	TracksFromGenParticles_subdetectorHitNumbers_end = new TTreeReaderArray<UInt_t>(*fReader, Form("%s.subdetectorHitNumbers_end", branch));
	_TracksFromGenParticles_trackStates_location = new TTreeReaderArray<Int_t>(*fReader, Form("_%s_trackStates.location", branch));
	_TracksFromGenParticles_trackStates_omega = new TTreeReaderArray<Float_t>(*fReader, Form("_%s_trackStates.omega", branch));
	_TracksFromGenParticles_trackStates_phi = new TTreeReaderArray<Float_t>(*fReader, Form("_%s_trackStates.phi", branch));
	_TracksFromGenParticles_trackStates_tanLambda = new TTreeReaderArray<Float_t>(*fReader, Form("_%s_trackStates.tanLambda", branch));
	_TracksFromGenParticles_trackStates_referencePoint_x = new TTreeReaderArray<Float_t>(*fReader, Form("_%s_trackStates.referencePoint.x", branch));
	_TracksFromGenParticles_trackStates_referencePoint_y = new TTreeReaderArray<Float_t>(*fReader, Form("_%s_trackStates.referencePoint.y", branch));
	_TracksFromGenParticles_trackStates_referencePoint_z = new TTreeReaderArray<Float_t>(*fReader, Form("_%s_trackStates.referencePoint.z", branch));
	_TracksFromGenParticles_subdetectorHitNumbers = new TTreeReaderArray<UInt_t>(*fReader, Form("_%s_subdetectorHitNumbers", branch));
      }
    }
  }

  // hits in vertex detector
  if (displayConfig.getBoolConfig("drawVertexHits"))
  {
    // check that branch name for inner barrel is set and branch exists
    std::string branchName = displayConfig.getStringConfig("vertexBarrelHits");
    if (branchName == "") {
      std::cout << "WARNING: vertexBarrelHits not set, disabling vertex hits" << std::endl;
      displayConfig.setBoolConfig("drawVertexHits", false);
      displayConfig.setStringConfig("vertexEndcapHits", "");
    }
    else {
      const char* branch = branchName.c_str();
      if (! fReader->GetTree()->FindBranch(Form("%s.cellID", branch)))
      {
        std::cout << "WARNING: branch " << branch << ".cellID not found, disabling vertex hits" << std::endl;
        displayConfig.setBoolConfig("drawVertexHits", false);
        displayConfig.setStringConfig("vertexBarrelHits", "");
        displayConfig.setStringConfig("vertexEndcapHits", "");
      }
    }

    // check that branch name for endcap is set and branch exists
    branchName = displayConfig.getStringConfig("vertexEndcapHits");
    if (branchName == "") {
      std::cout << "WARNING: vertexEndcapHits not set, disabling vertex hits" << std::endl;
      displayConfig.setBoolConfig("drawVertexHits", false);
      displayConfig.setStringConfig("vertexBarrelHits", "");
    }
    else {
      const char* branch = branchName.c_str();
      if (! fReader->GetTree()->FindBranch(Form("%s.cellID", branch)))
      {
        std::cout << "WARNING: branch " << branch << ".cellID not found, disabling vertex hits" << std::endl;
        displayConfig.setBoolConfig("drawVertexHits", false);
        displayConfig.setStringConfig("vertexBarrelHits", "");
        displayConfig.setStringConfig("vertexEndcapHits", "");
      }
    }

    // read the branches
    if (displayConfig.getBoolConfig("drawVertexHits")) {
      branchName = displayConfig.getStringConfig("vertexBarrelHits");
      const char* branch = branchName.c_str();
      VertexBarrelHits_cellID     = new TTreeReaderArray<ULong_t>(*fReader, Form("%s.cellID", branch));
      VertexBarrelHits_energy     = new TTreeReaderArray<Float_t>(*fReader, Form("%s.eDep", branch));
      VertexBarrelHits_position_x = new TTreeReaderArray<Double_t>(*fReader, Form("%s.position.x", branch));
      VertexBarrelHits_position_y = new TTreeReaderArray<Double_t>(*fReader, Form("%s.position.y", branch));
      VertexBarrelHits_position_z = new TTreeReaderArray<Double_t>(*fReader, Form("%s.position.z", branch));

      branchName = displayConfig.getStringConfig("vertexEndcapHits");
      branch = branchName.c_str();
      VertexEndcapHits_cellID     = new TTreeReaderArray<ULong_t>(*fReader, Form("%s.cellID", branch));
      VertexEndcapHits_energy     = new TTreeReaderArray<Float_t>(*fReader, Form("%s.eDep", branch));
      VertexEndcapHits_position_x = new TTreeReaderArray<Double_t>(*fReader, Form("%s.position.x", branch));
      VertexEndcapHits_position_y = new TTreeReaderArray<Double_t>(*fReader, Form("%s.position.y", branch));
      VertexEndcapHits_position_z = new TTreeReaderArray<Double_t>(*fReader, Form("%s.position.z", branch));
    }
  }

  // hits in drift chamber
  if (displayConfig.getBoolConfig("drawDriftChamberHits"))
  {
    std::string branchName = displayConfig.getStringConfig("driftChamberHits");
    const char* branch = branchName.c_str();
    if (branchName == "") {
      std::cout << "WARNING: driftChamberHits not set, disabling drift chamber hits" << std::endl;
      displayConfig.setBoolConfig("drawDriftChamberHits", false);
    }
    else {
      if (! fReader->GetTree()->FindBranch(Form("%s.cellID", branch)))
      {
        std::cout << "WARNING: branch " << branch << ".cellID not found, disabling drift chamber hits" << std::endl;
        displayConfig.setBoolConfig("drawDriftChamberHits", false);
        displayConfig.setStringConfig("driftChamberHits", "");
      }
      else
      {
        DriftChamberHits_cellID     = new TTreeReaderArray<ULong_t>(*fReader, Form("%s.cellID", branch));
        DriftChamberHits_energy     = new TTreeReaderArray<Float_t>(*fReader, Form("%s.eDep", branch));
        DriftChamberHits_position_x = new TTreeReaderArray<Double_t>(*fReader, Form("%s.position.x", branch));
        DriftChamberHits_position_y = new TTreeReaderArray<Double_t>(*fReader, Form("%s.position.y", branch));
        DriftChamberHits_position_z = new TTreeReaderArray<Double_t>(*fReader, Form("%s.position.z", branch));
      }
    }
  }

  // hits in si wrapper
  if (displayConfig.getBoolConfig("drawSiWrapperHits"))
  {
    // check that branch name for inner barrel is set and branch exists
    std::string branchName = displayConfig.getStringConfig("siWrapperBarrelHits");
    if (branchName == "") {
      std::cout << "WARNING: siWrapperBarrelHits not set, disabling silicon wrapper hits" << std::endl;
      displayConfig.setBoolConfig("drawSiWrapperHits", false);
      displayConfig.setStringConfig("siWrapperEndcapHits", "");
    }
    else {
      const char* branch = branchName.c_str();
      if (! fReader->GetTree()->FindBranch(Form("%s.cellID", branch)))
      {
        std::cout << "WARNING: branch " << branch << ".cellID not found, disabling silicon wrapper hits" << std::endl;
        displayConfig.setBoolConfig("drawSiWrapperHits", false);
        displayConfig.setStringConfig("siWrapperBarrelHits", "");
        displayConfig.setStringConfig("siWrapperEndcapHits", "");
      }
    }

    // check that branch name for endcap is set and branch exists
    branchName = displayConfig.getStringConfig("siWrapperEndcapHits");
    if (branchName == "") {
      std::cout << "WARNING: siWrapperEndcapHits not set, disabling silicon wrapper hits" << std::endl;
      displayConfig.setBoolConfig("drawSiWrapperHits", false);
      displayConfig.setStringConfig("siWrapperBarrelHits", "");
    }
    else {
      const char* branch = branchName.c_str();
      if (! fReader->GetTree()->FindBranch(Form("%s.cellID", branch)))
      {
        std::cout << "WARNING: branch " << branch << ".cellID not found, disabling silicon wrapper hits" << std::endl;
        displayConfig.setBoolConfig("drawSiWrapperHits", false);
        displayConfig.setStringConfig("siWrapperBarrelHits", "");
        displayConfig.setStringConfig("siWrapperEndcapHits", "");
      }
    }

    // read the branches
    if (displayConfig.getBoolConfig("drawSiWrapperHits")) {
      branchName = displayConfig.getStringConfig("siWrapperBarrelHits");
      const char* branch = branchName.c_str();
      SiWrapperBarrelHits_cellID     = new TTreeReaderArray<ULong_t>(*fReader, Form("%s.cellID", branch));
      SiWrapperBarrelHits_energy     = new TTreeReaderArray<Float_t>(*fReader, Form("%s.eDep", branch));
      SiWrapperBarrelHits_position_x = new TTreeReaderArray<Double_t>(*fReader, Form("%s.position.x", branch));
      SiWrapperBarrelHits_position_y = new TTreeReaderArray<Double_t>(*fReader, Form("%s.position.y", branch));
      SiWrapperBarrelHits_position_z = new TTreeReaderArray<Double_t>(*fReader, Form("%s.position.z", branch));

      branchName = displayConfig.getStringConfig("siWrapperEndcapHits");
      branch = branchName.c_str();
      SiWrapperEndcapHits_cellID     = new TTreeReaderArray<ULong_t>(*fReader, Form("%s.cellID", branch));
      SiWrapperEndcapHits_energy     = new TTreeReaderArray<Float_t>(*fReader, Form("%s.eDep", branch));
      SiWrapperEndcapHits_position_x = new TTreeReaderArray<Double_t>(*fReader, Form("%s.position.x", branch));
      SiWrapperEndcapHits_position_y = new TTreeReaderArray<Double_t>(*fReader, Form("%s.position.y", branch));
      SiWrapperEndcapHits_position_z = new TTreeReaderArray<Double_t>(*fReader, Form("%s.position.z", branch));
    }
  }

  // digis in vertex detector
  if (displayConfig.getBoolConfig("drawVertexDigis"))
  {
    // check that branch name for inner barrel is set and branch exists
    std::string branchName = displayConfig.getStringConfig("vertexBarrelDigis");
    if (branchName == "") {
      std::cout << "WARNING: vertexBarrelDigis not set, disabling vertex digis" << std::endl;
      displayConfig.setBoolConfig("drawVertexDigis", false);
      displayConfig.setStringConfig("vertexEndcapDigis", "");
    }
    else {
      const char* branch = branchName.c_str();
      if (! fReader->GetTree()->FindBranch(Form("%s.cellID", branch)))
      {
        std::cout << "WARNING: branch " << branch << ".cellID not found, disabling vertex digis" << std::endl;
        displayConfig.setBoolConfig("drawVertexDigis", false);
        displayConfig.setStringConfig("vertexBarrelDigis", "");
        displayConfig.setStringConfig("vertexEndcapDigis", "");
      }
    }

    // check that branch name for endcap is set and branch exists
    branchName = displayConfig.getStringConfig("vertexEndcapDigis");
    if (branchName == "") {
      std::cout << "WARNING: vertexEndcapDigis not set, disabling vertex digis" << std::endl;
      displayConfig.setBoolConfig("drawVertexDigis", false);
      displayConfig.setStringConfig("vertexBarrelDigis", "");
    }
    else {
      const char* branch = branchName.c_str();
      if (! fReader->GetTree()->FindBranch(Form("%s.cellID", branch)))
      {
        std::cout << "WARNING: branch " << branch << ".cellID not found, disabling vertex digis" << std::endl;
        displayConfig.setBoolConfig("drawVertexDigis", false);
        displayConfig.setStringConfig("vertexBarrelDigis", "");
        displayConfig.setStringConfig("vertexEndcapDigis", "");
      }
    }

    // read the branches
    if (displayConfig.getBoolConfig("drawVertexDigis")) {
      branchName = displayConfig.getStringConfig("vertexBarrelDigis");
      const char* branch = branchName.c_str();
      VertexBarrelDigis_cellID     = new TTreeReaderArray<ULong_t>(*fReader, Form("%s.cellID", branch));
      VertexBarrelDigis_energy     = new TTreeReaderArray<Float_t>(*fReader, Form("%s.eDep", branch));
      VertexBarrelDigis_position_x = new TTreeReaderArray<Double_t>(*fReader, Form("%s.position.x", branch));
      VertexBarrelDigis_position_y = new TTreeReaderArray<Double_t>(*fReader, Form("%s.position.y", branch));
      VertexBarrelDigis_position_z = new TTreeReaderArray<Double_t>(*fReader, Form("%s.position.z", branch));

      branchName = displayConfig.getStringConfig("vertexEndcapDigis");
      branch = branchName.c_str();
      VertexEndcapDigis_cellID     = new TTreeReaderArray<ULong_t>(*fReader, Form("%s.cellID", branch));
      VertexEndcapDigis_energy     = new TTreeReaderArray<Float_t>(*fReader, Form("%s.eDep", branch));
      VertexEndcapDigis_position_x = new TTreeReaderArray<Double_t>(*fReader, Form("%s.position.x", branch));
      VertexEndcapDigis_position_y = new TTreeReaderArray<Double_t>(*fReader, Form("%s.position.y", branch));
      VertexEndcapDigis_position_z = new TTreeReaderArray<Double_t>(*fReader, Form("%s.position.z", branch));
    }
  }

  // digis in drift chamber
  if (displayConfig.getBoolConfig("drawDriftChamberDigis"))
  {
    std::string branchName = displayConfig.getStringConfig("driftChamberDigis");
    const char* branch = branchName.c_str();
    if (branchName == "") {
      std::cout << "WARNING: driftChamberDigis not set, disabling drift chamber digis" << std::endl;
      displayConfig.setBoolConfig("drawDriftChamberDigis", false);
    }
    else {
      if (! fReader->GetTree()->FindBranch(Form("%s.cellID", branch)))
      {
        std::cout << "WARNING: branch " << branch << ".cellID not found, disabling drift chamber digis" << std::endl;
        displayConfig.setBoolConfig("drawDriftChamberDigis", false);
        displayConfig.setStringConfig("driftChamberDigis", "");
      }
      else
      {
        DriftChamberDigis_cellID     = new TTreeReaderArray<ULong_t>(*fReader, Form("%s.cellID", branch));
        DriftChamberDigis_energy     = new TTreeReaderArray<Float_t>(*fReader, Form("%s.eDep", branch));
        DriftChamberDigis_position_x = new TTreeReaderArray<Double_t>(*fReader, Form("%s.position.x", branch));
        DriftChamberDigis_position_y = new TTreeReaderArray<Double_t>(*fReader, Form("%s.position.y", branch));
        DriftChamberDigis_position_z = new TTreeReaderArray<Double_t>(*fReader, Form("%s.position.z", branch));
      }
    }
  }

  // digis in si wrapper
  if (displayConfig.getBoolConfig("drawSiWrapperDigis"))
  {
    // check that branch name for inner barrel is set and branch exists
    std::string branchName = displayConfig.getStringConfig("siWrapperBarrelDigis");
    if (branchName == "") {
      std::cout << "WARNING: siWrapperBarrelDigis not set, disabling silicon wrapper digis" << std::endl;
      displayConfig.setBoolConfig("drawSiWrapperDigis", false);
      displayConfig.setStringConfig("siWrapperEndcapDigis", "");
    }
    else {
      const char* branch = branchName.c_str();
      if (! fReader->GetTree()->FindBranch(Form("%s.cellID", branch)))
      {
        std::cout << "WARNING: branch " << branch << ".cellID not found, disabling silicon wrapper digis" << std::endl;
        displayConfig.setBoolConfig("drawSiWrapperDigis", false);
        displayConfig.setStringConfig("siWrapperBarrelDigis", "");
        displayConfig.setStringConfig("siWrapperEndcapDigis", "");
      }
    }

    // check that branch name for endcap is set and branch exists
    branchName = displayConfig.getStringConfig("siWrapperEndcapDigis");
    if (branchName == "") {
      std::cout << "WARNING: siWrapperEndcapDigis not set, disabling silicon wrapper digis" << std::endl;
      displayConfig.setBoolConfig("drawSiWrapperDigis", false);
      displayConfig.setStringConfig("siWrapperBarrelDigis", "");
    }
    else {
      const char* branch = branchName.c_str();
      if (! fReader->GetTree()->FindBranch(Form("%s.cellID", branch)))
      {
        std::cout << "WARNING: branch " << branch << ".cellID not found, disabling silicon wrapper digis" << std::endl;
        displayConfig.setBoolConfig("drawSiWrapperDigis", false);
        displayConfig.setStringConfig("siWrapperBarrelDigis", "");
        displayConfig.setStringConfig("siWrapperEndcapDigis", "");
      }
    }

    // read the branches
    if (displayConfig.getBoolConfig("drawSiWrapperDigis")) {
      branchName = displayConfig.getStringConfig("siWrapperBarrelDigis");
      const char* branch = branchName.c_str();
      SiWrapperBarrelDigis_cellID     = new TTreeReaderArray<ULong_t>(*fReader, Form("%s.cellID", branch));
      SiWrapperBarrelDigis_energy     = new TTreeReaderArray<Float_t>(*fReader, Form("%s.eDep", branch));
      SiWrapperBarrelDigis_position_x = new TTreeReaderArray<Double_t>(*fReader, Form("%s.position.x", branch));
      SiWrapperBarrelDigis_position_y = new TTreeReaderArray<Double_t>(*fReader, Form("%s.position.y", branch));
      SiWrapperBarrelDigis_position_z = new TTreeReaderArray<Double_t>(*fReader, Form("%s.position.z", branch));

      branchName = displayConfig.getStringConfig("siWrapperEndcapDigis");
      branch = branchName.c_str();
      SiWrapperEndcapDigis_cellID     = new TTreeReaderArray<ULong_t>(*fReader, Form("%s.cellID", branch));
      SiWrapperEndcapDigis_energy     = new TTreeReaderArray<Float_t>(*fReader, Form("%s.eDep", branch));
      SiWrapperEndcapDigis_position_x = new TTreeReaderArray<Double_t>(*fReader, Form("%s.position.x", branch));
      SiWrapperEndcapDigis_position_y = new TTreeReaderArray<Double_t>(*fReader, Form("%s.position.y", branch));
      SiWrapperEndcapDigis_position_z = new TTreeReaderArray<Double_t>(*fReader, Form("%s.position.z", branch));
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
      if (! fReader->GetTree()->FindBranch(Form("%s.energy", branch)))
      {
        std::cout << "WARNING: branch " << branch << ".energy not found, disabling ecal barrel hits" << std::endl;
        displayConfig.setBoolConfig("drawECalBarrelHits", false);
        displayConfig.setStringConfig("ecalBarrelHits", "");
      }
      else
      {
        ECalBarrelHits_energy     = new TTreeReaderArray<Float_t>(*fReader, Form("%s.energy", branch));
        ECalBarrelHits_position_x = new TTreeReaderArray<Float_t>(*fReader, Form("%s.stepPosition.x", branch));
        ECalBarrelHits_position_y = new TTreeReaderArray<Float_t>(*fReader, Form("%s.stepPosition.y", branch));
        ECalBarrelHits_position_z = new TTreeReaderArray<Float_t>(*fReader, Form("%s.stepPosition.z", branch));
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
        ECalBarrelCells_cellID     = new TTreeReaderArray<ULong_t>(*fReader, Form("%s.cellID", branch));
        ECalBarrelCells_energy     = new TTreeReaderArray<Float_t>(*fReader, Form("%s.energy", branch));
        ECalBarrelCells_position_x = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.x", branch));
        ECalBarrelCells_position_y = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.y", branch));
        ECalBarrelCells_position_z = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.z", branch));
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
        // ECalBarrelCells2_cellID     = new TTreeReaderArray<ULong_t>(*fReader, Form("%s.cellID", branch));    
        ECalBarrelCells2_energy     = new TTreeReaderArray<Float_t>(*fReader, Form("%s.energy", branch));    
        ECalBarrelCells2_position_x = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.x", branch));
        ECalBarrelCells2_position_y = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.y", branch));
        ECalBarrelCells2_position_z = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.z", branch));
      }
    }
  }
  else {
    displayConfig.setStringConfig("ecalBarrelMergedCells", "");
  }

  // hits in ECal endcap
  if (displayConfig.getBoolConfig("drawECalEndcapHits"))
  {
    std::string branchName = displayConfig.getStringConfig("ecalEndcapHits");
    const char* branch = branchName.c_str();
    if (branchName == "") {
      std::cout << "WARNING: ecalEndcapHits not set, disabling ecal endcap hits" << std::endl;
      displayConfig.setBoolConfig("drawECalEndcapHits", false);
    }
    else {
      if (! fReader->GetTree()->FindBranch(Form("%s.energy", branch)))
      {
        std::cout << "WARNING: branch " << branch << ".energy not found, disabling ecal endcap hits" << std::endl;
        displayConfig.setBoolConfig("drawECalEndcapHits", false);
        displayConfig.setStringConfig("ecalEndcapHits", "");
      }
      else
      {
        ECalEndcapHits_energy     = new TTreeReaderArray<Float_t>(*fReader, Form("%s.energy", branch));
        ECalEndcapHits_position_x = new TTreeReaderArray<Float_t>(*fReader, Form("%s.stepPosition.x", branch));
        ECalEndcapHits_position_y = new TTreeReaderArray<Float_t>(*fReader, Form("%s.stepPosition.y", branch));
        ECalEndcapHits_position_z = new TTreeReaderArray<Float_t>(*fReader, Form("%s.stepPosition.z", branch));
      }
    }
  }

  // cells in ECal endcap
  if (displayConfig.getBoolConfig("drawECalEndcapCells"))
  {
    std::string branchName = displayConfig.getStringConfig("ecalEndcapCells");
    const char* branch = branchName.c_str();
    if (branchName == "") {
      std::cout << "WARNING: ecalEndcapCells not set, disabling ecal endcap cells" << std::endl;
      displayConfig.setBoolConfig("drawECalEndcapCells", false);
    }
    else {
      if (! fReader->GetTree()->FindBranch(Form("%s.cellID", branch)))
      {
        std::cout << "WARNING: branch " << branch << ".cellID not found, disabling ecal endcap cells" << std::endl;
        displayConfig.setBoolConfig("drawECalEndcapCells", false);
        displayConfig.setStringConfig("ecalEndcapCells", "");
      }
      else
      {
        ECalEndcapCells_cellID     = new TTreeReaderArray<ULong_t>(*fReader, Form("%s.cellID", branch));
        ECalEndcapCells_energy     = new TTreeReaderArray<Float_t>(*fReader, Form("%s.energy", branch));
        ECalEndcapCells_position_x = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.x", branch));
        ECalEndcapCells_position_y = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.y", branch));
        ECalEndcapCells_position_z = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.z", branch));
      }
    }
  }
  
  if (m_doHCal) {
    // hits in HCal barrel
    if (displayConfig.getBoolConfig("drawHCalBarrelHits"))
    {
      std::string branchName = displayConfig.getStringConfig("hcalBarrelHits");
      const char* branch = branchName.c_str();
      if (branchName == "") {
        std::cout << "WARNING: hcalBarrelHits not set, disabling hcal barrel hits" << std::endl;
        displayConfig.setBoolConfig("drawHCalBarrelHits", false);
      }
      else {
        if (! fReader->GetTree()->FindBranch(Form("%s.energy", branch)))
        {
          std::cout << "WARNING: branch " << branch << ".energy not found, disabling hcal barrel hits" << std::endl;
          displayConfig.setBoolConfig("drawHCalBarrelHits", false);
          displayConfig.setStringConfig("hcalBarrelHits", "");
        }
        else
        {
          HCalBarrelHits_energy     = new TTreeReaderArray<Float_t>(*fReader, Form("%s.energy", branch));    
          HCalBarrelHits_position_x = new TTreeReaderArray<Float_t>(*fReader, Form("%s.stepPosition.x", branch));
          HCalBarrelHits_position_y = new TTreeReaderArray<Float_t>(*fReader, Form("%s.stepPosition.y", branch));
          HCalBarrelHits_position_z = new TTreeReaderArray<Float_t>(*fReader, Form("%s.stepPosition.z", branch));
        }
      }
    }
    
    // cells in HCal barrel
    if (displayConfig.getBoolConfig("drawHCalBarrelCells"))
    {
      std::string branchName = displayConfig.getStringConfig("hcalBarrelCells");
      const char* branch = branchName.c_str();
      if (branchName == "") {
       std::cout << "WARNING: hcalBarrelCells not set, disabling hcal barrel cells" << std::endl;
        displayConfig.setBoolConfig("drawHCalBarrelCells", false);
      }
      else {
        if (! fReader->GetTree()->FindBranch(Form("%s.cellID", branch)))
        {
          std::cout << "WARNING: branch " << branch << ".cellID not found, disabling hcal barrel cells" << std::endl;
          displayConfig.setBoolConfig("drawHCalBarrelCells", false);
          displayConfig.setStringConfig("hcalBarrelCells", "");
        }
        else
        {
          HCalBarrelCells_cellID     = new TTreeReaderArray<ULong_t>(*fReader, Form("%s.cellID", branch));
          HCalBarrelCells_energy     = new TTreeReaderArray<Float_t>(*fReader, Form("%s.energy", branch));
          HCalBarrelCells_position_x = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.x", branch));
          HCalBarrelCells_position_y = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.y", branch));
          HCalBarrelCells_position_z = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.z", branch));
        }
      }
    }

    // hits in HCal endcap
    if (displayConfig.getBoolConfig("drawHCalEndcapHits"))
    {
      std::string branchName = displayConfig.getStringConfig("hcalEndcapHits");
      const char* branch = branchName.c_str();
      if (branchName == "") {
        std::cout << "WARNING: hcalEndcapHits not set, disabling hcal endcap hits" << std::endl;
        displayConfig.setBoolConfig("drawHCalEndcapHits", false);
      }
      else {
        if (! fReader->GetTree()->FindBranch(Form("%s.energy", branch)))
        {
          std::cout << "WARNING: branch " << branch << ".energy not found, disabling hcal endcap hits" << std::endl;
          displayConfig.setBoolConfig("drawHCalEndcapHits", false);
          displayConfig.setStringConfig("hcalEndcapHits", "");
        }
        else
        {
          HCalEndcapHits_energy     = new TTreeReaderArray<Float_t>(*fReader, Form("%s.energy", branch));    
          HCalEndcapHits_position_x = new TTreeReaderArray<Float_t>(*fReader, Form("%s.stepPosition.x", branch));
          HCalEndcapHits_position_y = new TTreeReaderArray<Float_t>(*fReader, Form("%s.stepPosition.y", branch));
          HCalEndcapHits_position_z = new TTreeReaderArray<Float_t>(*fReader, Form("%s.stepPosition.z", branch));
        }
      }
    }
    
    // cells in HCal endcap
    if (displayConfig.getBoolConfig("drawHCalEndcapCells"))
    {
      std::string branchName = displayConfig.getStringConfig("hcalEndcapCells");
      const char* branch = branchName.c_str();
      if (branchName == "") {
       std::cout << "WARNING: hcalEndcapCells not set, disabling hcal endcap cells" << std::endl;
        displayConfig.setBoolConfig("drawHCalEndcapCells", false);
      }
      else {
        if (! fReader->GetTree()->FindBranch(Form("%s.cellID", branch)))
        {
          std::cout << "WARNING: branch " << branch << ".cellID not found, disabling hcal endcap cells" << std::endl;
          displayConfig.setBoolConfig("drawHCalEndcapCells", false);
          displayConfig.setStringConfig("hcalEndcapCells", "");
        }
        else
        {
          HCalEndcapCells_cellID     = new TTreeReaderArray<ULong_t>(*fReader, Form("%s.cellID", branch));
          HCalEndcapCells_energy     = new TTreeReaderArray<Float_t>(*fReader, Form("%s.energy", branch));
          HCalEndcapCells_position_x = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.x", branch));
          HCalEndcapCells_position_y = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.y", branch));
          HCalEndcapCells_position_z = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.z", branch));
        }
      }
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
  if (displayConfig.getBoolConfig("drawMuonHits"))
  {
    // check that branch name for inner barrel is set and branch exists
    std::string branchName = displayConfig.getStringConfig("muonBarrelHits");
    if (branchName == "") {
      std::cout << "WARNING: muonBarrelHits not set, disabling muon hits" << std::endl;
      displayConfig.setBoolConfig("drawMuonHits", false);
      displayConfig.setStringConfig("muonEndcapHits", "");
    }
    else {
      const char* branch = branchName.c_str();
      // tracker-based
      // if (! fReader->GetTree()->FindBranch(Form("%s.cellID", branch)))
      // calo-based
      if (! fReader->GetTree()->FindBranch(Form("%s.energy", branch)))
      {
        std::cout << "WARNING: branch " << branch << ".cellID not found, disabling muon hits" << std::endl;
        displayConfig.setBoolConfig("drawMuonHits", false);
        displayConfig.setStringConfig("muonBarrelHits", "");
        displayConfig.setStringConfig("muonEndcapHits", "");
      }
    }

    // check that branch name for endcap is set and branch exists
    branchName = displayConfig.getStringConfig("muonEndcapHits");
    if (branchName == "") {
      std::cout << "WARNING: muonEndcapHits not set, disabling muon hits" << std::endl;
      displayConfig.setBoolConfig("drawMuonHits", false);
      displayConfig.setStringConfig("muonBarrelHits", "");
    }
    else {
      const char* branch = branchName.c_str();
      if (! fReader->GetTree()->FindBranch(Form("%s.cellID", branch)))
      {
        std::cout << "WARNING: branch " << branch << ".cellID not found, disabling muon hits" << std::endl;
        displayConfig.setBoolConfig("drawMuonHits", false);
        displayConfig.setStringConfig("muonBarrelHits", "");
        displayConfig.setStringConfig("muonEndcapHits", "");
      }
    }

    // read the branches
    if (displayConfig.getBoolConfig("drawMuonHits")) {
      branchName = displayConfig.getStringConfig("muonBarrelHits");
      const char* branch = branchName.c_str();
      // for tracker-based SD
      // MuonBarrelHits_cellID     = new TTreeReaderArray<ULong_t>(*fReader, Form("%s.cellID", branch));
      // MuonBarrelHits_energy     = new TTreeReaderArray<Float_t>(*fReader, Form("%s.eDep", branch));
      // MuonBarrelHits_position_x = new TTreeReaderArray<Double_t>(*fReader, Form("%s.position.x", branch));
      // MuonBarrelHits_position_y = new TTreeReaderArray<Double_t>(*fReader, Form("%s.position.y", branch));
      // MuonBarrelHits_position_z = new TTreeReaderArray<Double_t>(*fReader, Form("%s.position.z", branch));
      // for calo-based SD
      MuonBarrelHits_energy     = new TTreeReaderArray<Float_t>(*fReader, Form("%s.energy", branch));
      MuonBarrelHits_position_x = new TTreeReaderArray<Float_t>(*fReader, Form("%s.stepPosition.x", branch));
      MuonBarrelHits_position_y = new TTreeReaderArray<Float_t>(*fReader, Form("%s.stepPosition.y", branch));
      MuonBarrelHits_position_z = new TTreeReaderArray<Float_t>(*fReader, Form("%s.stepPosition.z", branch));

      branchName = displayConfig.getStringConfig("muonEndcapHits");
      branch = branchName.c_str();
      MuonEndcapHits_cellID     = new TTreeReaderArray<ULong_t>(*fReader, Form("%s.cellID", branch));
      // for tracker-based SD
      MuonEndcapHits_energy     = new TTreeReaderArray<Float_t>(*fReader, Form("%s.eDep", branch));
      MuonEndcapHits_position_x = new TTreeReaderArray<Double_t>(*fReader, Form("%s.position.x", branch));
      MuonEndcapHits_position_y = new TTreeReaderArray<Double_t>(*fReader, Form("%s.position.y", branch));
      MuonEndcapHits_position_z = new TTreeReaderArray<Double_t>(*fReader, Form("%s.position.z", branch));
      // for calo-based SD
      // MuonEndcapHits_energy     = new TTreeReaderArray<Float_t>(*fReader, Form("%s.energy", branch));
      // MuonEndcapHits_position_x = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.x", branch));
      // MuonEndcapHits_position_y = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.y", branch));
      // MuonEndcapHits_position_z = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.z", branch));
    }
  }

  // cells in muon tagger
  if (displayConfig.getBoolConfig("drawMuonCells"))
  {
    // check that branch name for inner barrel is set and branch exists
    std::string branchName = displayConfig.getStringConfig("muonBarrelCells");
    if (branchName == "") {
      std::cout << "WARNING: muonBarrelCells not set, disabling muon cells" << std::endl;
      displayConfig.setBoolConfig("drawMuonCells", false);
      displayConfig.setStringConfig("muonEndcapCells", "");
    }
    else {
      const char* branch = branchName.c_str();
      if (! fReader->GetTree()->FindBranch(Form("%s.cellID", branch)))
      {
        std::cout << "WARNING: branch " << branch << ".cellID not found, disabling muon cells" << std::endl;
        displayConfig.setBoolConfig("drawMuonCells", false);
        displayConfig.setStringConfig("muonBarrelCells", "");
        displayConfig.setStringConfig("muonEndcapCells", "");
      }
    }

    // check that branch name for endcap is set and branch exists
    branchName = displayConfig.getStringConfig("muonEndcapCells");
    if (branchName == "") {
      std::cout << "WARNING: muonEndcapCells not set, disabling muon cells" << std::endl;
      displayConfig.setBoolConfig("drawMuonCells", false);
      displayConfig.setStringConfig("muonBarrelCells", "");
    }
    else {
      const char* branch = branchName.c_str();
      if (! fReader->GetTree()->FindBranch(Form("%s.cellID", branch)))
      {
        std::cout << "WARNING: branch " << branch << ".cellID not found, disabling muon cells" << std::endl;
        displayConfig.setBoolConfig("drawMuonCells", false);
        displayConfig.setStringConfig("muonBarrelCells", "");
        displayConfig.setStringConfig("muonEndcapCells", "");
      }
    }

    // read the branches
    if (displayConfig.getBoolConfig("drawMuonCells")) {
      branchName = displayConfig.getStringConfig("muonBarrelCells");
      const char* branch = branchName.c_str();
      MuonBarrelCells_cellID     = new TTreeReaderArray<ULong_t>(*fReader, Form("%s.cellID", branch));
      MuonBarrelCells_energy     = new TTreeReaderArray<Float_t>(*fReader, Form("%s.energy", branch));
      MuonBarrelCells_position_x = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.x", branch));
      MuonBarrelCells_position_y = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.y", branch));
      MuonBarrelCells_position_z = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.z", branch));

      branchName = displayConfig.getStringConfig("muonEndcapCells");
      branch = branchName.c_str();
      MuonEndcapCells_cellID     = new TTreeReaderArray<ULong_t>(*fReader, Form("%s.cellID", branch));
      MuonEndcapCells_energy     = new TTreeReaderArray<Float_t>(*fReader, Form("%s.energy", branch));
      MuonEndcapCells_position_x = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.x", branch));
      MuonEndcapCells_position_y = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.y", branch));
      MuonEndcapCells_position_z = new TTreeReaderArray<Float_t>(*fReader, Form("%s.position.z", branch));
    }
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
    delete genParticles_endpoint_x;
    delete genParticles_endpoint_y;
    delete genParticles_endpoint_z;
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
  if (TracksFromGenParticles_trackStates_begin) {
    delete TracksFromGenParticles_trackStates_begin;
    delete TracksFromGenParticles_trackStates_end;
    delete TracksFromGenParticles_subdetectorHitNumbers_begin;
    delete TracksFromGenParticles_subdetectorHitNumbers_end;
    delete _TracksFromGenParticles_trackStates_location;
    delete _TracksFromGenParticles_trackStates_omega;
    delete _TracksFromGenParticles_trackStates_phi;
    delete _TracksFromGenParticles_trackStates_tanLambda;
    delete _TracksFromGenParticles_trackStates_referencePoint_x;
    delete _TracksFromGenParticles_trackStates_referencePoint_y;
    delete _TracksFromGenParticles_trackStates_referencePoint_z;
    delete _TracksFromGenParticles_subdetectorHitNumbers;
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
  if (DriftChamberHits_cellID) {
    delete DriftChamberHits_cellID;
    delete DriftChamberHits_energy;
    delete DriftChamberHits_position_x;
    delete DriftChamberHits_position_y;
    delete DriftChamberHits_position_z;
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
  if (DriftChamberDigis_cellID) {
    delete DriftChamberDigis_cellID;
    delete DriftChamberDigis_energy;
    delete DriftChamberDigis_position_x;
    delete DriftChamberDigis_position_y;
    delete DriftChamberDigis_position_z;
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
