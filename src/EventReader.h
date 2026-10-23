/******************************************************************************/
// Simple event display for the ALLEGRO detector with ECAL with inclined modules
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
#include <TTreeReaderValue.h>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

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
  TTreeReaderArray<Double_t> *genParticles_endpoint_x = nullptr;
  TTreeReaderArray<Double_t> *genParticles_endpoint_y = nullptr;
  TTreeReaderArray<Double_t> *genParticles_endpoint_z = nullptr;
  TTreeReaderArray<Double_t> *genParticles_momentum_x = nullptr;
  TTreeReaderArray<Double_t> *genParticles_momentum_y = nullptr;
  TTreeReaderArray<Double_t> *genParticles_momentum_z = nullptr;
  TTreeReaderArray<UInt_t> *genParticles_parents_begin = nullptr;
  TTreeReaderArray<UInt_t> *genParticles_parents_end = nullptr;
  TTreeReaderArray<Int_t> *genParticles_parents_index = nullptr;
  TTreeReaderArray<UInt_t> *genParticles_daughters_begin = nullptr;
  TTreeReaderArray<UInt_t> *genParticles_daughters_end = nullptr;
  TTreeReaderArray<Int_t> *genParticles_daughters_index = nullptr;

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
  // these variables are double in latest versions
  // TTreeReaderArray<Float_t> *SimParticleSecondaries_momentum_x = nullptr;
  // TTreeReaderArray<Float_t> *SimParticleSecondaries_momentum_y = nullptr;
  // TTreeReaderArray<Float_t> *SimParticleSecondaries_momentum_z = nullptr;
  TTreeReaderArray<Double_t> *SimParticleSecondaries_momentum_x = nullptr;
  TTreeReaderArray<Double_t> *SimParticleSecondaries_momentum_y = nullptr;
  TTreeReaderArray<Double_t> *SimParticleSecondaries_momentum_z = nullptr;

  // the tracks
  TTreeReaderArray<UInt_t> *Tracks_trackStates_begin = nullptr;
  TTreeReaderArray<UInt_t> *Tracks_trackStates_end = nullptr;
  TTreeReaderArray<UInt_t> *Tracks_subdetectorHitNumbers_begin = nullptr;
  TTreeReaderArray<UInt_t> *Tracks_subdetectorHitNumbers_end = nullptr;
  TTreeReaderArray<Int_t> *_Tracks_trackStates_location = nullptr;
  TTreeReaderArray<Float_t> *_Tracks_trackStates_omega = nullptr;
  TTreeReaderArray<Float_t> *_Tracks_trackStates_phi = nullptr;
  TTreeReaderArray<Float_t> *_Tracks_trackStates_tanLambda = nullptr;
  TTreeReaderArray<Float_t> *_Tracks_trackStates_referencePoint_x = nullptr;
  TTreeReaderArray<Float_t> *_Tracks_trackStates_referencePoint_y = nullptr;
  TTreeReaderArray<Float_t> *_Tracks_trackStates_referencePoint_z = nullptr;
  TTreeReaderArray<Int_t> *_Tracks_subdetectorHitNumbers = nullptr;
  
  // the hits in the vertex (barrel)
  TTreeReaderArray<ULong_t> *VertexBarrelHits_cellID = nullptr;
  TTreeReaderArray<Float_t> *VertexBarrelHits_energy = nullptr;
  TTreeReaderArray<Double_t> *VertexBarrelHits_position_x = nullptr;
  TTreeReaderArray<Double_t> *VertexBarrelHits_position_y = nullptr;
  TTreeReaderArray<Double_t> *VertexBarrelHits_position_z = nullptr;

  // the hits in the vertex (endcap)
  TTreeReaderArray<ULong_t> *VertexEndcapHits_cellID = nullptr;
  TTreeReaderArray<Float_t> *VertexEndcapHits_energy = nullptr;
  TTreeReaderArray<Double_t> *VertexEndcapHits_position_x = nullptr;
  TTreeReaderArray<Double_t> *VertexEndcapHits_position_y = nullptr;
  TTreeReaderArray<Double_t> *VertexEndcapHits_position_z = nullptr;

  // the hits in the DCH/STT
  TTreeReaderArray<ULong_t> *MainTrackerHits_cellID = nullptr;
  TTreeReaderArray<Float_t> *MainTrackerHits_energy = nullptr;
  TTreeReaderArray<Double_t> *MainTrackerHits_position_x = nullptr;
  TTreeReaderArray<Double_t> *MainTrackerHits_position_y = nullptr;
  TTreeReaderArray<Double_t> *MainTrackerHits_position_z = nullptr;

  // the hits in the Si wrapper (barrel)
  TTreeReaderArray<ULong_t> *SiWrapperBarrelHits_cellID = nullptr;
  TTreeReaderArray<Float_t> *SiWrapperBarrelHits_energy = nullptr;
  TTreeReaderArray<Double_t> *SiWrapperBarrelHits_position_x = nullptr;
  TTreeReaderArray<Double_t> *SiWrapperBarrelHits_position_y = nullptr;
  TTreeReaderArray<Double_t> *SiWrapperBarrelHits_position_z = nullptr;

  // the hits in the Si wrapper (endcap)
  TTreeReaderArray<ULong_t> *SiWrapperEndcapHits_cellID = nullptr;
  TTreeReaderArray<Float_t> *SiWrapperEndcapHits_energy = nullptr;
  TTreeReaderArray<Double_t> *SiWrapperEndcapHits_position_x = nullptr;
  TTreeReaderArray<Double_t> *SiWrapperEndcapHits_position_y = nullptr;
  TTreeReaderArray<Double_t> *SiWrapperEndcapHits_position_z = nullptr;

  // the digis in the vertex (barrel)
  TTreeReaderArray<ULong_t> *VertexBarrelDigis_cellID = nullptr;
  TTreeReaderArray<Float_t> *VertexBarrelDigis_energy = nullptr;
  TTreeReaderArray<Double_t> *VertexBarrelDigis_position_x = nullptr;
  TTreeReaderArray<Double_t> *VertexBarrelDigis_position_y = nullptr;
  TTreeReaderArray<Double_t> *VertexBarrelDigis_position_z = nullptr;

  // the digis in the vertex (endcap)
  TTreeReaderArray<ULong_t> *VertexEndcapDigis_cellID = nullptr;
  TTreeReaderArray<Float_t> *VertexEndcapDigis_energy = nullptr;
  TTreeReaderArray<Double_t> *VertexEndcapDigis_position_x = nullptr;
  TTreeReaderArray<Double_t> *VertexEndcapDigis_position_y = nullptr;
  TTreeReaderArray<Double_t> *VertexEndcapDigis_position_z = nullptr;

  // the digis in the DCH/STT
  TTreeReaderArray<ULong_t> *MainTrackerDigis_cellID = nullptr;
  TTreeReaderArray<Float_t> *MainTrackerDigis_energy = nullptr;
  TTreeReaderArray<Double_t> *MainTrackerDigis_position_x = nullptr;
  TTreeReaderArray<Double_t> *MainTrackerDigis_position_y = nullptr;
  TTreeReaderArray<Double_t> *MainTrackerDigis_position_z = nullptr;

  // the digis in the Si wrapper (barrel)
  TTreeReaderArray<ULong_t> *SiWrapperBarrelDigis_cellID = nullptr;
  TTreeReaderArray<Float_t> *SiWrapperBarrelDigis_energy = nullptr;
  TTreeReaderArray<Double_t> *SiWrapperBarrelDigis_position_x = nullptr;
  TTreeReaderArray<Double_t> *SiWrapperBarrelDigis_position_y = nullptr;
  TTreeReaderArray<Double_t> *SiWrapperBarrelDigis_position_z = nullptr;

  // the digis in the Si wrapper (endcap)
  TTreeReaderArray<ULong_t> *SiWrapperEndcapDigis_cellID = nullptr;
  TTreeReaderArray<Float_t> *SiWrapperEndcapDigis_energy = nullptr;
  TTreeReaderArray<Double_t> *SiWrapperEndcapDigis_position_x = nullptr;
  TTreeReaderArray<Double_t> *SiWrapperEndcapDigis_position_y = nullptr;
  TTreeReaderArray<Double_t> *SiWrapperEndcapDigis_position_z = nullptr;

  // the digis in ECal barrel
  TTreeReaderArray<Float_t> *ECalBarrelHits_energy = nullptr;
  TTreeReaderArray<Float_t> *ECalBarrelHits_position_x = nullptr;
  TTreeReaderArray<Float_t> *ECalBarrelHits_position_y = nullptr;
  TTreeReaderArray<Float_t> *ECalBarrelHits_position_z = nullptr;

  // the cells in ECal barrel
  TTreeReaderArray<ULong_t> *ECalBarrelCells_cellID = nullptr;
  TTreeReaderArray<Float_t> *ECalBarrelCells_energy = nullptr;
  TTreeReaderArray<Float_t> *ECalBarrelCells_position_x = nullptr;
  TTreeReaderArray<Float_t> *ECalBarrelCells_position_y = nullptr;
  TTreeReaderArray<Float_t> *ECalBarrelCells_position_z = nullptr;

  // the cells in ECal barrel with coarser merging
  TTreeReaderArray<ULong_t> *ECalBarrelCells2_cellID = nullptr;
  TTreeReaderArray<Float_t> *ECalBarrelCells2_energy = nullptr;
  TTreeReaderArray<Float_t> *ECalBarrelCells2_position_x = nullptr;
  TTreeReaderArray<Float_t> *ECalBarrelCells2_position_y = nullptr;
  TTreeReaderArray<Float_t> *ECalBarrelCells2_position_z = nullptr;

  // the hits in ECal endcap
  TTreeReaderArray<Float_t> *ECalEndcapHits_energy = nullptr;
  TTreeReaderArray<Float_t> *ECalEndcapHits_position_x = nullptr;
  TTreeReaderArray<Float_t> *ECalEndcapHits_position_y = nullptr;
  TTreeReaderArray<Float_t> *ECalEndcapHits_position_z = nullptr;

  // the cells in ECal endcap
  TTreeReaderArray<ULong_t> *ECalEndcapCells_cellID = nullptr;
  TTreeReaderArray<Float_t> *ECalEndcapCells_energy = nullptr;
  TTreeReaderArray<Float_t> *ECalEndcapCells_position_x = nullptr;
  TTreeReaderArray<Float_t> *ECalEndcapCells_position_y = nullptr;
  TTreeReaderArray<Float_t> *ECalEndcapCells_position_z = nullptr;

  // the hits in HCal barrel
  TTreeReaderArray<Float_t> *HCalBarrelHits_energy = nullptr;
  TTreeReaderArray<Float_t> *HCalBarrelHits_position_x = nullptr;
  TTreeReaderArray<Float_t> *HCalBarrelHits_position_y = nullptr;
  TTreeReaderArray<Float_t> *HCalBarrelHits_position_z = nullptr;

  // the cells in HCal barrel
  TTreeReaderArray<ULong_t> *HCalBarrelCells_cellID = nullptr;
  TTreeReaderArray<Float_t> *HCalBarrelCells_energy = nullptr;
  TTreeReaderArray<Float_t> *HCalBarrelCells_position_x = nullptr;
  TTreeReaderArray<Float_t> *HCalBarrelCells_position_y = nullptr;
  TTreeReaderArray<Float_t> *HCalBarrelCells_position_z = nullptr;

  // the hits in HCal endcap
  TTreeReaderArray<Float_t> *HCalEndcapHits_energy = nullptr;
  TTreeReaderArray<Float_t> *HCalEndcapHits_position_x = nullptr;
  TTreeReaderArray<Float_t> *HCalEndcapHits_position_y = nullptr;
  TTreeReaderArray<Float_t> *HCalEndcapHits_position_z = nullptr;

  // the cells in HCal endcap
  TTreeReaderArray<ULong_t> *HCalEndcapCells_cellID = nullptr;
  TTreeReaderArray<Float_t> *HCalEndcapCells_energy = nullptr;
  TTreeReaderArray<Float_t> *HCalEndcapCells_position_x = nullptr;
  TTreeReaderArray<Float_t> *HCalEndcapCells_position_y = nullptr;
  TTreeReaderArray<Float_t> *HCalEndcapCells_position_z = nullptr;

  // the hits in muon barrel
  TTreeReaderArray<ULong_t> *MuonBarrelHits_cellID = nullptr;
  TTreeReaderArray<Float_t> *MuonBarrelHits_energy = nullptr;
  // TTreeReaderArray<Double_t> *MuonBarrelHits_position_x = nullptr;
  // TTreeReaderArray<Double_t> *MuonBarrelHits_position_y = nullptr;
  // TTreeReaderArray<Double_t> *MuonBarrelHits_position_z = nullptr;
  TTreeReaderArray<Float_t> *MuonBarrelHits_position_x = nullptr;
  TTreeReaderArray<Float_t> *MuonBarrelHits_position_y = nullptr;
  TTreeReaderArray<Float_t> *MuonBarrelHits_position_z = nullptr;

  // the cells in muon barrel
  TTreeReaderArray<ULong_t> *MuonBarrelCells_cellID = nullptr;
  TTreeReaderArray<Float_t> *MuonBarrelCells_energy = nullptr;
  TTreeReaderArray<Float_t> *MuonBarrelCells_position_x = nullptr;
  TTreeReaderArray<Float_t> *MuonBarrelCells_position_y = nullptr;
  TTreeReaderArray<Float_t> *MuonBarrelCells_position_z = nullptr;

  // the hits in muon endcap
  TTreeReaderArray<ULong_t> *MuonEndcapHits_cellID = nullptr;
  TTreeReaderArray<Float_t> *MuonEndcapHits_energy = nullptr;
  // TTreeReaderArray<Double_t> *MuonEndcapHits_position_x = nullptr;
  // TTreeReaderArray<Double_t> *MuonEndcapHits_position_y = nullptr;
  // TTreeReaderArray<Double_t> *MuonEndcapHits_position_z = nullptr;
  TTreeReaderArray<Float_t> *MuonEndcapHits_position_x = nullptr;
  TTreeReaderArray<Float_t> *MuonEndcapHits_position_y = nullptr;
  TTreeReaderArray<Float_t> *MuonEndcapHits_position_z = nullptr;

  // the cells in muon endcap
  TTreeReaderArray<ULong_t> *MuonEndcapCells_cellID = nullptr;
  TTreeReaderArray<Float_t> *MuonEndcapCells_energy = nullptr;
  TTreeReaderArray<Float_t> *MuonEndcapCells_position_x = nullptr;
  TTreeReaderArray<Float_t> *MuonEndcapCells_position_y = nullptr;
  TTreeReaderArray<Float_t> *MuonEndcapCells_position_z = nullptr;

  // IF WE ARE SURE THAT WE DO NOT PLOT BOTH SW AND TOPOCLUSTERS AT THE SAME TIME
  // we could simplify the code of the EventDisplay class by having a single
  // cluster reader, that points to either SW or topo clusters based on the clusters
  // that we want to draw
  // For the moment keep both (we might want to have SW clusters in ECAL for
  // EM objects and topoclusters for ECAL+HCAL for jets)
  
  // the (corrected) calo topo clusters
  TTreeReaderArray<Float_t> *CaloTopoClusters_energy = nullptr;
  TTreeReaderArray<Float_t> *CaloTopoClusters_position_x = nullptr;
  TTreeReaderArray<Float_t> *CaloTopoClusters_position_y = nullptr;
  TTreeReaderArray<Float_t> *CaloTopoClusters_position_z = nullptr;
  TTreeReaderArray<Float_t> *CaloTopoClusters_theta = nullptr;
  TTreeReaderArray<Float_t> *CaloTopoClusters_phi = nullptr;
  TTreeReaderArray<UInt_t> *CaloTopoClusters_hits_begin = nullptr;
  TTreeReaderArray<UInt_t> *CaloTopoClusters_hits_end = nullptr;  
  
  // the (positioned) cells in the topo clusters
  TTreeReaderArray<ULong_t> *CaloTopoClusterCells_cellID = nullptr;
  TTreeReaderArray<Float_t> *CaloTopoClusterCells_energy = nullptr;
  TTreeReaderArray<Float_t> *CaloTopoClusterCells_position_x = nullptr;
  TTreeReaderArray<Float_t> *CaloTopoClusterCells_position_y = nullptr;
  TTreeReaderArray<Float_t> *CaloTopoClusterCells_position_z = nullptr;

  // the (corrected) sliding window clusters
  TTreeReaderArray<Float_t> *CaloClusters_energy = nullptr;
  TTreeReaderArray<Float_t> *CaloClusters_position_x = nullptr;
  TTreeReaderArray<Float_t> *CaloClusters_position_y = nullptr;
  TTreeReaderArray<Float_t> *CaloClusters_position_z = nullptr;
  TTreeReaderArray<Float_t> *CaloClusters_theta = nullptr;
  TTreeReaderArray<Float_t> *CaloClusters_phi = nullptr;
  TTreeReaderArray<UInt_t> *CaloClusters_hits_begin = nullptr;
  TTreeReaderArray<UInt_t> *CaloClusters_hits_end = nullptr;

  // the (positioned) cells in the sliding window clusters
  TTreeReaderArray<ULong_t> *CaloClusterCells_cellID = nullptr;
  TTreeReaderArray<Float_t> *CaloClusterCells_energy = nullptr;
  TTreeReaderArray<Float_t> *CaloClusterCells_position_x = nullptr;
  TTreeReaderArray<Float_t> *CaloClusterCells_position_y = nullptr;
  TTreeReaderArray<Float_t> *CaloClusterCells_position_z = nullptr;

  // the links from the clusters to their cells (podio ObjectIDs), used when
  // there is no dedicated collection with the cells of the clusters
  TTreeReaderArray<Int_t> *CaloTopoClusters_hits_index = nullptr;
  TTreeReaderArray<UInt_t> *CaloTopoClusters_hits_collectionID = nullptr;
  TTreeReaderArray<Int_t> *CaloClusters_hits_index = nullptr;
  TTreeReaderArray<UInt_t> *CaloClusters_hits_collectionID = nullptr;

  // the cells of the clusters in the current event, in the order of the
  // cluster hits: the cells of cluster i are [hits_begin[i], hits_end[i]).
  // Filled by loadEvent, from the cluster cell collection or from the links
  struct ClusterCells
  {
    std::vector<ULong_t> cellID;
    std::vector<Float_t> energy;
    std::vector<Float_t> x, y, z;
    size_t size() const { return energy.size(); }
  };
  ClusterCells topoClusterCells;
  ClusterCells swClusterCells;

  // the calorimeter cell collections in the file, by podio collection ID,
  // used to follow the links
  struct CellCollection
  {
    std::string name;
    std::unique_ptr<TTreeReaderArray<ULong_t>> cellID;
    std::unique_ptr<TTreeReaderArray<Float_t>> energy, x, y, z;
  };
  std::map<UInt_t, CellCollection> cellCollections;
  // collection IDs pointed to by links but not found in the file (warned once)
  std::set<UInt_t> missingCellCollections;

  
 public:
  EventReader(TFile* f, bool doHCal=false);
  ~EventReader();
  bool checkCollection(const std::string& drawFlag,
		       const std::string& flagName,
		       const std::string& requiredBranch,
		       std::string& collectionName);
  bool checkCollections(const std::string& drawFlag,
			const std::string& flagName1,
			const std::string& flagName2,
			const std::string& requiredBranch,
			std::string& collectionName1,
			std::string& collectionName2);
  template<class T>
  TTreeReaderArray<T>* makeArray(const std::string& collection, const std::string& branch)
  {
    if (collection == "")
      return new TTreeReaderArray<T>(*fReader, branch.c_str());
    else
      return new TTreeReaderArray<T>(*fReader, Form("%s.%s", collection.c_str(), branch.c_str()));
  }
  void setBranches();
  // set up the reading of the cells of the clusters, from the cell collection
  // (if present) or from the links of the clusters to the cell collections
  // returns false if neither is available
  bool setClusterCells(const std::string& clustersKey,
                       const std::string& cellsKey,
                       TTreeReaderArray<ULong_t>*& cellID,
                       TTreeReaderArray<Float_t>*& energy,
                       TTreeReaderArray<Float_t>*& x,
                       TTreeReaderArray<Float_t>*& y,
                       TTreeReaderArray<Float_t>*& z,
                       TTreeReaderArray<Int_t>*& linkIndex,
                       TTreeReaderArray<UInt_t>*& linkCollectionID);
  // find all the calorimeter cell collections in the file
  void setCellCollections();
  // fill the cells of the clusters for the current event
  void fillClusterCells(ClusterCells& cells,
                        TTreeReaderArray<ULong_t>* cellID,
                        TTreeReaderArray<Float_t>* energy,
                        TTreeReaderArray<Float_t>* x,
                        TTreeReaderArray<Float_t>* y,
                        TTreeReaderArray<Float_t>* z,
                        TTreeReaderArray<Int_t>* linkIndex,
                        TTreeReaderArray<UInt_t>* linkCollectionID);
  // the CellIDEncoding strings stored in the metadata of the file, by collection name
  std::map<std::string, std::string> readCellIDEncodings();
  void loadEvent(int event);
  bool m_doHCal;
  bool m_doSW;
  bool m_doTopo;
  bool m_readMergedCells;
};

#endif
