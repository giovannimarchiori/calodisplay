/******************************************************************************/
// Event display for the ALLEGRO detector for FCC-ee
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
  // a branch of the events tree, read for the current event
  // (null if the collection is not read)
  template <typename T>
  using Array = std::unique_ptr<TTreeReaderArray<T>>;

  // the reader of the events tree: declared first, so that it is destroyed
  // after the arrays, which deregister from it
  std::unique_ptr<TTreeReader> fReader;

  Int_t nEvents = 0; // Number of events in file

  // a collection of hits, digis or cells: their cellID (not read for the
  // simulated calorimeter hits and the merged ECal cells), energy and position
  // (double for the tracker hits and digis, float for the calorimeter ones)
  template <typename PositionT>
  struct HitCollection
  {
    Array<ULong_t> cellID;
    Array<Float_t> energy;
    Array<PositionT> x, y, z;
    bool isSet() const { return x != nullptr; }
    unsigned int size() const { return x ? x->GetSize() : 0; }
  };
  using TrackerHits = HitCollection<Double_t>;
  using CaloHits = HitCollection<Float_t>;

  // MC particles
  struct Particles
  {
    Array<Int_t> PDG, generatorStatus, simulatorStatus;
    Array<Float_t> charge, time;
    Array<Double_t> mass;
    Array<Double_t> vertex_x, vertex_y, vertex_z;
    Array<Double_t> endpoint_x, endpoint_y, endpoint_z;
    Array<Double_t> momentum_x, momentum_y, momentum_z;
    Array<UInt_t> parents_begin, parents_end, daughters_begin, daughters_end;
    Array<Int_t> parents_index, daughters_index;
    bool isSet() const { return PDG != nullptr; }
    unsigned int size() const { return PDG ? PDG->GetSize() : 0; }
  };

  // reconstructed tracks: the track states of track i are
  // [trackStates_begin[i], trackStates_end[i]), and its numbers of hits per
  // subdetector [subdetectorHitNumbers_begin[i], subdetectorHitNumbers_end[i])
  struct Tracks
  {
    Array<UInt_t> trackStates_begin, trackStates_end;
    Array<UInt_t> subdetectorHitNumbers_begin, subdetectorHitNumbers_end;
    Array<Int_t> trackState_location;
    Array<Float_t> trackState_omega, trackState_phi, trackState_tanLambda, trackState_d0, trackState_z0;
    Array<Float_t> trackState_referencePoint_x, trackState_referencePoint_y, trackState_referencePoint_z;
    Array<Int_t> subdetectorHitNumbers;
    bool isSet() const { return trackStates_begin != nullptr; }
    unsigned int size() const { return trackStates_begin ? trackStates_begin->GetSize() : 0; }
  };

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

  // calorimeter clusters (topo or sliding window) and their cells
  struct Clusters
  {
    std::string name; // name of the collection
    Array<Float_t> energy;
    Array<Float_t> x, y, z;
    Array<Float_t> theta, phi;
    Array<UInt_t> hits_begin, hits_end;
    // the dedicated collection with the cells of the clusters, if any
    CaloHits cellCollection;
    // otherwise, the links from the clusters to their cells (podio ObjectIDs)
    Array<Int_t> hits_index;
    Array<UInt_t> hits_collectionID;
    // the cells of the clusters in the current event
    ClusterCells cells;
    bool isSet() const { return energy != nullptr; }
    unsigned int size() const { return x ? x->GetSize() : 0; }
  };

  // particles and tracks
  Particles genParticles;
  Particles simParticles;
  Tracks tracks;

  // hits and digis in the tracker
  TrackerHits vertexBarrelHits, vertexEndcapHits;
  TrackerHits mainTrackerHits;
  TrackerHits siWrapperBarrelHits, siWrapperEndcapHits;
  TrackerHits vertexBarrelDigis, vertexEndcapDigis;
  TrackerHits mainTrackerDigis;
  TrackerHits siWrapperBarrelDigis, siWrapperEndcapDigis;

  // hits and cells in the calorimeters and the muon tagger
  CaloHits ecalBarrelHits, ecalEndcapHits;
  CaloHits hcalBarrelHits, hcalEndcapHits;
  CaloHits muonBarrelHits, muonEndcapHits;
  CaloHits ecalBarrelCells, ecalEndcapCells;
  // the cells in ECal barrel with coarser merging
  CaloHits ecalBarrelMergedCells;
  CaloHits hcalBarrelCells, hcalEndcapCells;
  CaloHits muonBarrelCells, muonEndcapCells;

  // the clusters: one or more collections of topo clusters (e.g. for
  // ECAL+HCAL jets) and of sliding-window clusters (e.g. in ECAL for EM
  // objects), e.g. in the barrel and in the endcap, can be drawn at the same time
  // (the vectors are filled by setBranches and not changed afterwards)
  std::vector<Clusters> topoClusters;
  std::vector<Clusters> swClusters;

  // the calorimeter cell collections in the file, by podio collection ID,
  // used to follow the links
  struct CellCollection
  {
    std::string name;
    CaloHits cells;
  };
  std::map<UInt_t, CellCollection> cellCollections;
  // collection IDs pointed to by links but not found in the file (warned once)
  std::set<UInt_t> missingCellCollections;

 public:
  EventReader(TFile* f);
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
  Array<T> makeArray(const std::string& collection, const std::string& branch)
  {
    if (collection == "")
      return std::make_unique<TTreeReaderArray<T>>(*fReader, branch.c_str());
    else
      return std::make_unique<TTreeReaderArray<T>>(*fReader, Form("%s.%s", collection.c_str(), branch.c_str()));
  }
  // set up the reading of a collection of hits, digis or cells, with the
  // energy in branch energyBranch and the position in positionBranch.x/y/z
  template<typename PositionT>
  void setHits(HitCollection<PositionT>& hits, const std::string& collection,
               const std::string& energyBranch, const std::string& positionBranch,
               bool readCellID)
  {
    if (readCellID)
      hits.cellID = makeArray<ULong_t>(collection, "cellID");
    hits.energy = makeArray<Float_t>(collection, energyBranch);
    hits.x = makeArray<PositionT>(collection, positionBranch + ".x");
    hits.y = makeArray<PositionT>(collection, positionBranch + ".y");
    hits.z = makeArray<PositionT>(collection, positionBranch + ".z");
  }
  void setParticles(Particles& particles, const std::string& collection);
  void setClusters(Clusters& clusters, const std::string& collection);
  // set up the reading of the collections of clusters listed (comma-separated)
  // in the configuration parameter clustersKey, and of their cells (the
  // collections of cells listed in cellsKey, in the same order), if drawFlag is set
  void setClusterCollections(std::vector<Clusters>& clustersList, const std::string& drawFlag,
                             const std::string& clustersKey, const std::string& cellsKey);
  void setBranches();
  // set up the reading of the cells of the clusters, from the cell collection
  // (if present) or from the links of the clusters to the cell collections
  // returns the cell collection used ("" if the links are used, or neither
  // is available)
  std::string setClusterCells(Clusters& clusters, std::string collection);
  // find all the calorimeter cell collections in the file
  void setCellCollections();
  // fill the cells of the clusters for the current event
  void fillClusterCells(Clusters& clusters);
  // the CellIDEncoding strings stored in the metadata of the file, by collection name
  std::map<std::string, std::string> readCellIDEncodings();
  // the options of the Gaudi job that produced the file, stored in the
  // configuration metadata, as "<algorithm or tool>.<property>" -> value
  std::map<std::string, std::string> readGaudiOptions();
  void loadEvent(int event);
};

#endif
