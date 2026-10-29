/******************************************************************************/
// Event display for the ALLEGRO detector for FCC-ee
// author: Giovanni Marchiori (giovanni.marchiori@cern.ch)
//
// class EventReader: sets up the TTreeReader branches for reading event data
//
/******************************************************************************/

#include "EventReader.h"
#include "Globals.h"
#include "PodioUtils.h"

#include <TBranchElement.h>

#include <sstream>

namespace
{
// the names in a comma-separated list (ignoring the spaces around them)
std::vector<std::string> splitList(const std::string& list)
{
  std::vector<std::string> names;
  std::stringstream ss(list);
  std::string name;
  while (std::getline(ss, name, ','))
  {
    const size_t first = name.find_first_not_of(" \t");
    if (first == std::string::npos)
      continue;
    names.push_back(name.substr(first, name.find_last_not_of(" \t") - first + 1));
  }
  return names;
}

// the comma-separated list of names (empty names are kept, to keep the positions)
std::string joinList(const std::vector<std::string>& names)
{
  std::string list;
  for (size_t i = 0; i < names.size(); i++)
    list += (i == 0 ? "" : ",") + names[i];
  return list;
}
}

EventReader::EventReader(TFile* f)
{
  fReader = std::make_unique<TTreeReader>("events", f);
  nEvents = fReader->GetEntries();
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

void EventReader::setParticles(Particles& particles, const std::string& collection)
{
  particles.PDG = makeArray<Int_t>(collection, "PDG");
  particles.generatorStatus = makeArray<Int_t>(collection, "generatorStatus");
  particles.simulatorStatus = makeArray<Int_t>(collection, "simulatorStatus");
  particles.charge = makeArray<Float_t>(collection, "charge");
  particles.time = makeArray<Float_t>(collection, "time");
  particles.mass = makeArray<Double_t>(collection, "mass");
  particles.vertex_x = makeArray<Double_t>(collection, "vertex.x");
  particles.vertex_y = makeArray<Double_t>(collection, "vertex.y");
  particles.vertex_z = makeArray<Double_t>(collection, "vertex.z");
  particles.endpoint_x = makeArray<Double_t>(collection, "endpoint.x");
  particles.endpoint_y = makeArray<Double_t>(collection, "endpoint.y");
  particles.endpoint_z = makeArray<Double_t>(collection, "endpoint.z");
  particles.momentum_x = makeArray<Double_t>(collection, "momentum.x");
  particles.momentum_y = makeArray<Double_t>(collection, "momentum.y");
  particles.momentum_z = makeArray<Double_t>(collection, "momentum.z");
  particles.parents_begin = makeArray<UInt_t>(collection, "parents_begin");
  particles.parents_end = makeArray<UInt_t>(collection, "parents_end");
  particles.parents_index = makeArray<Int_t>("_" + collection + "_parents", "index");
  particles.daughters_begin = makeArray<UInt_t>(collection, "daughters_begin");
  particles.daughters_end = makeArray<UInt_t>(collection, "daughters_end");
  particles.daughters_index = makeArray<Int_t>("_" + collection + "_daughters", "index");
}

void EventReader::setClusters(Clusters& clusters, const std::string& collection)
{
  clusters.energy     = makeArray<Float_t>(collection, "energy");
  clusters.x          = makeArray<Float_t>(collection, "position.x");
  clusters.y          = makeArray<Float_t>(collection, "position.y");
  clusters.z          = makeArray<Float_t>(collection, "position.z");
  clusters.theta      = makeArray<Float_t>(collection, "iTheta");
  clusters.phi        = makeArray<Float_t>(collection, "iPhi");
  clusters.hits_begin = makeArray<UInt_t> (collection, "hits_begin");
  clusters.hits_end   = makeArray<UInt_t> (collection, "hits_end");
}

void EventReader::setBranches()
{
  std::string collection, collection2;

  // primary and secondary particles
  if (checkCollection("drawGenParticles", "genParticles", "PDG", collection))
    setParticles(genParticles, collection);
  if (checkCollection("drawSimParticles", "simParticles", "PDG", collection))
    setParticles(simParticles, collection);

  // reconstructed tracks
  if (checkCollection("drawTracks", "tracks", "trackStates_begin", collection))
  {
    std::string trackStateColl = "_" + collection + "_trackStates";
    tracks.trackStates_begin = makeArray<UInt_t>(collection, "trackStates_begin");
    tracks.trackStates_end = makeArray<UInt_t>(collection, "trackStates_end");
    tracks.subdetectorHitNumbers_begin = makeArray<UInt_t>(collection, "subdetectorHitNumbers_begin");
    tracks.subdetectorHitNumbers_end = makeArray<UInt_t>(collection, "subdetectorHitNumbers_end");
    tracks.trackState_location = makeArray<Int_t>(trackStateColl, "location");
    tracks.trackState_omega = makeArray<Float_t>(trackStateColl, "omega");
    tracks.trackState_phi = makeArray<Float_t>(trackStateColl, "phi");
    tracks.trackState_tanLambda = makeArray<Float_t>(trackStateColl, "tanLambda");
    tracks.trackState_d0 = makeArray<Float_t>(trackStateColl, "D0");
    tracks.trackState_z0 = makeArray<Float_t>(trackStateColl, "Z0");
    tracks.trackState_referencePoint_x = makeArray<Float_t>(trackStateColl, "referencePoint.x");
    tracks.trackState_referencePoint_y = makeArray<Float_t>(trackStateColl, "referencePoint.y");
    tracks.trackState_referencePoint_z = makeArray<Float_t>(trackStateColl, "referencePoint.z");
    tracks.subdetectorHitNumbers = makeArray<Int_t>("", "_" + collection + "_subdetectorHitNumbers");
  }

  // hits and digis in the tracker (SimTrackerHits and TrackerHits)
  if (checkCollections("drawVertexHits", "vertexBarrelHits", "vertexEndcapHits", "cellID", collection, collection2))
  {
    setHits(vertexBarrelHits, collection, "eDep", "position", true);
    setHits(vertexEndcapHits, collection2, "eDep", "position", true);
  }
  if (checkCollection("drawMainTrackerHits", "mainTrackerHits", "cellID", collection))
    setHits(mainTrackerHits, collection, "eDep", "position", true);
  if (checkCollections("drawSiWrapperHits", "siWrapperBarrelHits", "siWrapperEndcapHits", "cellID", collection, collection2))
  {
    setHits(siWrapperBarrelHits, collection, "eDep", "position", true);
    setHits(siWrapperEndcapHits, collection2, "eDep", "position", true);
  }
  if (checkCollections("drawVertexDigis", "vertexBarrelDigis", "vertexEndcapDigis", "cellID", collection, collection2))
  {
    setHits(vertexBarrelDigis, collection, "eDep", "position", true);
    setHits(vertexEndcapDigis, collection2, "eDep", "position", true);
  }
  if (checkCollection("drawMainTrackerDigis", "mainTrackerDigis", "cellID", collection))
    setHits(mainTrackerDigis, collection, "eDep", "position", true);
  if (checkCollections("drawSiWrapperDigis", "siWrapperBarrelDigis", "siWrapperEndcapDigis", "cellID", collection, collection2))
  {
    setHits(siWrapperBarrelDigis, collection, "eDep", "position", true);
    setHits(siWrapperEndcapDigis, collection2, "eDep", "position", true);
  }

  // hits (CaloHitContributions, at the position of the step) and cells
  // (CalorimeterHits) in the ECal
  if (checkCollection("drawECalBarrelHits", "ecalBarrelHits", "energy", collection))
    setHits(ecalBarrelHits, collection, "energy", "stepPosition", false);
  if (checkCollection("drawECalBarrelCells", "ecalBarrelCells", "energy", collection))
    setHits(ecalBarrelCells, collection, "energy", "position", true);
  if (checkCollection("drawECalBarrelMergedCells", "ecalBarrelMergedCells", "energy", collection))
    setHits(ecalBarrelMergedCells, collection, "energy", "position", false);
  if (checkCollection("drawECalEndcapHits", "ecalEndcapHits", "energy", collection))
    setHits(ecalEndcapHits, collection, "energy", "stepPosition", false);
  if (checkCollection("drawECalEndcapCells", "ecalEndcapCells", "energy", collection))
    setHits(ecalEndcapCells, collection, "energy", "position", true);

  // hits and cells in the HCal
  if (checkCollection("drawHCalBarrelHits", "hcalBarrelHits", "energy", collection))
    setHits(hcalBarrelHits, collection, "energy", "stepPosition", false);
  if (checkCollection("drawHCalBarrelCells", "hcalBarrelCells", "energy", collection))
    setHits(hcalBarrelCells, collection, "energy", "position", true);
  if (checkCollection("drawHCalEndcapHits", "hcalEndcapHits", "energy", collection))
    setHits(hcalEndcapHits, collection, "energy", "stepPosition", false);
  if (checkCollection("drawHCalEndcapCells", "hcalEndcapCells", "energy", collection))
    setHits(hcalEndcapCells, collection, "energy", "position", true);

  // hits and cells in the muon tagger
  if (checkCollections("drawMuonHits", "muonBarrelHits", "muonEndcapHits", "energy", collection, collection2))
  {
    setHits(muonBarrelHits, collection, "energy", "stepPosition", false);
    setHits(muonEndcapHits, collection2, "energy", "stepPosition", false);
  }
  if (checkCollections("drawMuonCells", "muonBarrelCells", "muonEndcapCells", "energy", collection, collection2))
  {
    setHits(muonBarrelCells, collection, "energy", "position", true);
    setHits(muonEndcapCells, collection2, "energy", "position", true);
  }

  // SW and topo clusters and their cells
  setClusterCollections(swClusters, "drawCaloClusters", "caloClusters", "caloClusterCells");
  setClusterCollections(topoClusters, "drawTopoClusters", "topoClusters", "topoClusterCells");
}

void EventReader::setClusterCollections(std::vector<Clusters>& clustersList, const std::string& drawFlag,
                                        const std::string& clustersKey, const std::string& cellsKey)
{
  if (!displayConfig.getBoolConfig(drawFlag))
    return;
  const std::vector<std::string> names = splitList(displayConfig.getStringConfig(clustersKey));
  const std::vector<std::string> cellNames = splitList(displayConfig.getStringConfig(cellsKey));
  if (cellNames.size() > names.size())
    std::cout << "WARNING: more collections in " << cellsKey << " than in " << clustersKey
              << ", ignoring the last ones" << std::endl;

  // (reserve, so that the clusters are not moved while their branches are set)
  clustersList.reserve(names.size());
  std::vector<std::string> usedNames, usedCellNames;
  for (size_t i = 0; i < names.size(); i++)
  {
    const std::string& name = names[i];
    if (!fReader->GetTree()->FindBranch(Form("%s.energy", name.c_str())))
    {
      std::cout << "WARNING: branch " << name << ".energy not found, not drawing the clusters " << name << std::endl;
      continue;
    }
    Clusters& clusters = clustersList.emplace_back();
    clusters.name = name;
    setClusters(clusters, name);
    usedNames.push_back(name);
    usedCellNames.push_back(setClusterCells(clusters, i < cellNames.size() ? cellNames[i] : ""));
  }
  if (usedNames.empty())
  {
    std::cout << "WARNING: no collection of " << clustersKey << " found, disabling " << clustersKey << std::endl;
    displayConfig.setBoolConfig(drawFlag, false);
  }
  displayConfig.setStringConfig(clustersKey, joinList(usedNames));
  displayConfig.setStringConfig(cellsKey, joinList(usedCellNames));
}

std::string EventReader::setClusterCells(Clusters& clusters, std::string collection)
{
  const std::string& clustersName = clusters.name;
  if (collection == "") {
    // try to use default collection
    collection = clustersName;
    collection.pop_back();
    collection += "Cells";
    std::cout << "WARNING: cell collection of " << clustersName << " not set, trying with " << collection << std::endl;
  }

  // 1. a dedicated collection with the cells of the clusters
  if (fReader->GetTree()->FindBranch(Form("%s.cellID", collection.c_str()))) {
    setHits(clusters.cellCollection, collection, "energy", "position", true);
    return collection;
  }

  // 2. the links of the clusters to the cells in the cell collections
  std::string links = "_" + clustersName + "_hits";
  if (fReader->GetTree()->FindBranch(Form("%s.collectionID", links.c_str()))) {
    std::cout << "Branch " << collection << ".cellID not found, using the links " << links
              << " of the clusters to the cell collections" << std::endl;
    clusters.hits_index        = makeArray<Int_t>(links, "index");
    clusters.hits_collectionID = makeArray<UInt_t>(links, "collectionID");
    if (cellCollections.empty())
      setCellCollections();
    return "";
  }

  std::cout << "WARNING: neither branch " << collection << ".cellID nor " << links
            << " found, drawing the clusters without their cells" << std::endl;
  return "";
}

void EventReader::setCellCollections()
{
  // the podio collection ID is a hash of the collection name, so the
  // collections can be identified without reading the podio metadata
  for (TObject* obj : *fReader->GetTree()->GetListOfBranches()) {
    TBranchElement* branch = dynamic_cast<TBranchElement*>(obj);
    if (!branch || std::string(branch->GetClassName()) != "vector<edm4hep::CalorimeterHitData>")
      continue;
    std::string name = branch->GetName();
    CellCollection& collection = cellCollections[PodioUtils::collectionID(name)];
    collection.name = name;
    setHits(collection.cells, name, "energy", "position", true);
  }
  std::cout << "Calorimeter cell collections available for the cluster links:";
  for (const auto& [id, collection] : cellCollections)
    std::cout << " " << collection.name;
  std::cout << std::endl;
}

void EventReader::fillClusterCells(Clusters& clusters)
{
  ClusterCells& cells = clusters.cells;
  cells.cellID.clear();
  cells.energy.clear();
  cells.x.clear();
  cells.y.clear();
  cells.z.clear();

  auto addCell = [&cells](ULong_t id, Float_t e, Float_t xx, Float_t yy, Float_t zz) {
    cells.cellID.push_back(id);
    cells.energy.push_back(e);
    cells.x.push_back(xx);
    cells.y.push_back(yy);
    cells.z.push_back(zz);
  };

  const CaloHits& cellCollection = clusters.cellCollection;
  if (cellCollection.isSet()) {
    // dedicated collection with the cells of the clusters
    for (unsigned int i = 0; i < cellCollection.size(); i++)
      addCell((*cellCollection.cellID)[i], (*cellCollection.energy)[i],
              (*cellCollection.x)[i], (*cellCollection.y)[i], (*cellCollection.z)[i]);
  }
  else if (clusters.hits_index) {
    // follow the links to the cell collections
    for (unsigned int i = 0; i < clusters.hits_index->GetSize(); i++) {
      UInt_t collectionID = (*clusters.hits_collectionID)[i];
      Int_t index = (*clusters.hits_index)[i];
      auto it = cellCollections.find(collectionID);
      if (it == cellCollections.end() || index < 0 || (unsigned int) index >= it->second.cells.size()) {
        if (it == cellCollections.end() && missingCellCollections.insert(collectionID).second)
          std::cout << "WARNING: cell collection with podio ID " << collectionID
                    << " used by the cluster links is not in the file, its cells are not drawn" << std::endl;
        // keep the cell (with cellID 0) so that hits_begin/hits_end stay aligned
        addCell(0, 0., 0., 0., 0.);
        continue;
      }
      const CaloHits& c = it->second.cells;
      addCell((*c.cellID)[index], (*c.energy)[index], (*c.x)[index], (*c.y)[index], (*c.z)[index]);
    }
  }
}

std::map<std::string, std::string> EventReader::readCellIDEncodings()
{
  // podio stores the file metadata as parameters (keys and values) in the
  // "metadata" tree, with keys <collection>__CellIDEncoding
  std::map<std::string, std::string> encodings;
  TFile* file = fReader->GetTree()->GetCurrentFile();
  TTree* tree = file ? dynamic_cast<TTree*>(file->Get("metadata")) : nullptr;
  if (!tree || !tree->GetBranch("GPStringKeys") || !tree->GetBranch("GPStringValues"))
    return encodings;

  TTreeReader reader(tree);
  TTreeReaderValue<std::vector<std::string>> keys(reader, "GPStringKeys");
  TTreeReaderValue<std::vector<std::vector<std::string>>> values(reader, "GPStringValues");
  if (!reader.Next())
    return encodings;

  const std::string suffix = "__CellIDEncoding";
  for (size_t i = 0; i < keys->size() && i < values->size(); i++) {
    const std::string& key = (*keys)[i];
    if (key.size() > suffix.size() && key.compare(key.size() - suffix.size(), suffix.size(), suffix) == 0
        && !(*values)[i].empty())
      encodings[key.substr(0, key.size() - suffix.size())] = (*values)[i][0];
  }
  return encodings;
}

std::map<std::string, std::string> EventReader::readGaudiOptions()
{
  // k4FWCore stores the options in the "configuration_metadata" tree, under the
  // key gaudiConfigOptions, as strings <name>.<property> = <value>, with the
  // value between double quotes
  std::map<std::string, std::string> options;
  TFile* file = fReader->GetTree()->GetCurrentFile();
  TTree* tree = file ? dynamic_cast<TTree*>(file->Get("configuration_metadata")) : nullptr;
  if (!tree || !tree->GetBranch("GPStringKeys") || !tree->GetBranch("GPStringValues"))
    return options;

  TTreeReader reader(tree);
  TTreeReaderValue<std::vector<std::string>> keys(reader, "GPStringKeys");
  TTreeReaderValue<std::vector<std::vector<std::string>>> values(reader, "GPStringValues");
  if (!reader.Next())
    return options;

  for (size_t i = 0; i < keys->size() && i < values->size(); i++) {
    if ((*keys)[i] != "gaudiConfigOptions")
      continue;
    for (std::string option : (*values)[i]) {
      // each option ends with a newline
      option.erase(option.find_last_not_of(" \n\r") + 1);
      const size_t eq = option.find(" = ");
      if (eq == std::string::npos)
        continue;
      std::string name = option.substr(0, eq);
      std::string value = option.substr(eq + 3);
      name.erase(0, name.find_first_not_of(' '));
      if (value.size() >= 2 && value.front() == '"' && value.back() == '"')
        value = value.substr(1, value.size() - 2);
      options[name] = value;
    }
  }
  return options;
}

void EventReader::loadEvent(int event) {
  if (event<0 || event>=nEvents) {
    std::cout << "Event number " << event << " out of range : [0, " << nEvents-1 << "]" << std::endl;
  }
  fReader->SetEntry(event);

  for (std::vector<Clusters>* clustersList : {&topoClusters, &swClusters})
    for (Clusters& clusters : *clustersList)
      fillClusterCells(clusters);
}
