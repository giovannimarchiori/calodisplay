/******************************************************************************/
// Event display for the ALLEGRO detector for FCC-ee
// author: Giovanni Marchiori (giovanni.marchiori@cern.ch)
//
// class EventDisplay: creates and populates the event display
//
/******************************************************************************/

#ifndef EVENTDISPLAY_H
#define EVENTDISPLAY_H

/******************************************************************************/
// dependencies
/******************************************************************************/
#include <TEveViewer.h>
#include <TEveBrowser.h>
#include <TEveGeoShape.h>
#include <TEveScene.h>
#include <TEveProjectionManager.h>
#include <TEveElement.h>
#include <TEveTrackPropagator.h>
#include <TEveTrack.h>
#include <TEvePointSet.h>
#include <TEveQuadSet.h>
#include <TEveRGBAPalette.h>
#include <TEveBoxSet.h>
#include <TEveStraightLineSet.h>

#include <TGLViewer.h>

#include <TGNumberEntry.h>
#include <TGLabel.h>

#include <TVector3.h>
#include <TPRegexp.h>

#include <memory>

#include "DetectorGeometry.h"
#include "GeometryCache.h"
#include "EventReader.h"
#include "CaloCluster.h"
#include "GraphicUtils.h"
#include "Globals.h"


/******************************************************************************/
// auxiliary classes/functions
/******************************************************************************/

// struct to hold start/end/momentum info to calculate e- brem from MCParticles
struct ElectronSegment
{
  TVector3 start;
  TVector3 end;
  TVector3 momentum;
};

// derived TEvePointSet class that shows point info in console when selected
class MyPointSet : public TEvePointSet {
public:
  std::vector<std::string> fTooltips;  // per-point info
  bool fShowTooltips;                  // flag to enable/disable tooltips

  MyPointSet(const char* name = "MyPoints", bool showTooltips = true)
    : TEvePointSet(name), fShowTooltips(showTooltips) {
    if (showTooltips) {
        SetPickable(kTRUE);
        SetOwnIds(kTRUE);
      }
  }

  void Reset(Int_t n_points = 0, Int_t n_int_ids = 0) {
    // 1. call base class Reset
    TEvePointSet::Reset(n_points, n_int_ids);

    // 2. restore your desired behaviour
    if (fShowTooltips) {
      SetOwnIds(kTRUE);
      SetPickable(kTRUE);
    }

    // 3. resets the tooltips vector
    fTooltips.clear();
  }

  void PointSelected(Int_t idx) override {
    if (!fShowTooltips) return;
    std::cout << "tooltips size: " << fTooltips.size() << std::endl;
    if (idx < 0 || idx >= (Int_t)fTooltips.size()) return;
    std::cout << "Point selected: idx = " << idx << ", " << fTooltips[idx] << std::endl;
  }
};


/******************************************************************************/
// main class
/******************************************************************************/

class EventDisplay
{

public:
  /******************************************************************************/
  // SETTINGS FOR THE DISPLAY
  // - geometry file (and corresponding merging of theta cells/modules vs layer
  // - data file
  // - flags to draw or not elements of the event (particles, hits, ..)
  // - minimum energy of particles, hits, cells, clusters
  /******************************************************************************/

  bool debug = false;
  bool doHCal = false;
  bool doEndcaps = false;
  bool showFullDetector = false;
  bool drawClustersBarycenterVsLayer = true;

  bool useTransparencies = false;
  
  // will read this from config file
  // tried to make them private but code was behaving weirdly..
  // min particle energy (GeV)
  float ParticleEnergyThreshold = 0.0;
  // min hit and cell energy (GeV)
  float HitEnergyThreshold = 0.0;
  float CellEnergyThreshold = 0.0;
  // min cluster energy (in GeV)
  float ClusterEnergyThreshold = 0.0;
  // max radius for propagator
  float rMax = 0.0;

  // G4 geometry file
  bool useG4geom = true;
  std::string geomFile = "data/ALLEGRO_o1_v03_noDCHcells.root";
  int detectorVersion = 3;
  // read the G4 geometry without the volumes that are never drawn from a
  // cache (see GeometryCache.h), writing it the first time
  bool useGeometryCache = true;

  // merging along theta and module directions and drawing options
  // std::string evtFile = "output_fullCalo_SimAndDigi_withTopoCluster_MagneticField_False_pMin_10000_MeV_ThetaMinMax_40_140_pdgId_11_pythiaFalse_NoiseFalse.root";
  // std::string evtFile = "test.root";
  // const std::vector<int> mergedCells_Theta = {4, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};
  // const std::vector<int> mergedModules = {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
  // bool drawMergedCells = false;

  // std::string evtFile = "output_fullCalo_SimAndDigi_withTopoCluster_MagneticField_False_pMin_10000_MeV_ThetaMinMax_40_140_pdgId_11_pythiaFalse_NoiseFalse_testmerge.root";
  // const std::vector<int> mergedCells_Theta = {4, 8, 2, 1, 8, 4, 2, 1, 4, 2, 1, 8};
  // const std::vector<int> mergedModules = {2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1};
  // bool drawMergedCells = true;

  std::string evtFile = "";

  Int_t nEvents = 0; // Number of events in file
  Int_t eventId = 0; // Current event id

  /******************************************************************************/
  // ELEMENTS OF THE EVENT DISPLAY
  /******************************************************************************/
  // Magnetic field: owned by the display and shared by the track propagators
  // (which must not own it, or it would be deleted by each of them)
  std::unique_ptr<TEveMagField> magField;

  // MC particles
  TEveTrackList *particles = nullptr;
  
  // Sim hits
  MyPointSet *vtxHits = nullptr;
  MyPointSet *dchHits = nullptr;
  MyPointSet *siwrHits = nullptr;
  MyPointSet *ecalHits = nullptr;
  MyPointSet *hcalHits = nullptr;
  MyPointSet *muonHits = nullptr;
  TEveElementList *hits = nullptr;

  // Digitised hits
  MyPointSet *vtxDigis = nullptr;
  MyPointSet *dchDigis = nullptr;
  MyPointSet *siwrDigis = nullptr;
  MyPointSet *ecalCells = nullptr;
  MyPointSet *hcalCells = nullptr;
  MyPointSet *muonCells = nullptr;
  MyPointSet *cells_merged = nullptr;
  TEveElementList *digis = nullptr;

  // Reconstructed tracks
  TEveTrackList *tracks = nullptr;

  // clusters: the graphic elements of the clusters of a type (topo or
  // sliding window); the cells are drawn only for the ECAL and HCAL barrel
  struct ClusterSet
  {
    std::string type;    // used in the names of the elements
    std::string label;   // used in the printouts
    Color_t centerColor; // colour of the barycenters
    bool drawWindow;     // draw the window of the sliding-window clusters
    // the clusters in the event reader (set when the event file is opened)
    const EventReader::Clusters *clusters = nullptr;
    // the clusters above threshold, with their energy and barycenter per layer
    std::vector<std::unique_ptr<CaloCluster>> data;
    // the lists of the clusters in 3D and in the projections
    TEveElementList *list3D = nullptr;
    TEveElementList *listRhoZ = nullptr;
    TEveElementList *listRhoPhi = nullptr;
    // the barycenters and directions of the clusters in 3D
    TEvePointSet *center = nullptr;
    TEveStraightLineSet *direction = nullptr;
  };
  ClusterSet topoClusterSet{"topo", "topo", kWhite, false};
  ClusterSet swClusterSet{"sw", "SW", kGray, true};

  //
  TEveGeoShape *ecalbarrel = nullptr;
  TEveElementList *ecalendcap = nullptr;
  TEveGeoShape *hcalbarrel = nullptr;
  TEveElementList *hcalendcap = nullptr;

  // the various views/scenes/projection managers
  TEveViewer *rhoPhiView = nullptr;
  TEveViewer *rhoZView = nullptr;
  TGLViewer *rhoPhiGLView = nullptr;
  TGLViewer *rhoZGLView = nullptr;
  TGLViewer *mainGLView = nullptr;
  TEveScene *rhoPhiScene = nullptr;
  TEveScene *rhoPhiEventScene = nullptr;
  TEveScene *rhoPhiEventSceneManual = nullptr;
  TEveScene *rhoZScene = nullptr;
  TEveScene *rhoZEventScene = nullptr;
  TEveScene *rhoZEventSceneManual = nullptr;
  TEveProjectionManager *rhoPhiProjManager = nullptr;
  TEveProjectionManager *rhoZProjManager = nullptr;

  // the list of the geometry
  TEveElementList* geom = nullptr;
  // the full G4 geometry, kept until the display is set up to write the
  // geometry cache, and the elements created from it
  TEveGeoShapeExtract* geometryExtract = nullptr;
  GeometryCache::ExtractMap geometryExtracts;
  // identifies the display options that change which volumes are drawn
  std::string geometryCacheTag() const;

  // the list of the readout elements
  TEveElementList *readout = nullptr;
  TGLTextLabel *processLabel = nullptr;
  // palettes for the energy of the cells in the clusters (values in MeV) in
  // each view, shared by all events, with a legend in the projections; their
  // range is set from the largest cell energy of the event in the view
  TEveRGBAPalette *cellPalette3D = nullptr;
  TEveRGBAPalette *cellPaletteRhoPhi = nullptr;
  TEveRGBAPalette *cellPaletteRhoZ = nullptr;
  int maxCellValue3D = 0;
  int maxCellValueRhoPhi = 0;
  int maxCellValueRhoZ = 0;
  // the legends of the palettes and their titles, that can be hidden
  std::vector<TGLOverlayElement *> paletteLegends;
  bool showPaletteLegends = true;
  TGLTextLabel *eventLabel = nullptr;
  // the event number entries and the labels with the range of the event
  // numbers in the navigation bar (one per tab)
  std::vector<TGNumberEntry *> eventEntries;
  std::vector<TGLabel *> eventRangeLabels;

  std::unique_ptr<DetectorGeometry> geomReader;
  // the file must be declared before the reader, so that the reader is destroyed first
  std::unique_ptr<TFile> eventFile;
  std::unique_ptr<EventReader> eventReader;

  TGLViewer* activeGLViewer = nullptr;

  TEveBrowser* browser = nullptr;

  bool initRhoPhiView = false;
  bool initRhoZView = false;
  // the geometry and the readout have been projected in the views
  bool rhoPhiProjectionBuilt = false;
  bool rhoZProjectionBuilt = false;

  bool showReadout = true;

  /******************************************************************************/
  // class methods
  /******************************************************************************/

  // constructor - not needed, use compiler-generated one
  // EventDisplay();

  // start event display
  void startDisplay(int initialEvent = 0);

  // fill cluster objects for given collection (sw or topo) and loaded event
  // invoked by drawEvent
  void fillClusters(ClusterSet &set);

  // draw the MC particles
  void drawParticles();

  // draw the hits
  void drawHits();

  // draw the digis
  void drawDigis();

  // draw the reconstructed tracks
  void drawTracks();
  
  // draw the clusters created by fillClusters
  void drawClusters(ClusterSet &set);

  // draw the event
  void drawEvent(int event);

  // move to next event
  void fwd();

  // move to previous event
  void bck();

  // move to the event typed in the event number entry
  void goToEvent();

  // move to the next/previous event with the up/down arrows of the event
  // number entry (val: the argument of TGNumberEntry::ValueChanged)
  void stepEvent(Long_t val);

  // show the current event and the range of the event numbers in the navigation bar
  void updateEventEntries();

  // show/hide the readouts
  void toggleReadout();

  // show/hide the legends of the cluster cell energies
  void togglePaletteLegends();

  // create the buttons for navigating through the events
  void makeGui();

  // keep track of active GL viewer
  // void onRhoPhiViewActivated();
  // void onRhoZViewActivated();
  // void on3DViewActivated();
  void onTabSelected(Int_t tab);

  // take a screenshot
  void takeScreenshot();

  // find path marks for electron brem to make more realistic rendering of electron trajectories
  std::vector<TEvePathMarkD> getElectronBremsstrahlungMarks(
    unsigned int electronIndex,
    const TTreeReaderArray<Int_t>& daughters,
    const TTreeReaderArray<Int_t>& pdg,
    const TTreeReaderArray<Double_t>& vx,
    const TTreeReaderArray<Double_t>& vy,
    const TTreeReaderArray<Double_t>& vz,
    const TTreeReaderArray<Double_t>& px,
    const TTreeReaderArray<Double_t>& py,
    const TTreeReaderArray<Double_t>& pz,
    const TTreeReaderArray<UInt_t>& daughters_begin,
    const TTreeReaderArray<UInt_t>& daughters_end);

  std::vector<ElectronSegment> getElectronBremsstrahlungSegments(
    unsigned int electronIndex,
    const TTreeReaderArray<Int_t>& daughters,
    const TTreeReaderArray<Int_t>& pdg,
    const TTreeReaderArray<Double_t>& vx,
    const TTreeReaderArray<Double_t>& vy,
    const TTreeReaderArray<Double_t>& vz,
    const TTreeReaderArray<Double_t>& px,
    const TTreeReaderArray<Double_t>& py,
    const TTreeReaderArray<Double_t>& pz,
    const TTreeReaderArray<UInt_t>& daughters_begin,
    const TTreeReaderArray<UInt_t>& daughters_end);

  void createElementList(TEveElementList*& list, const std::string& listName);

  // draw a collection of hits, digis or cells with energy above Ethr
  // (negative for no threshold in the name of the collection)
  template <typename T> void drawHitCollection(MyPointSet*& hitSet, const std::string& name,
                                               const EventReader::HitCollection<T>& hitCollection,
                                               const float Ethr,
                                               bool digi = false,    // if true, use style for digis instead of hits
                                               bool addTooltips = false,  // if true, add tooltips (do it only for debug!)
                                               bool append = false)  // if append is true, and the pointer is not null, the collection is not reset
  {
    if (hitSet == nullptr)
    {
      hitSet = new MyPointSet("MyPointSet", addTooltips && digi && debug);
      if (Ethr>=0.0)
        hitSet->SetName(Form("%s (E>%.2f GeV)", name.c_str(), Ethr));
      else
        hitSet->SetName(name.c_str());

      if (digi) {
        hitSet->SetMarkerStyle(digiMarkerStyle);
        hitSet->SetMarkerSize(digiMarkerSize);
        hitSet->SetMarkerColor(digiMarkerColor);
        digis->AddElement(hitSet);
      }
      else {
        hitSet->SetMarkerStyle(hitsMarkerStyle);
        hitSet->SetMarkerSize(hitsMarkerSize);
        hitSet->SetMarkerColor(hitsMarkerColor);
        hits->AddElement(hitSet);
      }
    }
    else {
      if (!append)
        hitSet->Reset();
    }

    for (unsigned int i = 0; i < hitCollection.size(); i++)
    {
      // skip hit if E<threshold
      if (hitCollection.energy && (*hitCollection.energy)[i] < Ethr)
        continue;
      hitSet->SetNextPoint((*hitCollection.x)[i] * Units::mm,
                           (*hitCollection.y)[i] * Units::mm,
                           (*hitCollection.z)[i] * Units::mm);
    }
  }

private:
  struct SubDetector
  {
    std::string name;
    std::string regex;  // regular expression for the volume name (empty: use contains)
    std::vector<std::string> contains;
    TEveElementList* list = nullptr;
    bool visibleByDefault = true;
    std::function<void(TEveElement*, SubDetector&)> handler;
    // optionally
    // DetectorType type;
  };
  std::vector<SubDetector> m_subDetectors;

  SubDetector* findSubDetector(TString const& s);

  // the steps of startDisplay
  // read the energy thresholds and the detector version
  void readConfiguration();
  // create the magnetic field (ideal solenoid or field map, see bField)
  void createMagneticField();
  // create the eve manager, the main window and the 3D view
  void createMainWindow();
  // create the geometry (returns false if it cannot be read)
  bool buildGeometry();
  // the G4 geometry, from the geometry cache or the geometry file
  bool buildG4Geometry();
  // the full detector (--fulldet), with the handlers of the subdetectors
  void buildFullDetector(TEveElement *world);
  // the calorimeters only (returns false if the ECAL barrel is not found)
  bool buildCalorimeters(TEveElement *world);
  // the simplified geometry, from the calculated dimensions (--nog4)
  void buildSimplifiedGeometry();
  // the rho-phi and rho-z views: the projection of the geometry and the readout
  void createRhoPhiView();
  // project the geometry and the readout in the views (once, when they are
  // first shown), showing the progress in the status bar
  void buildRhoPhiProjection();
  void buildRhoZProjection();
  void buildProjection(const char *view, TEveScene *scene, const std::function<void()> &build);
  void projectGeometryRhoPhi();
  void drawRhoPhiReadout();
  void createRhoZView();
  void projectGeometryRhoZ();
  void drawRhoZReadout();
  // write the geometry cache (when the geometry was read from the geometry file)
  void writeGeometryCache();
  // open the event file and set up the event reader
  void openEventFile();
  // legends of the cluster cell energies in the projections
  void addPaletteLegends();
  // set the 3D view as the active tab, rotate the camera and set the clipping
  void setInitialCamera();

  // decode the cellIDs with the CellIDEncoding stored in the event file
  void setCellIDEncodings();

  // final window of the sliding-window clusters: nTheta x nPhi towers of size
  // deltaTheta x deltaPhi, centred on the tower containing the cluster
  // position, with tower indices (thetaMax - theta)/deltaTheta and
  // (phi - phiMin)/deltaPhi; read from the options of the job that produced
  // the clusters, stored in the event file (nTheta = 0 if not found)
  struct SlidingWindow
  {
    int nTheta = 0;
    int nPhi = 0;
    double deltaTheta = 0.;
    double deltaPhi = 0.;
    double thetaMax = 0.;
    double phiMin = 0.;
    // radial range of the barrel calorimeters used for the towers
    double rMin = 0.;
    double rMax = 0.;
  };
  SlidingWindow swWindow;
  void setSlidingWindow();
  // theta and phi range of the final window of a cluster
  void getSlidingWindowRange(float xcl, float ycl, float zcl,
                             double &thetaLow, double &thetaHigh,
                             double &phiLow, double &phiHigh) const;
  // draw the final window of a sliding-window cluster in the projections
  void drawSlidingWindow(float xcl, float ycl, float zcl, const char *name,
                         TEveElement *rhoz, TEveElement *rhophi);

  void addSimpleElement(TEveElement* element, SubDetector&);
  void addVertex(TEveElement* element, SubDetector&);
  void addDCH(TEveElement* element, SubDetector&);
  void addWrapper(TEveElement* element, SubDetector&);
  void addECalB(TEveElement* element, SubDetector&);
  void addECalEC(TEveElement* element, SubDetector&);
  void addHCalB(TEveElement* element, SubDetector&);
  void addHCalEC(TEveElement* element, SubDetector&);
  void addMuonB(TEveElement* element, SubDetector&);
  void addMuonEC(TEveElement* element, SubDetector&);

};

#endif
